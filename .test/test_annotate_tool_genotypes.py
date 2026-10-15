"""Known-truth tests for workflow/scripts/annotate_tool_genotypes.py.

Synthetic VCFs only: they check the matching logic (SUPP_VEC gating, ID collisions across
callsets, position/type checks, failure on unmatched calls), not Jasmine's real output
format, which has not been inspected in this repository.
"""

import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "annotate_tool_genotypes.py"

FORMATS = (
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depth">\n'
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">\n'
)


def callset_vcf(path, records):
    lines = ["##fileformat=VCFv4.2\n", FORMATS, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n"]
    for pos, svtype, id_, sample in records:
        lines.append(f"chr1\t{pos}\t{id_}\tN\t<{svtype}>\t.\tPASS\tSVTYPE={svtype}\tGT:DP:AD:PL\t{sample}\n")
    path.write_text("".join(lines))


def merged_vcf(path, records):
    lines = [
        "##fileformat=VCFv4.2\n",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n",
    ]
    for pos, svtype, supp_vec, idlist, gts in records:
        info = f"SVTYPE={svtype};SUPP_VEC={supp_vec};IDLIST={idlist}"
        lines.append(f"chr1\t{pos}\t1\tN\t<{svtype}>\t.\tPASS\t{info}\tGT\t{gts[0]}\t{gts[1]}\n")
    path.write_text("".join(lines))


def run(tmp_path, merged, a, b, extra=()):
    out = tmp_path / "out.vcf"
    return out, subprocess.run(
        [
            sys.executable, "-I", str(SCRIPT), str(merged), str(out),
            "--callset", "a", str(a), "--callset", "b", str(b),
            "--max-dist", "100", *extra,
        ],
        capture_output=True, text=True, check=False,
    )


def sample_columns(out):
    rows = [line.split("\t") for line in out.read_text().splitlines() if not line.startswith("#")]
    return [(r[8], r[9], r[10]) for r in rows]


def test_fields_copied_only_for_supporting_callsets(tmp_path):
    a, b, merged = tmp_path / "a.vcf", tmp_path / "b.vcf", tmp_path / "m.vcf"
    callset_vcf(a, [(1000, "DEL", "x.1", "0/1:20:10,10:50,0,60")])
    callset_vcf(b, [(1000, "DEL", "x.1", "1/1:30:2,28:90,20,0")])
    merged_vcf(merged, [(1005, "DEL", "10", "x.1,x.1", ("0/1", "./."))])
    out, result = run(tmp_path, merged, a, b)
    assert result.returncode == 0, result.stderr
    (fmt, col_a, col_b), = sample_columns(out)
    assert fmt == "GT:DP:AD:PL"
    assert col_a == "0/1:20:10,10:50,0,60"
    # b is not in SUPP_VEC: left as Jasmine wrote it, never filled from a nearby record.
    assert col_b == "./.:.:.:."


def test_id_collision_across_callsets_resolved_by_position(tmp_path):
    # Two callsets reuse the ID "x.1" for different loci; the one in range must win.
    a, b, merged = tmp_path / "a.vcf", tmp_path / "b.vcf", tmp_path / "m.vcf"
    callset_vcf(a, [(1000, "DEL", "x.1", "0/1:20:10,10:50,0,60"), (5000, "DEL", "x.2", "1/1:9:0,9:99,9,0")])
    callset_vcf(b, [(5000, "DEL", "x.1", "0/0:7:7,0:0,9,99")])
    merged_vcf(merged, [(5000, "DEL", "01", "x.1", ("./.", "0/0"))])
    out, result = run(tmp_path, merged, a, b)
    assert result.returncode == 0, result.stderr
    (_, _, col_b), = sample_columns(out)
    assert col_b == "0/0:7:7,0:0,9,99"


def test_type_mismatch_is_unmatched_and_fails(tmp_path):
    a, b, merged = tmp_path / "a.vcf", tmp_path / "b.vcf", tmp_path / "m.vcf"
    callset_vcf(a, [(1000, "INS", "x.1", "0/1:20:10,10:50,0,60")])
    callset_vcf(b, [])
    merged_vcf(merged, [(1000, "DEL", "10", "x.1", ("0/1", "./."))])
    _, result = run(tmp_path, merged, a, b)
    assert result.returncode != 0
    assert "Unmatched fraction" in result.stderr


@pytest.mark.parametrize("supp_vec", ["1", "111"])
def test_supp_vec_length_must_match_callsets(tmp_path, supp_vec):
    a, b, merged = tmp_path / "a.vcf", tmp_path / "b.vcf", tmp_path / "m.vcf"
    callset_vcf(a, [])
    callset_vcf(b, [])
    merged_vcf(merged, [(1000, "DEL", supp_vec, "x.1", ("0/1", "0/1"))])
    _, result = run(tmp_path, merged, a, b)
    assert result.returncode != 0
    assert "SUPP_VEC" in result.stderr
