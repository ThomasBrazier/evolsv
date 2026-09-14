"""Unit tests for workflow/scripts/check_reference_names.py.

Scope: the check that the sequence report and the reference FASTA name the same
contigs. A mismatch there empties the final VCF without any error (see the script
docstring), so the pass, warn and fail branches are each tested on small files.
These tests run the script as a subprocess, as rule check_reference_names does.
"""

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "workflow" / "scripts" / "check_reference_names.py"


def molecule(name, accession, length, assembly="GCA_000000001.1", **extra):
    row = {
        "assemblyAccession": assembly,
        "chrName": name,
        "role": "assembled-molecule",
        "length": length,
        "genbankAccession": accession,
    }
    row.update(extra)
    return row


def run_check(tmp_path, rows, fai, raw_lines=None):
    report = tmp_path / "sequence_report.jsonl"
    lines = raw_lines if raw_lines is not None else [json.dumps(row) for row in rows]
    report.write_text("\n".join(lines) + "\n")
    fai_path = tmp_path / "ref.fna.fai"
    fai_path.write_text("".join(f"{name}\t{length}\t0\t60\t61\n" for name, length in fai))
    out = tmp_path / "check.txt"
    result = subprocess.run(
        [
            sys.executable, str(SCRIPT),
            "--sequence-report", str(report), "--fai", str(fai_path), "--report", str(out),
        ],
        capture_output=True, text=True, check=False,
    )
    return result, out


def test_matching_report_passes(tmp_path):
    rows = [molecule("1", "OX1.1", 100), molecule("Z", "OX2.1", 50)]
    result, out = run_check(tmp_path, rows, [("OX1.1", 100), ("OX2.1", 50)])
    assert result.returncode == 0, result.stderr
    assert "assembled molecules found in the FASTA: 2" in out.read_text()
    assert "WARNING" not in out.read_text()


def test_no_name_in_fasta_fails(tmp_path):
    """GenBank report against a RefSeq-named FASTA: the final VCF would be empty."""
    rows = [molecule("1", "OX1.1", 100)]
    result, out = run_check(tmp_path, rows, [("NC_1.1", 100)])
    assert result.returncode == 1
    assert "none of the 1 assembled molecules" in result.stderr
    assert not out.exists()


def test_length_mismatch_fails(tmp_path):
    rows = [molecule("1", "OX1.1", 100)]
    result, _ = run_check(tmp_path, rows, [("OX1.1", 99)])
    assert result.returncode == 1
    assert "different assembly versions" in result.stderr


def test_partial_mismatch_warns(tmp_path):
    rows = [molecule("1", "OX1.1", 100), molecule("2", "OX2.1", 80)]
    result, out = run_check(tmp_path, rows, [("OX1.1", 100)])
    assert result.returncode == 0, result.stderr
    assert "WARNING: 1 of 2 assembled molecule(s) are not in the FASTA index: OX2.1" in (
        out.read_text()
    )


def test_refseq_assembly_uses_refseq_accession(tmp_path):
    """Same column choice as autosomes_sexchromosomes.R: GCF_ -> refseqAccession."""
    rows = [molecule("1", "OX1.1", 100, assembly="GCF_000000001.1", refseqAccession="NC_1.1")]
    result, _ = run_check(tmp_path, rows, [("NC_1.1", 100)])
    assert result.returncode == 0, result.stderr

    result, _ = run_check(tmp_path, rows, [("OX1.1", 100)])
    assert result.returncode == 1


def test_two_assembly_accessions_fail(tmp_path):
    rows = [molecule("1", "OX1.1", 100), molecule("2", "OX2.1", 80, assembly="GCA_2.1")]
    result, _ = run_check(tmp_path, rows, [("OX1.1", 100), ("OX2.1", 80)])
    assert result.returncode == 1
    assert "mixes 2 assembly accessions" in result.stderr


def test_scaffold_level_assembly_warns(tmp_path):
    rows = [molecule("Un", "OX1.1", 100, role="unplaced-scaffold")]
    result, out = run_check(tmp_path, rows, [("OX1.1", 100)])
    assert result.returncode == 0, result.stderr
    assert "scaffold-level assembly" in out.read_text()


def test_missing_fields_and_bad_json_fail(tmp_path):
    rows = [{"assemblyAccession": "GCA_1.1", "role": "assembled-molecule", "length": 5}]
    result, _ = run_check(tmp_path, rows, [("OX1.1", 5)])
    assert result.returncode == 1
    assert "without 'chrName' or 'genbankAccession'" in result.stderr

    result, _ = run_check(tmp_path, [], [("OX1.1", 5)], raw_lines=["{not json"])
    assert result.returncode == 1
    assert "line 1 is not valid JSON" in result.stderr
