"""Validate a user-supplied, pre-aligned BAM before it enters the workflow.

Used by rule stage_bam (start_from_bam mode). The checks below all guard against
failures that would otherwise surface much later, and far less legibly:

* an unsorted BAM is rejected by the SV callers;
* a missing BAI index makes every caller and mosdepth fail;
* contigs that disagree with the reference produce silent nonsense in
  `bcftools reheader --fai` (workflow/rules/merging.smk, workflow/rules/genotyping.smk)
  and raise KeyError inside workflow/scripts/vcf_sv_specification.py;
* an @RG SM tag that differs from sample_name makes samplot emit empty plots,
  because rule samplot_plot passes `--sample_ids {sample_id}`.

Every problem found is reported, not just the first, so that a mis-specified BAM
can be fixed in one pass.
"""

import argparse
import sys

import pysam


def read_header(bam_path):
    """Return the BAM header as a plain dict, and whether a BAI index is present.

    pysam exposes the header as a dict in old versions and as an AlignmentHeader
    object from 0.14 onwards; both are handled here.
    """
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        header = bam.header
        header = header.to_dict() if hasattr(header, "to_dict") else dict(header)
        has_index = bam.has_index()
    return header, has_index


def read_fai(fai_path):
    """Parse a samtools faidx index into an ordered {contig: length} mapping."""
    lengths = {}
    with open(fai_path) as fai:
        for line in fai:
            if not line.strip():
                continue
            fields = line.split("\t")
            lengths[fields[0]] = int(fields[1])
    return lengths


def check_sorted(header):
    sort_order = header.get("HD", {}).get("SO")
    if sort_order != "coordinate":
        return [
            "BAM is not coordinate-sorted (@HD SO:{}). "
            "Sort it with `samtools sort` before running in start_from_bam mode.".format(
                sort_order or "absent"
            )
        ]
    return []


def check_index(has_index):
    if not has_index:
        return ["BAM index could not be opened. Create it with `samtools index`."]
    return []


def check_contigs(header, reference_lengths):
    """Every BAM contig must exist in the reference with an identical length."""
    problems = []
    bam_contigs = [(sq["SN"], int(sq["LN"])) for sq in header.get("SQ", [])]

    if not bam_contigs:
        return ["BAM header declares no @SQ contigs."]

    missing = [name for name, _ in bam_contigs if name not in reference_lengths]
    if missing:
        problems.append(
            "{} contig(s) in the BAM are absent from the reference: {}{}. "
            "The BAM was most likely aligned against a differently named copy of the "
            "assembly (RefSeq/UCSC vs GenBank naming). Set `reference_fasta` in "
            "config.yaml to the exact FASTA used for the alignment.".format(
                len(missing),
                ", ".join(missing[:10]),
                ", ..." if len(missing) > 10 else "",
            )
        )

    mismatched = [
        "{} (BAM {} bp, reference {} bp)".format(
            name, length, reference_lengths[name]
        )
        for name, length in bam_contigs
        if name in reference_lengths and reference_lengths[name] != length
    ]
    if mismatched:
        problems.append(
            "{} contig(s) have a different length in the BAM and in the reference: "
            "{}{}. These are different assemblies, not just different names.".format(
                len(mismatched),
                "; ".join(mismatched[:10]),
                "; ..." if len(mismatched) > 10 else "",
            )
        )

    return problems


def check_read_group(header, sample_id):
    """All @RG SM tags must match the sample_name given in the sample sheet."""
    read_groups = header.get("RG", [])
    if not read_groups:
        return [
            "BAM header declares no @RG line, so no sample name is attached to the "
            "reads. Add one with "
            "`samtools addreplacerg -r 'ID:{0}' -r 'SM:{0}'`.".format(sample_id)
        ]

    wrong = sorted({rg.get("SM", "<absent>") for rg in read_groups} - {sample_id})
    if wrong:
        return [
            "@RG SM tag(s) {} do not match the sample_name '{}' from the sample sheet. "
            "rule samplot_plot selects reads by sample id and would produce empty "
            "plots. Fix the tag with `samtools addreplacerg`, or change sample_name "
            "to match the BAM.".format(
                ", ".join("'{}'".format(s) for s in wrong), sample_id
            )
        ]
    return []


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True, help="user-supplied pre-aligned BAM")
    parser.add_argument("--fai", required=True, help="samtools faidx index of the reference")
    parser.add_argument("--sample-id", required=True, help="sample_name from the sample sheet")
    parser.add_argument("--aligner", required=True, help="aligner the BAM was produced with")
    parser.add_argument("--report", required=True, help="path of the provenance report to write")
    args = parser.parse_args()

    header, has_index = read_header(args.bam)
    reference_lengths = read_fai(args.fai)

    problems = (
        check_sorted(header)
        + check_index(has_index)
        + check_contigs(header, reference_lengths)
        + check_read_group(header, args.sample_id)
    )

    if problems:
        sys.stderr.write(
            "\n{} check(s) failed for the {} BAM '{}':\n\n".format(
                len(problems), args.aligner, args.bam
            )
        )
        for problem in problems:
            sys.stderr.write("  * {}\n\n".format(problem))
        sys.exit(1)

    bam_contigs = header.get("SQ", [])
    # Contigs present in the reference but not in the BAM are legitimate (the BAM may
    # have been aligned against a subset), but worth recording.
    absent_from_bam = len(reference_lengths) - len(bam_contigs)

    with open(args.report, "w") as report:
        report.write("Pre-aligned BAM accepted by EvolSV (start_from_bam mode)\n")
        report.write("=" * 56 + "\n\n")
        report.write("aligner:          {}\n".format(args.aligner))
        report.write("source BAM:       {}\n".format(args.bam))
        report.write("reference index:  {}\n".format(args.fai))
        report.write("sample id:        {}\n".format(args.sample_id))
        report.write("sort order:       coordinate\n")
        report.write("BAM contigs:      {}\n".format(len(bam_contigs)))
        report.write("reference contigs: {}\n".format(len(reference_lengths)))
        report.write(
            "reference contigs absent from the BAM: {}\n\n".format(absent_from_bam)
        )
        report.write(
            "CAVEAT: in start_from_bam mode the chopper read filters\n"
            "(chopper_quality, chopper_minlength, chopper_maxlength, chopper_headcrop,\n"
            "chopper_tailcrop) are NOT applied. The alignment is used as supplied, and\n"
            "the reads passed to SVJedi-graph for genotyping are the ones declared in\n"
            "the sample sheet, unfiltered. Read-level QC (FastQC, NanoPlot) is skipped;\n"
            "alignment QC is still produced in mapping_QC/ and callability/.\n"
        )


if __name__ == "__main__":
    main()
