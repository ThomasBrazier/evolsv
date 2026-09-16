"""Check that the sequence report and the reference FASTA describe the same assembly.

Used by rule check_reference_names, before rule autosomes_sexchromosomes. That rule
(workflow/scripts/autosomes_sexchromosomes.R) builds the autosome and sex-chromosome
BED files from the sequence report, and rule final_vcf keeps only the variants on
those contigs. If the contig names of the report are not the names in the FASTA, the
BED files select nothing and the final VCF is empty, with no error.

The contig name is chosen as in autosomes_sexchromosomes.R: `genbankAccession` when
`assemblyAccession` contains "GCA_", else `refseqAccession`. Only the rows with
`role == "assembled-molecule"` go into the BED files.

Errors (exit 1): unreadable report, missing field, more than one assemblyAccession,
no assembled molecule found in the FASTA index, or a length that differs from the
FASTA index. Warnings (report only): some assembled molecules absent from the FASTA
index, or no assembled molecule at all (scaffold-level assembly).

Every problem found is reported, not just the first.
"""

import argparse
import json
import sys

ASSEMBLED = "assembled-molecule"


def read_fai(fai_path):
    """Parse a samtools faidx index into an ordered {contig: length} mapping.

    Same parser as workflow/scripts/check_bam_reference.py, which cannot be imported
    here because it needs pysam.
    """
    lengths = {}
    with open(fai_path) as fai:
        for line in fai:
            if not line.strip():
                continue
            fields = line.split("\t")
            lengths[fields[0]] = int(fields[1])
    return lengths


def read_report(report_path):
    """Return (rows, errors). One JSON object per non-empty line (JSON Lines)."""
    rows = []
    errors = []
    with open(report_path, encoding="utf-8") as report:
        for number, line in enumerate(report, start=1):
            if not line.strip():
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError as error:
                errors.append(f"line {number} is not valid JSON ({error.msg}).")
                continue
            if not isinstance(row, dict):
                errors.append(f"line {number} is not a JSON object.")
                continue
            rows.append((number, row))
    if not rows and not errors:
        errors.append("the sequence report contains no sequence.")
    return rows, errors


def accession_field(rows):
    """Return (field name, assembly accession, errors), as autosomes_sexchromosomes.R."""
    missing = [number for number, row in rows if not row.get("assemblyAccession")]
    if missing:
        return None, None, [f"assemblyAccession is missing on line(s) {_first(missing)}."]
    accessions = sorted({str(row["assemblyAccession"]) for _, row in rows})
    if len(accessions) > 1:
        return None, None, [
            "the report mixes {} assembly accessions ({}). It must describe one "
            "assembly.".format(len(accessions), ", ".join(accessions))
        ]
    field = "genbankAccession" if "GCA_" in accessions[0] else "refseqAccession"
    return field, accessions[0], []


def check(rows, reference_lengths):
    """Return (errors, warnings, assembled contig names) for the parsed report."""
    field, _, errors = accession_field(rows)
    if errors:
        return errors, [], []

    warnings = []
    assembled = []
    for number, row in rows:
        if "role" not in row or "length" not in row:
            errors.append(f"line {number} has no 'role' or no 'length' field.")
            continue
        if row["role"] != ASSEMBLED:
            continue
        if not row.get("chrName") or not row.get(field):
            errors.append(
                f"line {number} is an assembled molecule without 'chrName' or '{field}'."
            )
            continue
        try:
            length = int(row["length"])
        except (TypeError, ValueError):
            errors.append(f"line {number} has a non-integer length {row['length']!r}.")
            continue
        assembled.append((row[field], row["chrName"], length))

    if errors:
        return errors, warnings, []

    if not assembled:
        warnings.append(
            f"no sequence has role '{ASSEMBLED}' (scaffold-level assembly?). The "
            "autosome and sex-chromosome BED files will be empty, so the final VCF "
            "files will contain no variant."
        )
        return errors, warnings, []

    found = [contig for contig in assembled if contig[0] in reference_lengths]
    missing = [contig for contig in assembled if contig[0] not in reference_lengths]
    if not found:
        expected = _first([name for name, _, _ in assembled])
        in_fasta = _first(list(reference_lengths))
        errors.append(
            f"none of the {len(assembled)} assembled molecules of the report is in the "
            f"FASTA index (names from '{field}': {expected}; FASTA has: {in_fasta}). The "
            "report and the FASTA use different contig names (GenBank vs RefSeq vs UCSC) "
            "or describe different assemblies."
        )
    elif missing:
        names = _first([name for name, _, _ in missing])
        warnings.append(
            f"{len(missing)} of {len(assembled)} assembled molecule(s) are not in the "
            f"FASTA index: {names}. Variants on these molecules cannot be kept."
        )

    mismatched = [
        f"{name} (report {length} bp, FASTA {reference_lengths[name]} bp)"
        for name, _, length in found
        if reference_lengths[name] != length
    ]
    if mismatched:
        errors.append(
            f"{len(mismatched)} contig(s) have a different length in the report and in "
            f"the FASTA: {_first(mismatched, separator='; ')}. These are different "
            "assembly versions."
        )

    return errors, warnings, [name for name, _, _ in found]


def _first(items, limit=10, separator=", "):
    items = [str(item) for item in items]
    return separator.join(items[:limit]) + (separator + "..." if len(items) > limit else "")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequence-report", required=True, help="sequence_report.jsonl")
    parser.add_argument("--fai", required=True, help="samtools faidx index of the reference")
    parser.add_argument("--report", required=True, help="path of the check report to write")
    args = parser.parse_args()

    rows, errors = read_report(args.sequence_report)
    reference_lengths = read_fai(args.fai)
    warnings = []
    found = []
    if not errors:
        errors, warnings, found = check(rows, reference_lengths)

    if errors:
        sys.stderr.write(
            f"\n{len(errors)} check(s) failed for the sequence report "
            f"'{args.sequence_report}' against '{args.fai}':\n\n"
        )
        for error in errors:
            sys.stderr.write(f"  * {error}\n\n")
        sys.exit(1)

    with open(args.report, "w") as report:
        report.write("Sequence report checked against the reference FASTA\n")
        report.write("=" * 52 + "\n\n")
        report.write(f"sequence report:   {args.sequence_report}\n")
        report.write(f"reference index:   {args.fai}\n")
        report.write(f"report sequences:  {len(rows)}\n")
        report.write(f"FASTA contigs:     {len(reference_lengths)}\n")
        report.write(f"assembled molecules found in the FASTA: {len(found)}\n\n")
        report.writelines(f"WARNING: {warning}\n" for warning in warnings)
    for warning in warnings:
        sys.stderr.write(f"WARNING: {warning}\n")


if __name__ == "__main__":
    main()
