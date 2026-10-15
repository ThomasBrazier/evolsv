"""Re-attach per-tool FORMAT fields (DP, AD, PL) to the Jasmine-merged VCF.

Jasmine (--output_genotypes) keeps only GT in the per-callset sample columns, although
each filtered callset VCF it merged was genotyped by SVjedi-graph and carries
GT:DP:AD:PL. This script puts those fields back, so that the callset columns of the
final VCF are as informative as the SVjedi-graph column on the merged calls.

For each merged call and each callset whose SUPP_VEC bit is 1, the call is matched to
the record of that callset it was built from:
  1. candidates are the callset's records on the same CHROM, with the same SVTYPE and
     POS within --max-dist of the merged call (Jasmine's jasmine_max_dist);
  2. among them, one whose ID is listed in the merged call's INFO/IDLIST is preferred,
     otherwise the nearest by POS. IDs alone are not trusted: the same caller run on two
     aligners can produce the same ID, so the position check always applies.
Jasmine does not merge two records of one callset (--allow_intrasample is off), so at
most one record per callset is expected.

Callsets that do not support a call are left untouched (missing data, not absence of
variation). A supporting callset with no matching record is also left untouched and
counted; the script fails if more than --max-unmatched-fraction of the supporting
callset/call pairs cannot be matched, because that means the matching assumptions above
do not hold for this data.

Coordinates: VCF POS, 1-based, on both sides. Only the standard library is used.
"""

import argparse
import bisect
import sys

KEYS = ("GT", "DP", "AD", "PL")
MISSING = {"GT": "./."}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("merged", help="Jasmine merged VCF")
    parser.add_argument("output", help="Output VCF")
    parser.add_argument(
        "--callset",
        nargs=2,
        action="append",
        metavar=("NAME", "VCF"),
        required=True,
        help="Callset name and its genotyped VCF, in SUPP_VEC / sample-column order",
    )
    parser.add_argument("--max-dist", type=int, required=True)
    parser.add_argument("--max-unmatched-fraction", type=float, default=0.05)
    return parser.parse_args()


def info_value(info, key):
    for field in info.split(";"):
        if field.startswith(key + "="):
            return field[len(key) + 1 :]
    return None


def read_callset(path):
    """Return (header FORMAT lines, {chrom: sorted [(pos, id, svtype, fields)]})."""
    formats = {}
    records = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("##FORMAT="):
                formats[line.split("ID=", 1)[1].split(",", 1)[0]] = line.rstrip("\n")
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                continue
            fields = dict(zip(cols[8].split(":"), cols[9].split(":")))
            records.setdefault(cols[0], []).append(
                (int(cols[1]), cols[2], info_value(cols[7], "SVTYPE"), fields)
            )
    for chrom_records in records.values():
        chrom_records.sort(key=lambda r: r[0])
    return formats, records


def find_match(records, chrom, pos, svtype, idlist, max_dist):
    chrom_records = records.get(chrom, [])
    positions = [r[0] for r in chrom_records]
    lo = bisect.bisect_left(positions, pos - max_dist)
    hi = bisect.bisect_right(positions, pos + max_dist)
    candidates = [r for r in chrom_records[lo:hi] if r[2] == svtype]
    if not candidates:
        return None
    listed = [r for r in candidates if r[1] in idlist]
    return min(listed or candidates, key=lambda r: abs(r[0] - pos))


def main():
    args = parse_args()
    names = [name for name, _ in args.callset]
    callsets = [read_callset(path) for _, path in args.callset]

    # FORMAT header definitions are taken from the callset VCFs, not invented here.
    format_lines = {}
    for formats, _ in callsets:
        for key in KEYS:
            if key in formats:
                format_lines.setdefault(key, formats[key])
    missing_defs = [k for k in KEYS if k not in format_lines and k != "GT"]
    if missing_defs:
        sys.exit(f"No ##FORMAT definition for {missing_defs} in any callset VCF.")

    supported = 0
    unmatched = 0
    per_callset_unmatched = dict.fromkeys(names, 0)
    reused = 0
    used = set()

    with open(args.merged) as merged, open(args.output, "w") as out:
        header_done = False
        seen_formats = []
        for line in merged:
            if line.startswith("##"):
                if line.startswith("##FORMAT="):
                    seen_formats.append(line.split("ID=", 1)[1].split(",", 1)[0])
                out.write(line)
                continue
            if line.startswith("#CHROM"):
                existing = set(seen_formats)
                for key in KEYS:
                    if key not in existing and key in format_lines:
                        out.write(format_lines[key] + "\n")
                n_samples = len(line.rstrip("\n").split("\t")) - 9
                if n_samples != len(names):
                    sys.exit(
                        f"Merged VCF has {n_samples} sample columns but "
                        f"{len(names)} callsets were given."
                    )
                out.write(line)
                header_done = True
                continue
            if not header_done:
                sys.exit("Merged VCF has no #CHROM header line.")

            cols = line.rstrip("\n").split("\t")
            chrom, pos, info = cols[0], int(cols[1]), cols[7]
            svtype = info_value(info, "SVTYPE")
            supp_vec = info_value(info, "SUPP_VEC")
            if supp_vec is None or len(supp_vec) != len(names):
                sys.exit(
                    f"{chrom}:{pos}: SUPP_VEC={supp_vec} does not have one bit per "
                    f"callset ({len(names)})."
                )
            idlist = set((info_value(info, "IDLIST") or "").split(","))

            old_keys = cols[8].split(":")
            extra = [k for k in old_keys if k not in KEYS]
            new_keys = list(KEYS) + extra
            new_samples = []
            for i, sample in enumerate(cols[9:]):
                values = dict(zip(old_keys, sample.split(":")))
                if supp_vec[i] == "1":
                    supported += 1
                    match = find_match(
                        callsets[i][1], chrom, pos, svtype, idlist, args.max_dist
                    )
                    if match is None:
                        unmatched += 1
                        per_callset_unmatched[names[i]] += 1
                    else:
                        if (i, match[1], match[0]) in used:
                            reused += 1
                        used.add((i, match[1], match[0]))
                        for key in KEYS:
                            if key in match[3]:
                                values[key] = match[3][key]
                new_samples.append(
                    ":".join(values.get(k, MISSING.get(k, ".")) for k in new_keys)
                )
            out.write(
                "\t".join(cols[:8] + [":".join(new_keys)] + new_samples) + "\n"
            )

    fraction = unmatched / supported if supported else 0.0
    sys.stderr.write(
        f"supporting callset/call pairs: {supported}; unmatched: {unmatched} "
        f"({fraction:.4f}); source records matched more than once: {reused}\n"
    )
    for name, count in per_callset_unmatched.items():
        sys.stderr.write(f"  unmatched {name}: {count}\n")
    if fraction > args.max_unmatched_fraction:
        sys.exit(
            f"Unmatched fraction {fraction:.4f} exceeds "
            f"--max-unmatched-fraction {args.max_unmatched_fraction}."
        )


if __name__ == "__main__":
    main()
