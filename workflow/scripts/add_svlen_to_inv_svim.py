import argparse
import warnings

parser = argparse.ArgumentParser()
parser.add_argument("input", help="Input VCF to parse")
parser.add_argument("output", help="Output VCF with SVLEN")
args = parser.parse_args()

input = args.input
output = args.output

# Import VCF
from pysam import VariantFile

bcf_in = VariantFile(input)  # auto-detect input format
bcf_out = VariantFile(output, 'w', header=bcf_in.header)

# Iterate over the VCF to add sequence to ALT field
with VariantFile(output, "w", header=bcf_in.header) as out:
    # out.header.info.add("END", ".", "Integer", "End position of the variant described in this record")
    # out.header.add_meta(key="INFO", items=[("ID", "END"), ("Number", "1"), ("Type", "Integer"), ("Description", "End coordinate")])

    for rec in bcf_in.fetch():
        # Coerce INFO/END, required for SVjedi graph
        chr = rec.chrom
        start = rec.pos
        # END=POS + LEN(REF) - 1
        # end = start + rec.rlen - 1
        end = rec.info["SVLEN"]
        
        if "<INV>" in rec.alleles:
            svlen = abs(end - start)\
            rec.info["SVLEN"]=svlen
            out.write(rec)
        else:
            out.write(rec)
        
        



