import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input", help="Input VCF to parse")
parser.add_argument("output", help="Output VCF with sequences added to ALT")
args = parser.parse_args()

input = args.input
output = args.output

# input="data/GCA_947369205.1/GCA_947369205.1_minimap2_debreak_normalize.vcf"
# output="data/GCA_947369205.1/test.vcf"

# Import VCF
from pysam import VariantFile

bcf_in = VariantFile(input)  # auto-detect input format
bcf_out = VariantFile(output, 'w', header=bcf_in.header)

# Iterate over the VCF to add symbolic type to REF/ALT field
with VariantFile(output, "w", header=bcf_in.header) as out:

    for rec in bcf_in.fetch():
        if "INS" in rec.info["OLDTYPE"]:
            rec.alleles[1] = "<INS>"
        
        if "DEL" in rec.info["OLDTYPE"]:
            rec.alleles[0] = "<DEL>"
        
        if "DUP" in rec.info["OLDTYPE"]:
            rec.alleles[1] = "<DUP>"
        
        if "INV" in rec.info["OLDTYPE"]:
            rec.alleles[1] = "<INV>"

        if "TRA" in rec.info["OLDTYPE"]:
            rec.alleles[1] = "<TRA>"

        out.write(rec)

        
        



