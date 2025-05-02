import argparse
import warnings

parser = argparse.ArgumentParser()
parser.add_argument("input", help="Input VCF to parse")
parser.add_argument("output", help="Output VCF with SVLEN")
parser.add_argument("fasta", help="Genome FASTA file")

args = parser.parse_args()

input = args.input
output = args.output
fasta = args.fasta

# Import FASTA
from pysam import FastaFile
 # read FASTA file
sequences = FastaFile(fasta)


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
        
        if rec.info["SVTYPE"] == "DEL":
            # shift START
            start = start - 1
            
            svlen = rec.info["SVLEN"]
            end = start + svlen

            # Get new sequence from FASTA
            # Add the first base to beginning of seq and ALT
            tmp_sequence = sequences.fetch(chr, start - 1, end)
            single_base = sequences.fetch(chr, start - 1, start)

            rec.alleles=(tmp_sequence, single_base)

            # end = rec.info["AVG_END"]
            # svlen = round(svlen, 0)
            # end = round(end, 0)
            rec.info["END"] = end
            rec.info["SVLEN"] = -abs(svlen)
            rec.pos = start
            out.write(rec)
        else:
            out.write(rec)
        
        



