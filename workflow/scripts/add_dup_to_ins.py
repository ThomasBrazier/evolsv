import argparse
import warnings

parser = argparse.ArgumentParser()
parser.add_argument("input", help="Input VCF to parse")
parser.add_argument("output", help="Output VCF with sequences added to ALT")
parser.add_argument("fasta", help="Genome FASTA file")
args = parser.parse_args()

input = args.input
output = args.output
fasta = args.fasta

# input="../../data/GCA_947369205.1/GCA_947369205.1_minimap2_debreak_normalize.vcf"
# output="../../data/GCA_947369205.1/test.vcf"
# fasta="../../data/GCA_947369205.1/GCA_947369205.1.fna"

# Import FASTA
from pysam import FastaFile
 # read FASTA file
sequences = FastaFile(fasta)

# Import VCF
from pysam import VariantFile

bcf_in = VariantFile(input)  # auto-detect input format
bcf_out = VariantFile(output, 'w', header=bcf_in.header)

# for rec in bcf_in.fetch():
#     if "<DUP>" in rec.alleles:
#         print(rec.alts)
#         # Get sequence from FASTA
#         chr = rec.chrom
#         start = rec.pos
#         end = start + rec.rlen
#         rec.info.__setitem__('END', end)
#         # print(chr)
#         # print(start)
#         # print(end)
#         tmp_sequence = sequences.fetch(chr, start, end)
#         # print(tmp_sequence)
#         # Add sequence to ALT field
#         rec.alleles=('N', tmp_sequence)
#     # Save new line in new vcf
#     bcf_out.write(rec)

# Iterate over the VCF to add sequence to ALT field
with VariantFile(output, "w", header=bcf_in.header) as out:
    # out.header.info.add("END", ".", "Integer", "End position of the variant described in this record")
    # out.header.add_meta(key="INFO", items=[("ID", "END"), ("Number", "1"), ("Type", "Integer"), ("Description", "End coordinate")])

    for rec in bcf_in.fetch():
        # Coerce INFO/END, required for SVjedi graph
        chr = rec.chrom
        start = rec.pos
        # END=POS + LEN(REF) - 1
        end = start + rec.rlen - 1

        # if end >= len(sequences[chr]):
        #     warnings.warn("Out of index END " + end + ">=" + len(sequences[chr]))

        # rec.stop = end
        
        if "<DUP>" in rec.alleles:
            # print(rec.alts)
            # Get sequence from FASTA
            tmp_sequence = sequences.fetch(chr, start, end)
            # print(tmp_sequence)
            # Add sequence to ALT field
            rec.alleles=('N', tmp_sequence)
            # Save new line in new vcf
            # Remove DUP with POS < 1 
            # errors in SVjedi-graph when too close to chromosome start (!!!)
            if start > 10:
                out.write(rec)
        else:
            out.write(rec)
        
        



