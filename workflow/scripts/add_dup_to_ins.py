import argparse
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

# Iterate over the VCF to add sequence to ALT field
for rec in bcf_in.fetch():
    if "<DUP>" in rec.alts:
        print(rec.alts)
        # Get sequence from FASTA
        chr = rec.chrom
        start = rec.start
        end = rec.stop
        print(chr)
        print(start)
        print(end)
        tmp_sequence = sequences.fetch(chr, start, end)
        print(tmp_sequence)
        # Add sequence to ALT field
        rec.alts=(tmp_sequence,)
    # Save new line in new vcf
    bcf_out.write(rec)

