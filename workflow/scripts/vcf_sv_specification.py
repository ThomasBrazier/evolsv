import argparse
import warnings

parser = argparse.ArgumentParser()
parser.add_argument("input", help="Input VCF to parse")
parser.add_argument("output", help="Output VCF")
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
# bcf_in.header.info.add("OLDTYPE", ".", "String", "Type before DUP to INS")

bcf_out = VariantFile(output, 'w', header=bcf_in.header)
# bcf_out.header.add_meta(key="INFO", items=[("ID", "OLDTYPE"), ("Number", "1"), ("Type", "Integer"), ("Description", "Type before DUP to INS")])
# bcf_out.header.info.add("OLDTYPE", ".", "Integer", "Type before DUP to INS")

dist_to_chromosome_end = []

chrom_lengths = {name: contig.length for name, contig in bcf_in.header.contigs.items()}

output_ignored = VariantFile(output + ".ignored", "w", header=bcf_in.header)

# Iterate over the VCF
with VariantFile(output, "w", header=bcf_in.header) as out:
    # out.header.info.add("OLDTYPE", ".", "Integer", "Type before DUP to INS")
    # out.header.add_meta(key="INFO", items=[("ID", "OLDTYPE"), ("Number", "1"), ("Type", "Integer"), ("Description", "Type before DUP to INS")])

    for rec in bcf_in.fetch():
        chr = rec.chrom
        start = rec.pos

        if rec.info["SVTYPE"] == "INS":
            rec.info["OLDTYPE"] = "INS"
            # ```INFO/END``` field must be equal to ```START```
            rec.info["END"] = start + 1
            seqalt = rec.alleles[1]
            seqref = rec.ref
            # The ALT sequence contains the insertion sequence + the START nucleotide in reference genome at the beginning
            # With Debreak REF = N
            # Add the single base
            if seqref == "N" :
                single_base = sequences.fetch(chr, start - 1, start)
                seqalt = single_base + seqalt
                seqref = seqalt[0]
                rec.alleles = (seqref, seqalt)
            # ```INFO/SVLEN``` is the length of the sequence in ```ALT``` field
            if seqalt != "<INS>" :
                svlen = len(seqalt) - 1
                rec.info["SVLEN"] = svlen

        if rec.info["SVTYPE"] == "DUP":
            # HOMOGENIZE WITH INS
            # Regions with Copy Number Variation
            # The main difference is that ALT sequence of the duplication is the actual sequence in ref genome at start:end
            # END is not START + 1
            # And REF is first base of DUP sequence i.e. ALT[0]
            # ```INFO``` field must contain ```SVTYPE=INS```
            rec.info["SVTYPE"] = "INS"
            rec.info["OLDTYPE"] = "DUP"
            # ```ALT``` field must contain the sequence of the duplication
            # ```REF``` field must be single base at START i.e. fetch(start -1, start)
            # REF = ALT[0]
            end = rec.info["END"]
            # seqalt = rec.alleles[1]
            # seqref = rec.ref
            # if seqref == "N" :
            #     # Shift start by one base
            #     start = start - 1 # Include the first base before
            #     rec.pos = start
            #     # Get this first base in fasta
            #     # 0-based half-open interval
            # seqref = sequences.fetch(chr, start - 1, start)
            seqalt = sequences.fetch(chr, start - 1, end)
            seqref = seqalt[0]
            rec.alleles = (seqref, seqalt)
            # ```INFO/SVLEN``` is the length of the sequence in ```ALT``` field
            svlen = len(seqalt) - 1
            rec.info["SVLEN"] = svlen

        if rec.info["SVTYPE"] == "DEL":
            # DEL are regions where the sequence is missing in the ALT haplotype
            # START and END are coordinates of the deleted sequence in REF haplotype
            # ALT is a single base at start position (fetch(start -1, start)), REF[0]
            rec.info["OLDTYPE"] = "DEL"
            # ```INFO/END``` field must contain ```END=pos``` (with `pos` being the end position of the deleted segment)
            # ```INFO/END``` is ```START``` + (length of the sequence in ```REF``` - 1)
            seqalt = rec.alleles[1]
            seqref = rec.ref
            svlen = rec.info["SVLEN"]
            # svlen = len(seqref) - 1
            # If ALT is N, then add the a single first base (shift start)
            if seqalt == "N" :
                start = start - 1
                rec.pos = start
                # svlen = svlen + 1
            # svlen = rec.info["SVLEN"]
            # start + svlen instead of start + svlen -1 because
            end = start + abs(svlen)
            seqref = sequences.fetch(chr, start - 1, end)
            seqalt = sequences.fetch(chr, start - 1, start)
            rec.alleles = (seqref, seqalt)
            rec.info["END"] = end
            rec.info["SVLEN"] = -abs(svlen)

        if rec.info["SVTYPE"] == "INV":
            rec.info["OLDTYPE"] = "INV"
            seqref = rec.alleles[0]
            seqalt = rec.alleles[1]
            if seqalt != "<INV>":
                svlen = len(seqref)
                end = start + svlen
                rec.info["END"] = end
                rec.info["SVLEN"] = svlen
                rec.alleles = (seqref[0], "<INV>")
        
        end = rec.info["END"]
        chr_length = chrom_lengths[chr]
        dist_chr_end = chr_length - end

        dist_to_chromosome_end.append(dist_chr_end)

        # Filter out SVs starting at pos 1
        # Mostly <INV> called by CuteSV, most likely erroneous calls
        # Produces errors in SVjedi-graph
        if rec.pos < 100 or dist_chr_end < 100:
            print("SV starting or ending at less than 100bp of chromosome end")
            output_ignored.write(rec)
        else:
            out.write(rec)


with open(output + '.dist_to_chromosome_end.txt', 'w') as f:
    for line in dist_to_chromosome_end:
        f.write(f"{line}\n")
