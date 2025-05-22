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
bcf_in.header.info.add("OLDTYPE", ".", "String", "Type before DUP to INS")

bcf_out = VariantFile(output, 'w', header=bcf_in.header)
# bcf_out.header.add_meta(key="INFO", items=[("ID", "OLDTYPE"), ("Number", "1"), ("Type", "Integer"), ("Description", "Type before DUP to INS")])
# bcf_out.header.info.add("OLDTYPE", ".", "Integer", "Type before DUP to INS")

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
            # With Debreak REF = N
            if seqref == "N" :
                seqref = seqalt[0]
                rec.alleles = (seqref, seqalt)
            # ```INFO/SVLEN``` is the length of the sequence in ```ALT``` field
            svlen = len(seqalt)
            rec.info["SVLEN"] = svlen

        if rec.info["SVTYPE"] == "DUP":
            # HOMOGENIZE WITH INS
            # ```INFO``` field must contain ```SVTYPE=INS```
            rec.info["SVTYPE"] = "INS"
            rec.info["OLDTYPE"] = "DUP"
            # ```ALT``` field must contain the sequence of the duplication
            # ```REF``` field must be single base at START
            end = rec.info["END"]
            seqalt = rec.alleles[1]
            seqref = rec.ref
            if seqref == "N" :
                # Shift start by one base
                start = start - 1 # Include the first base before
                rec.pos = start
                # Get this first base in fasta
                # 0-based half-open interval
            seqref = sequences.fetch(chr, start - 1, start)
            seqalt = sequences.fetch(chr, start - 1, end)
            rec.alleles = (seqref, seqalt)
            # ```INFO/SVLEN``` is the length of the sequence in ```ALT``` field
            svlen = len(seqalt) - 1
            rec.info["SVLEN"] = svlen

        if rec.info["SVTYPE"] == "DEL":
            rec.info["OLDTYPE"] = "DEL"
            # ```INFO/END``` field must contain ```END=pos``` (with `pos` being the end position of the deleted segment)
            # ```INFO/END``` is ```START``` + (length of the sequence in ```REF``` - 1)
            seqalt = rec.alleles[1]
            seqref = rec.ref
            svlen = rec.info["SVLEN"]
            end = start + abs(svlen)
            if seqalt == "N" :
                start = start - 1
                seqref = sequences.fetch(chr, start - 1, end)
                seqalt = sequences.fetch(chr, start - 1, start)
                rec.alleles = (seqref, seqalt)
            svlen = -(len(seqref) - 1)         
            rec.info["END"] = end
            rec.info["SVLEN"] = svlen

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
        
        # Filter out SVs starting at pos 1
        # Mostly <INV> called by CuteSV, most likely erroneous calls
        # Produces errors in SVjedi-graph
        if rec.pos == 1:
            print("SV starting at position 1")
        else:
            out.write(rec)

