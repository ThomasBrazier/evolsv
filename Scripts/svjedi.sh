#!/bin/bash
sra=$1

# Format DUP as INS for SVjedi
# Duplication must be defined as an insertion event whith CHR and POS corresponding to the position of insertion of the novel copy
# INFO field must contain SVTYPE=INS
# ALT field must contain the sequence of the duplication
sed 's/SVTYPE=DUP[:A-Z]*/SVTYPE=INS/g'  $sra/${sra}_merged.vcf >  $sra/${sra}_merged_formatSVjedi.vcf
# Extract INFO/DP into a tab-delimited annotation file
bcftools query -f '%CHROM\t%POS\t%DP\n' test.vcf | bgzip -c > annot.txt.gz

# Index the file with tabix
tabix -s1 -b2 -e2 annot.txt.gz

# Create a header line for the new annotation
echo -e '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">' >> hdr.txt

# Transfer the annotation to sample 'smpl1'
bcftools annotate -s smpl1 -a annot.txt.gz -h hdr.txt -c CHROM,POS,FORMAT/DP test.vcf


# Genotyping
svjedi -v $sra/${sra}_merged_formatSVjedi.vcf -r $sra/ref_${sra}.fna -i $sra/${sra}.fastq.gz -o $sra/${sra}_genotypes.vcf
