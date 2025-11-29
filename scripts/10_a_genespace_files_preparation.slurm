#!/usr/bin/env bash
# =============================================================================
# Script: 10_a_genespace_files_preparation.slurm
# Purpose: Prepare input files (BED and FASTA) for GENESPACE comparative
#          genomics analysis from multiple accessions
# =============================================================================

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=3:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=genespace_files_preparation
#SBATCH --output=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/genespace_files_preparation_output_%j.o
#SBATCH --error=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/genespace_files_preparation_error_%j.o

COURSE_DIR=/data/courses/assembly-annotation-course/CDS_annotation
WORK_DIR=/data/users/aschmid/Euk-Genome-Organization-and-Annotation
OUT_DIR=$WORK_DIR/comparative_genomics
PEPTIDE_DIR=$OUT_DIR/peptide
BED_DIR=$OUT_DIR/bed

# Input files for focal genome (Ms-0)
GFF_GENES=$WORK_DIR/gene_annotation/final/filtered.genes.renamed.gff
PROT_SEQUENCES=$WORK_DIR/gene_annotation/final/assembly.all.maker.proteins.renamed.filtered.longest.fasta

# Output files for focal genome
BED_FILE=$BED_DIR/Ms_0.bed
FA_FILE=$PEPTIDE_DIR/Ms_0.fa

# Reference genome files (TAIR10)
REF_BED_FILE_SOURCE=$COURSE_DIR/data/TAIR10.bed
REF_BED_FILE=$BED_DIR/TAIR10.bed
REF_FA_FILE_SOURCE=$COURSE_DIR/data/TAIR10.fa
REF_FA_FILE=$PEPTIDE_DIR/TAIR10.fa

# Additional accessions to include in analysis
ACCESSIONS=(Altai-5 Are-6 Ice-1)
ACCESSIONS_FA_DIR=$COURSE_DIR/data/Lian_et_al/protein/selected
ACCESSIONS_GFF_DIR=$COURSE_DIR/data/Lian_et_al/gene_gff/selected

# Create Output Directories
mkdir -p $OUT_DIR $PEPTIDE_DIR $BED_DIR

# Extract positions of each gene from the filtered gff file
# Also sort by position as the GFF is not sorted for this accession
grep -P "\tgene\t" $GFF_GENES \
| awk 'BEGIN{OFS="\t"} {split($9,a,";"); split(a[1],b,"="); print $1, $4-1, $5, b[2]}' \
| sort -k1,1 -k2,2n \
> $BED_FILE

# Clean protein FASTA headers for GENESPACE
# Remove transcript suffixes (e.g., -RA, -RB) and keep only gene ID
# Example: >Ms-00029094-RA protein AED:0.36 ... -> >Ms-00029094
awk 'BEGIN{FS=" "; OFS=" "}
  /^>/ {
    # take only the first "word" of the header
    id = $1
    # remove leading ">"
    sub(/^>/, "", id)
    # drop transcript suffix like -RA, -RB, -RC
    sub(/-R[A-Z]$/, "", id)
    # print back with ">"
    print ">" id
    next
  }
  {print}
' $PROT_SEQUENCES > $FA_FILE

# Copy Reference Genome Files (TAIR10)
cp $REF_BED_FILE_SOURCE $REF_BED_FILE
cp $REF_FA_FILE_SOURCE $REF_FA_FILE

# Process each additional accession
for acc in "${ACCESSIONS[@]}"; do
    # Replace - with _ in accession name (GENESPACE requirement)
    acc_clean=${acc//-/_}

    # Define input and output paths
    gff_in=$ACCESSIONS_GFF_DIR/${acc}.EVM.v3.5.ann.protein_coding_genes.gff
    bed_out=$BED_DIR/${acc_clean}.bed
    fa_in=$ACCESSIONS_FA_DIR/${acc}.protein.faa
    fa_out=$PEPTIDE_DIR/${acc_clean}.fa

    echo "Processing $acc ..."

    # Extract gene positions and convert to BED format
    grep -P "\tgene\t" $gff_in | \
    awk 'BEGIN{OFS="\t"} {split($9,a,";"); split(a[1],b,"="); print $1, $4-1, $5, b[2]}' \
    > $bed_out

    # Copy protein FASTA file
    cp $fa_in $fa_out
done

echo "All files prepared in: $OUT_DIR"