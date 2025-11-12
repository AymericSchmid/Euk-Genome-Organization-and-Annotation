#!/usr/bin/env bash

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

GFF_GENES=$WORK_DIR/gene_annotation/final/filtered.genes.renamed.gff
PROT_SEQUENCES=$WORK_DIR/gene_annotation/final/assembly.all.maker.proteins.renamed.filtered.longest.fasta
BED_FILE=$BED_DIR/Ms_0.bed
FA_FILE=$PEPTIDE_DIR/Ms_0.fa
REF_BED_FILE_SOURCE=$COURSE_DIR/data/TAIR10.bed
REF_BED_FILE=$BED_DIR/TAIR10.bed
REF_FA_FILE_SOURCE=$COURSE_DIR/data/TAIR10.fa
REF_FA_FILE=$PEPTIDE_DIR/TAIR10.fa

ACCESSIONS=(Altai-5 Are-6 Ice-1)
ACCESSIONS_FA_DIR=$COURSE_DIR/data/Lian_et_al/protein/selected
ACCESSIONS_GFF_DIR=$COURSE_DIR/data/Lian_et_al/gene_gff/selected

mkdir -p $OUT_DIR $PEPTIDE_DIR $BED_DIR

# Extract positions of each gene from the filtered gff file
# Also sort by position as the GFF is not sorted for this accession
grep -P "\tgene\t" $GFF_GENES \
| awk 'BEGIN{OFS="\t"} {split($9,a,";"); split(a[1],b,"="); print $1, $4-1, $5, b[2]}' \
| sort -k1,1 -k2,2n \
> $BED_FILE

# Remove the suffix we find in xxx.longest.fasta. Example: >Ms-00029094-RA protein AED:0.36 eAED:0.36 QI:0|-1|0|1|-1|1|1|0|508 -> >Ms-00029094
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

cp $REF_BED_FILE_SOURCE $REF_BED_FILE
cp $REF_FA_FILE_SOURCE $REF_FA_FILE

for acc in "${ACCESSIONS[@]}"; do
    acc_clean=${acc//-/_}       # Replace the - in the accesion with a _

    gff_in=$ACCESSIONS_GFF_DIR/${acc}.EVM.v3.5.ann.protein_coding_genes.gff
    bed_out=$BED_DIR/${acc_clean}.bed
    fa_in=$ACCESSIONS_FA_DIR/${acc}.protein.faa
    fa_out=$PEPTIDE_DIR/${acc_clean}.fa

    echo "Processing $acc ..."

    grep -P "\tgene\t" $gff_in | \
    awk 'BEGIN{OFS="\t"} {split($9,a,";"); split(a[1],b,"="); print $1, $4-1, $5, b[2]}' \
    > $bed_out

    cp $fa_in $fa_out
done

echo All files prepared in: $OUT_DIR