#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=32G
#SBATCH --time=10:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=busco_QC
#SBATCH --output=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/busco_QC_output_%j.o
#SBATCH --error=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/busco_QC_error_%j.o

COURSE_DIR="/data/courses/assembly-annotation-course/CDS_annotation"
WORK_DIR="/data/users/aschmid/Euk-Genome-Organization-and-Annotation"

TRANSCRIPT_FASTA="${WORK_DIR}/gene_annotation/final/assembly.all.maker.transcripts.renamed.filtered.fasta"
PROTEIN_FASTA="${WORK_DIR}/gene_annotation/final/assembly.all.maker.proteins.renamed.filtered.fasta"
TRANSCRIPT_FASTA_LONGEST="${WORK_DIR}/gene_annotation/final/assembly.all.maker.transcripts.renamed.filtered.longest.fasta"
PROTEIN_FASTA_LONGEST="${WORK_DIR}/gene_annotation/final/assembly.all.maker.proteins.renamed.filtered.longest.fasta"

IDS_LONGEST_PROTS=$SCRATCH/ids.prots.longest.txt
IDS_LONGEST_TRANS=$SCRATCH/ids.trans.longest.txt

OUT_BUSCO_DIR=$WORK_DIR/annotation_qc/busco

SEQKIT_MODULE=SeqKit/2.6.1
BUSCO_MODULE=BUSCO/5.4.2-foss-2021a

mkdir -p $OUT_BUSCO_DIR

module load $SEQKIT_MODULE
module load $BUSCO_MODULE

# Extract Longest Protein per Gene
seqkit fx2tab -n -l $PROTEIN_FASTA \
    | awk -F'\t' '{split($1,a,/-R/); gene=a[1]; print gene"\t"$1"\t"$2}' \
    | sort -k1,1 -k3,3nr -k2,2 \
    | awk '!seen[$1]++{print $2}' > $IDS_LONGEST_PROTS
seqkit grep -f $IDS_LONGEST_PROTS $PROTEIN_FASTA > $PROTEIN_FASTA_LONGEST

# Extract Longest Transcript per Gene
seqkit fx2tab -n -l $TRANSCRIPT_FASTA \
    | awk -F'\t' '{split($1,a,/-R/); gene=a[1]; print gene"\t"$1"\t"$2}' \
    | sort -k1,1 -k3,3nr -k2,2 \
    | awk '!seen[$1]++{print $2}' > $IDS_LONGEST_TRANS
seqkit grep -f $IDS_LONGEST_TRANS $TRANSCRIPT_FASTA > $TRANSCRIPT_FASTA_LONGEST

cd $OUT_BUSCO_DIR

# Run BUSCO
busco -i $PROTEIN_FASTA_LONGEST -l brassicales_odb10 -o proteins -m proteins
busco -i $TRANSCRIPT_FASTA_LONGEST -l brassicales_odb10 -o transcriptome -m transcriptome