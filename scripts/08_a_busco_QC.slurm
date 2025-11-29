#!/usr/bin/env bash
# =============================================================================
# Script: 08_a_busco_QC.slurm
# Purpose: Assess annotation quality using BUSCO (Benchmarking Universal
#          Single-Copy Orthologs) on protein and transcript sequences
# =============================================================================

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
OUT_BUSCO_PLOT_DIR=$OUT_BUSCO_DIR/plots
ASSEMBLY_BUSCO_RESULT=/data/users/aschmid/Genome-and-Transcriptome-Assembly/assembly_QC/busco/hifiasm/short_summary.specific.brassicales_odb10.hifiasm.txt

SEQKIT_MODULE=SeqKit/2.6.1
BUSCO_MODULE=BUSCO/5.4.2-foss-2021a

# Create Output Directory
mkdir -p $OUT_BUSCO_DIR

# Load Required Modules
module load $SEQKIT_MODULE
module load $BUSCO_MODULE

# For BUSCO analysis, use only the longest protein/transcript per gene
# to avoid inflating completeness scores with multiple isoforms

# Extract longest protein per gene
# 1. Get sequence names and lengths
# 2. Split ID to get gene name (remove transcript suffix like -RA, -RB)
# 3. Sort by gene, then by length (descending), then by ID
# 4. Keep only first occurrence per gene (longest)
seqkit fx2tab -n -l $PROTEIN_FASTA \
    | awk -F'\t' '{split($1,a,/-R/); gene=a[1]; print gene"\t"$1"\t"$2}' \
    | sort -k1,1 -k3,3nr -k2,2 \
    | awk '!seen[$1]++{print $2}' > $IDS_LONGEST_PROTS
seqkit grep -f $IDS_LONGEST_PROTS $PROTEIN_FASTA > $PROTEIN_FASTA_LONGEST

# Extract longest transcript per gene
seqkit fx2tab -n -l $TRANSCRIPT_FASTA \
    | awk -F'\t' '{split($1,a,/-R/); gene=a[1]; print gene"\t"$1"\t"$2}' \
    | sort -k1,1 -k3,3nr -k2,2 \
    | awk '!seen[$1]++{print $2}' > $IDS_LONGEST_TRANS
seqkit grep -f $IDS_LONGEST_TRANS $TRANSCRIPT_FASTA > $TRANSCRIPT_FASTA_LONGEST

cd $OUT_BUSCO_DIR

# Run BUSCO on protein sequences
busco -i $PROTEIN_FASTA_LONGEST -l brassicales_odb10 -o proteins -m proteins

# Run BUSCO on transcript sequences
busco -i $TRANSCRIPT_FASTA_LONGEST -l brassicales_odb10 -o transcriptome -m transcriptome

mkdir -p $OUT_BUSCO_PLOT_DIR

# Copy BUSCO summary files for comparison plotting
cp $ASSEMBLY_BUSCO_RESULT $OUT_BUSCO_PLOT_DIR
cp $OUT_BUSCO_DIR/proteins/short_summary.specific.brassicales_odb10.proteins.txt $OUT_BUSCO_PLOT_DIR
cp $OUT_BUSCO_DIR/transcriptome/short_summary.specific.brassicales_odb10.transcriptome.txt $OUT_BUSCO_PLOT_DIR

# Generate BUSCO comparison plot
PLOT_SCRIPT=$(dirname $(which busco))/generate_plot.py
python3 $PLOT_SCRIPT -wd $OUT_BUSCO_PLOT_DIR