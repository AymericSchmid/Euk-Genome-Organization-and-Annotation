#!/usr/bin/env bash
# =============================================================================
# Script: 04_visualize_TE_distribution.slurm
# Purpose: Visualize the distribution of transposable elements across the
#          genome using circular plots (circlize)
# =============================================================================

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=1:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=visualize_TE_distribution
#SBATCH --output=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/visualize_TE_distribution_output_%j.o
#SBATCH --error=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/visualize_TE_distribution_error_%j.o

WORK_DIR=/data/users/aschmid/Euk-Genome-Organization-and-Annotation
ASSEMBLY=$WORK_DIR/assembly/assembly.fa
ASSEMBLY_FAI=$WORK_DIR/assembly/assembly.fa.fai
TE_ANNOTATION=$WORK_DIR/TE_annotation/TEanno.gff3
PLOT_DIR=$WORK_DIR/plots
R_MODULE=R-bundle-IBU/2023072800-foss-2021a-R-4.2.1
SAMTOOLS_IMG=/containers/apptainer/samtools-1.19.sif 
OUT_FILE=04_TE_distribution.pdf

# Create Output Directory
mkdir -p $PLOT_DIR

# Index Assembly (Required for R Script)
# Create FASTA index file (.fai) needed for the visualization script
apptainer exec --bind /data $SAMTOOLS_IMG \
samtools faidx $ASSEMBLY

module load $R_MODULE 
cd $PLOT_DIR

# Run R script to create circular plot showing TE distribution
# Arguments: assembly index, TE annotation GFF, output PDF filename
Rscript $WORK_DIR/scripts/RAnalysis/04_annotation_circlize.R \
$ASSEMBLY_FAI $TE_ANNOTATION $OUT_FILE