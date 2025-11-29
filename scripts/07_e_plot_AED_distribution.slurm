#!/usr/bin/env bash
# =============================================================================
# Script: 07_e_plot_AED_distribution.slurm
# Purpose: Visualize the distribution of Annotation Edit Distance (AED) scores
#          for gene annotations, showing annotation quality
# =============================================================================

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=1:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=AED_distribution_plot
#SBATCH --output=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/AED_distribution_plot_output_%j.o
#SBATCH --error=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/AED_distribution_plot_error_%j.o

WORK_DIR=/data/users/aschmid/Euk-Genome-Organization-and-Annotation
AED_DATA=$WORK_DIR/gene_annotation/final/filtered.genes.renamed.gff.AED.txt
PLOT_DIR=$WORK_DIR/plots
R_MODULE=R-bundle-IBU/2023072800-foss-2021a-R-4.2.1

# Create Output Directory
mkdir -p $PLOT_DIR
cd $PLOT_DIR

module load $R_MODULE 

# Run R script to create AED distribution plots (histogram and cumulative)
# Arguments: AED data file, output PDF filename
Rscript $WORK_DIR/scripts/RAnalysis/07_genes_AED_distribution.R \
$AED_DATA \
"07_AED_distribution.pdf"