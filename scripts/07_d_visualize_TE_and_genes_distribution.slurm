#!/usr/bin/env bash
# =============================================================================
# Script: 07_d_visualize_TE_and_genes_distribution.slurm
# Purpose: Visualize the distribution of both transposable elements and genes
#          across the genome using circular plots
# =============================================================================

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=1:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=visualize_TE_and_genes_distribution
#SBATCH --output=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/visualize_TE_and_genes_distribution_output_%j.o
#SBATCH --error=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/visualize_TE_and_genes_distribution_error_%j.o

WORK_DIR=/data/users/aschmid/Euk-Genome-Organization-and-Annotation
ASSEMBLY_FAI=$WORK_DIR/assembly/assembly.fa.fai
TE_ANNOTATION=$WORK_DIR/TE_annotation/TEanno.gff3
GENE_ANNOTATION=$WORK_DIR/gene_annotation/final/filtered.genes.renamed.gff
PLOT_DIR=$WORK_DIR/plots
R_MODULE=R-bundle-IBU/2023072800-foss-2021a-R-4.2.1
OUT_FILE=07_TE_and_genes_distribution.pdf

# Create Output Directory
mkdir -p $PLOT_DIR

module load $R_MODULE 
cd $PLOT_DIR

# Run R script to create circular plot showing both TE and gene distribution
Rscript $WORK_DIR/scripts/RAnalysis/04_annotation_circlize.R \
$ASSEMBLY_FAI $TE_ANNOTATION $GENE_ANNOTATION $OUT_FILE