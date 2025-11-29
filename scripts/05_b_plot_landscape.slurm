#!/usr/bin/env bash
# =============================================================================
# Script: 05_b_plot_landscape.slurm
# Purpose: Generate TE divergence landscape plot from parsed RepeatMasker data
#          Shows the distribution of TE sequences by divergence percentage
# =============================================================================

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=1:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=landscape_plot
#SBATCH --output=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/landscape_plot_output_%j.o
#SBATCH --error=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/landscape_plot_error_%j.o

WORK_DIR=/data/users/aschmid/Euk-Genome-Organization-and-Annotation
DIV_DATA=$WORK_DIR/divergence/assembly.fa.mod.out.landscape.Div.Rname.tab
PLOT_DIR=$WORK_DIR/plots
R_MODULE=R-bundle-IBU/2023072800-foss-2021a-R-4.2.1

# Create Output Directory
mkdir -p $PLOT_DIR
cd $PLOT_DIR

module load $R_MODULE 

# Run R script to create divergence landscape plot
# Argument: path to divergence data table
Rscript $WORK_DIR/scripts/RAnalysis/05_plot_div.R \
$DIV_DATA