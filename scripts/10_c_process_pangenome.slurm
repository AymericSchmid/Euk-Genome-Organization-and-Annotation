#!/usr/bin/env bash
# =============================================================================
# Script: 10_c_process_pangenome.slurm
# Purpose: Process pangenome data from GENESPACE to calculate core, accessory,
#          and species-specific genes and generate frequency plots
# =============================================================================

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=3:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=process_pangenome
#SBATCH --output=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/process_pangenome_output_%j.o
#SBATCH --error=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/process_pangenome_error_%j.o

COURSE_DIR=/data/courses/assembly-annotation-course/CDS_annotation
WORK_DIR=/data/users/aschmid/Euk-Genome-Organization-and-Annotation
GENESPACE_WORK_DIR=$WORK_DIR/comparative_genomics
OUT_DIR=$WORK_DIR/plots
OUT_PLOT_FILE=10_c_pangenome_frequency_plot.pdf

R_SCRIPT=$WORK_DIR/scripts/RAnalysis/10_c_process_pangenome.R
R_MODULE=R-bundle-IBU/2023072800-foss-2021a-R-4.2.1

# Load R Module
module load $R_MODULE

cd $OUT_DIR

# Run R script to process pangenome data and generate frequency plots
Rscript $R_SCRIPT $GENESPACE_WORK_DIR Ms_0 $OUT_PLOT_FILE