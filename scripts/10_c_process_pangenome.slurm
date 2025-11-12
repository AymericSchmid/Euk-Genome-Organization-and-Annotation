#!/usr/bin/env bash

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

module load $R_MODULE

cd $OUT_DIR
Rscript $R_SCRIPT $GENESPACE_WORK_DIR Ms_0 $OUT_PLOT_FILE