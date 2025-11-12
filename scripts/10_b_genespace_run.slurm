#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=32G
#SBATCH --time=3:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=genespace_run
#SBATCH --output=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/genespace_run_output_%j.o
#SBATCH --error=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/genespace_run_error_%j.o

COURSE_DIR=/data/courses/assembly-annotation-course/CDS_annotation
WORK_DIR=/data/users/aschmid/Euk-Genome-Organization-and-Annotation
GENESPACE_WORK_DIR=$WORK_DIR/comparative_genomics

GENESPACE_IMG=$COURSE_DIR/containers/genespace_latest.sif
R_SCRIPT=$WORK_DIR/scripts/RAnalysis/10_b_run_genespace.R
MCSCANX_PATH=$COURSE_DIR/softwares/MCScanX

cd $GENESPACE_WORK_DIR
apptainer exec --bind /data --bind $SCRATCH:/temp \
    $GENESPACE_IMG Rscript \
    $R_SCRIPT $GENESPACE_WORK_DIR $MCSCANX_PATH