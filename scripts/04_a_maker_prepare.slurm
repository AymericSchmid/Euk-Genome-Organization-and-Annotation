#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=1:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=maker_control_file
#SBATCH --output=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/maker_control_file_output_%j.o
#SBATCH --error=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/maker_control_file_error_%j.o

WORK_DIR=/data/users/aschmid/Euk-Genome-Organization-and-Annotation
OUT_DIR=$WORK_DIR/gene_annotation
MAKER_IMG=/data/courses/assembly-annotation-course/CDS_annotation/containers/MAKER_3.01.03.sif

mkdir -p $OUT_DIR
cd $OUT_DIR

apptainer exec --bind /data $MAKER_IMG \
maker -CTL