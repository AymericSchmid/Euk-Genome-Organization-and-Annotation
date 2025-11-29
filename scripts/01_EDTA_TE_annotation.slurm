#!/usr/bin/env bash
# =============================================================================
# Script: 01_EDTA_TE_annotation.slurm
# Purpose: Annotate transposable elements (TEs) in the genome using EDTA
#          (Extensive de-novo TE Annotator)
# =============================================================================

#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=16G
#SBATCH --time=24:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=EDTA_TE_annotation
#SBATCH --output=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/EDTA_TE_annotation_output_%j.o
#SBATCH --error=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/EDTA_TE_annotation_error_%j.o

WORK_DIR=/data/users/aschmid/Euk-Genome-Organization-and-Annotation
OUT_DIR=$WORK_DIR/TE_annotation/EDTA
EDTA_IMG=/data/courses/assembly-annotation-course/CDS_annotation/containers/EDTA2.2.sif
ASSEMBLY=$WORK_DIR/assembly/assembly.fa
CDS_FILE=/data/courses/assembly-annotation-course/CDS_annotation/data/TAIR10_cds_20110103_representative_gene_model_updated

# Create Output Directory
mkdir -p $OUT_DIR

# Change to output directory (EDTA writes output to current directory)
cd $OUT_DIR

apptainer exec --bind /data $EDTA_IMG \
EDTA.pl \
--genome $ASSEMBLY \
--species others \
--step all \
--sensitive 1 \
--cds $CDS_FILE \
--anno 1 \
--threads $SLURM_CPUS_PER_TASK