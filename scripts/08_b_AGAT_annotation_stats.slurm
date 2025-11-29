#!/usr/bin/env bash
# =============================================================================
# Script: 08_b_AGAT_annotation_stats.slurm
# Purpose: Generate comprehensive statistics for gene annotation using AGAT
#          (Another GFF Analysis Toolkit)
# =============================================================================

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=1:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=AGAT_annotation_stats
#SBATCH --output=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/AGAT_annotation_stats_output_%j.o
#SBATCH --error=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/AGAT_annotation_stats_error_%j.o

COURSE_DIR="/data/courses/assembly-annotation-course/CDS_annotation"
WORK_DIR="/data/users/aschmid/Euk-Genome-Organization-and-Annotation"
GFF_FILE=$WORK_DIR/gene_annotation/final/filtered.genes.renamed.gff
OUT_DIR=$WORK_DIR/annotation_qc/agat
OUT_FILE=$OUT_DIR/filtered.genes.renamed.gff.stat

AGAT_IMG=/containers/apptainer/agat-1.2.0.sif

# Create Output Directory
mkdir -p $OUT_DIR

# Run AGAT to generate comprehensive annotation statistics
apptainer exec \
    --bind /data \
    $AGAT_IMG \
    agat_sp_statistics.pl -i $GFF_FILE -o $OUT_FILE