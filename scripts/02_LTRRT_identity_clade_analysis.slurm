#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=1:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=LTRRT_identity_clade
#SBATCH --output=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/LTRRT_identity_clade_output_%j.o
#SBATCH --error=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/LTRRT_identity_clade_error_%j.o

WORK_DIR=/data/users/aschmid/Euk-Genome-Organization-and-Annotation
OUT_DIR=$WORK_DIR/TE_annotation/clades/LTRRT
PLOT_DIR=$WORK_DIR/plots
LTRRT_ANNOTATION=$WORK_DIR/TE_annotation/EDTA/assembly.fa.mod.EDTA.raw/assembly.fa.mod.LTR.intact.raw.gff3
LTRRT_LIBRARY=$WORK_DIR/TE_annotation/EDTA/assembly.fa.mod.EDTA.raw/assembly.fa.mod.LTR.raw.fa
TESORTER_IMG=/data/courses/assembly-annotation-course/CDS_annotation/containers/TEsorter_1.3.0.sif
R_MODULE=R-bundle-IBU/2023072800-foss-2021a-R-4.2.1
TSV_FILE=$OUT_DIR/assembly.fa.mod.LTR.raw.fa.rexdb-plant.cls.tsv

mkdir -p $PLOT_DIR
mkdir -p $OUT_DIR

cd $OUT_DIR

# Refining TE classification for full length LTR-RTs and split them into Clades
apptainer exec --bind /data $TESORTER_IMG \
TEsorter \
$LTRRT_LIBRARY -db rexdb-plant

cd $PLOT_DIR

# Plot the number of LTR-RTs in each clade with their corresponding percent identity
module load $R_MODULE 
Rscript $WORK_DIR/scripts/RAnalysis/02_full_length_LTRs_identity.R \
$LTRRT_ANNOTATION $TSV_FILE \