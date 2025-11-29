#!/usr/bin/env bash
# =============================================================================
# Script: 06_c_maker_output_prep.slurm
# Purpose: Merge MAKER output files from distributed datastore into single
#          GFF and FASTA files for downstream analysis
# =============================================================================

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=1:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=maker_output_prep
#SBATCH --output=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/maker_output_prep_output_%j.o
#SBATCH --error=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/maker_output_prep_error_%j.o

COURSE_DIR=/data/courses/assembly-annotation-course/CDS_annotation
WORK_DIR=/data/users/aschmid/Euk-Genome-Organization-and-Annotation
MAKER_OUTPUT_DIR=$WORK_DIR/gene_annotation/assembly.maker.output
MAKER_BIN=$COURSE_DIR/softwares/Maker_v3.01.03/src/bin
DATASTORE_INDEX_LOG=$MAKER_OUTPUT_DIR/assembly_master_datastore_index.log

OUT_DIR=$WORK_DIR/gene_annotation
OUT_ALL_GFF=$OUT_DIR/assembly.all.maker.gff
OUT_NOSEQ_ALL_GFF=$OUT_DIR/assembly.all.maker.noseq.gff

# Create Output Directory
mkdir -p $OUT_DIR
cd $OUT_DIR

# Merge GFF files with sequences included
$MAKER_BIN/gff3_merge -s -d $DATASTORE_INDEX_LOG > $OUT_ALL_GFF

# Merge GFF files without sequences (smaller file, faster for some analyses)
$MAKER_BIN/gff3_merge -n -s -d $DATASTORE_INDEX_LOG > $OUT_NOSEQ_ALL_GFF

# Merge FASTA files (proteins and transcripts)
$MAKER_BIN/fasta_merge -d $DATASTORE_INDEX_LOG -o assembly