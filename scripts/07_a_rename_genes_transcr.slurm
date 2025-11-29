#!/usr/bin/env bash
# =============================================================================
# Script: 07_a_rename_genes_transcr.slurm
# Purpose: Rename gene and transcript IDs in MAKER output files to use a
#          consistent naming scheme with accession prefix
# =============================================================================

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=1:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=rename_genes_transcr
#SBATCH --output=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/rename_genes_transcr_output_%j.o
#SBATCH --error=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/rename_genes_transcr_error_%j.o

COURSE_DIR=/data/courses/assembly-annotation-course/CDS_annotation
WORK_DIR=/data/users/aschmid/Euk-Genome-Organization-and-Annotation
ACCESSION=Ms-0
MAKER_BIN=$COURSE_DIR/softwares/Maker_v3.01.03/src/bin
IN_DIR=$WORK_DIR/gene_annotation
OUT_DIR=$IN_DIR/final

# Input files
protein=$IN_DIR/assembly.all.maker.proteins.fasta
transcript=$IN_DIR/assembly.all.maker.transcripts.fasta
gff=$IN_DIR/assembly.all.maker.noseq.gff

# Output files (renamed)
protein_ren=$OUT_DIR/assembly.all.maker.proteins.renamed.fasta
transcript_ren=$OUT_DIR/assembly.all.maker.transcripts.renamed.fasta
gff_ren=$OUT_DIR/assembly.all.maker.noseq.renamed.gff

# Create Output Directory and Copy Files
mkdir -p $OUT_DIR

# Copy files to output directory for renaming
cp $gff $gff_ren
cp $protein $protein_ren
cp $transcript $transcript_ren

cd $OUT_DIR

# Generate ID mapping file
$MAKER_BIN/maker_map_ids --prefix $ACCESSION --justify 7 $gff_ren > id.map

# Apply ID mapping to GFF file
$MAKER_BIN/map_gff_ids id.map $gff_ren

# Apply ID mapping to protein FASTA file
$MAKER_BIN/map_fasta_ids id.map $protein_ren

# Apply ID mapping to transcript FASTA file
$MAKER_BIN/map_fasta_ids id.map $transcript_ren