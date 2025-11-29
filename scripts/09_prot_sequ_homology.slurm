#!/usr/bin/env bash
# =============================================================================
# Script: 09_prot_sequ_homology.slurm
# Purpose: Perform BLAST searches against protein databases (UniProt Viridiplantae
#          and TAIR10) to identify homologous proteins and annotate gene functions
# =============================================================================

#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=16G
#SBATCH --time=3:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=BLAST_prot_sequ_homology
#SBATCH --output=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/BLAST_prot_sequ_homology_output_%j.o
#SBATCH --error=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/BLAST_prot_sequ_homology_error_%j.o

COURSE_DIR="/data/courses/assembly-annotation-course/CDS_annotation"
WORK_DIR="/data/users/aschmid/Euk-Genome-Organization-and-Annotation"
OUT_DIR=$WORK_DIR/gene_annotation/functional_annotation
BLAST_OUTPUT=blastp.viridiplantae
BLAST_TAIR10_OUTPUT=blastp.tair10

# BLAST databases
DATABASE=$COURSE_DIR/data/uniprot/uniprot_viridiplantae_reviewed.fa 
DATABASE_TAIR10=$COURSE_DIR/data/TAIR10_pep_20110103_representative_gene_model
DATABASE_TAIR10_INDEXED=$WORK_DIR/gene_annotation/functional_annotation/tair10

# Input files
MAKER_PROTEINS=$WORK_DIR/gene_annotation/final/assembly.all.maker.proteins.renamed.filtered.fasta
MAKER_PROTEINS_UNIPROT=${MAKER_PROTEINS}.uniprot
GFF_GENES=$WORK_DIR/gene_annotation/final/filtered.genes.renamed.gff
GFF_GENES_UNIPROT=${GFF_GENES}.uniprot

BLAST_MODULE=BLAST+/2.15.0-gompi-2021a
MAKER_BIN=$COURSE_DIR/softwares/Maker_v3.01.03/src/bin

# Create Output Directory
mkdir -p $OUT_DIR
cd $OUT_DIR

# Load BLAST Module
module load $BLAST_MODULE

# BLAST Against UniProt Viridiplantae Database
echo "Blastp using viridiplantae db"
# Search against UniProt Viridiplantae (green plants) reviewed proteins
blastp -query $MAKER_PROTEINS -db $DATABASE -num_threads $SLURM_CPUS_PER_TASK -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -out $BLAST_OUTPUT

# Extract best hit per query sequence
# Sort by query ID, then by E-value (ascending), then keep unique queries
sort -k1,1 -k12,12g $BLAST_OUTPUT | sort -u -k1,1 --merge > ${BLAST_OUTPUT}.besthits

# Map Functional Annotations to GFF and FASTA Files
echo "Map protein putative functions to the maker GFF and FASTA files"
# Copy original files
cp $MAKER_PROTEINS $MAKER_PROTEINS_UNIPROT
cp $GFF_GENES $GFF_GENES_UNIPROT

# Add functional annotations from UniProt to FASTA headers
$MAKER_BIN/maker_functional_fasta $DATABASE ${BLAST_OUTPUT}.besthits $MAKER_PROTEINS > $MAKER_PROTEINS_UNIPROT

# Add functional annotations from UniProt to GFF attributes
$MAKER_BIN/maker_functional_gff $DATABASE $BLAST_OUTPUT $GFF_GENES > $GFF_GENES_UNIPROT

# Prepare TAIR10 Database
echo "Prepare TAIR10 db"
# Create BLAST database from TAIR10 protein sequences
makeblastdb -in $DATABASE_TAIR10 -dbtype prot -out $DATABASE_TAIR10_INDEXED

# BLAST Against TAIR10 Database
echo "Blastp using TAIR10 db"
# Search against Arabidopsis thaliana TAIR10 representative gene models
# This helps identify orthologs in a well-annotated reference genome
blastp -query $MAKER_PROTEINS -db $DATABASE_TAIR10_INDEXED -num_threads $SLURM_CPUS_PER_TASK -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -out $BLAST_TAIR10_OUTPUT

# Extract best hit per query sequence
sort -k1,1 -k12,12g $BLAST_TAIR10_OUTPUT | sort -u -k1,1 --merge > ${BLAST_TAIR10_OUTPUT}.besthits