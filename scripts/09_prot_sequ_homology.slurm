#!/usr/bin/env bash

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

DATABASE=$COURSE_DIR/data/uniprot/uniprot_viridiplantae_reviewed.fa 
DATABASE_TAIR10=$COURSE_DIR/data/TAIR10_pep_20110103_representative_gene_model
DATABASE_TAIR10_INDEXED=$WORK_DIR/gene_annotation/functional_annotation/tair10
MAKER_PROTEINS=$WORK_DIR/gene_annotation/final/assembly.all.maker.proteins.renamed.filtered.fasta
MAKER_PROTEINS_UNIPROT=${MAKER_PROTEINS}.uniprot
GFF_GENES=$WORK_DIR/gene_annotation/final/filtered.genes.renamed.gff
GFF_GENES_UNIPROT=${GFF_GENES}.uniprot

BLAST_MODULE=BLAST+/2.15.0-gompi-2021a
MAKER_BIN=$COURSE_DIR/softwares/Maker_v3.01.03/src/bin

mkdir -p $OUT_DIR
cd $OUT_DIR

module load $BLAST_MODULE

echo "Blastp using viridiplantae db"
blastp -query $MAKER_PROTEINS -db $DATABASE -num_threads $SLURM_CPUS_PER_TASK -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -out $BLAST_OUTPUT
# Now sort the blast output to keep only the best hit per query sequence
sort -k1,1 -k12,12g $BLAST_OUTPUT | sort -u -k1,1 --merge > ${BLAST_OUTPUT}.besthits

echo "Map protein putative functions to the maker GFF and FASTA files"
# map the protein putative functions to the MAKER produced GFF3 and FASTA files
cp $MAKER_PROTEINS $MAKER_PROTEINS_UNIPROT
cp $GFF_GENES $GFF_GENES_UNIPROT
$MAKER_BIN/maker_functional_fasta $DATABASE ${BLAST_OUTPUT}.besthits $MAKER_PROTEINS > $MAKER_PROTEINS_UNIPROT
$MAKER_BIN/maker_functional_gff $DATABASE $BLAST_OUTPUT $GFF_GENES > $GFF_GENES_UNIPROT

echo "Prepare TAIR10 db"
# Prepare tair10 db
makeblastdb -in $DATABASE_TAIR10 -dbtype prot -out $DATABASE_TAIR10_INDEXED

echo "Blastp using TAIR10 db"
# get the best blast hit with Arabidopsis thaliana TAIR10 representative gene models
blastp -query $MAKER_PROTEINS -db $DATABASE_TAIR10_INDEXED -num_threads $SLURM_CPUS_PER_TASK -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -out $BLAST_TAIR10_OUTPUT
# Now sort the blast output to keep only the best hit per query sequence
sort -k1,1 -k12,12g $BLAST_TAIR10_OUTPUT | sort -u -k1,1 --merge > ${BLAST_TAIR10_OUTPUT}.besthits