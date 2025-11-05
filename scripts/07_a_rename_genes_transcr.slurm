#!/usr/bin/env bash

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

mkdir -p $OUT_DIR

protein=$IN_DIR/assembly.all.maker.proteins.fasta
transcript=$IN_DIR/assembly.all.maker.transcripts.fasta
gff=$IN_DIR/assembly.all.maker.noseq.gff

protein_ren=$OUT_DIR/assembly.all.maker.proteins.renamed.fasta
transcript_ren=$OUT_DIR/assembly.all.maker.transcripts.renamed.fasta
gff_ren=$OUT_DIR/assembly.all.maker.noseq.renamed.gff

cp $gff $gff_ren
cp $protein $protein_ren
cp $transcript $transcript_ren

cd $OUT_DIR

$MAKER_BIN/maker_map_ids --prefix $ACCESSION --justify 7 $gff_ren > id.map
$MAKER_BIN/map_gff_ids id.map $gff_ren
$MAKER_BIN/map_fasta_ids id.map $protein_ren
$MAKER_BIN/map_fasta_ids id.map $transcript_ren