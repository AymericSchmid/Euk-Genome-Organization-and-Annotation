#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=1:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=iprscan_annotation
#SBATCH --output=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/iprscan_annotation_output_%j.o
#SBATCH --error=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/iprscan_annotation_error_%j.o

COURSE_DIR=/data/courses/assembly-annotation-course/CDS_annotation
WORK_DIR=/data/users/aschmid/Euk-Genome-Organization-and-Annotation

IPRSCAN_IMG=$COURSE_DIR/containers/interproscan_latest.sif
IPRSCAN_DATA_DIR=$COURSE_DIR/data/interproscan-5.70-102.0/data
OUT_IPRSCAN_DIR=$WORK_DIR/gene_annotation/iprscan.output

MAKER_BIN=$COURSE_DIR/softwares/Maker_v3.01.03/src/bin

PROTEIN_FILE=$WORK_DIR/gene_annotation/final/assembly.all.maker.proteins.renamed.fasta
GFF_FILE_IN=$WORK_DIR/gene_annotation/final/assembly.all.maker.noseq.renamed.gff
GFF_FILE_OUT=$WORK_DIR/gene_annotation/final/assembly.all.maker.noseq.renamed.iprscan.gff

mkdir -p $OUT_IPRSCAN_DIR

# InterProScan to annotate protein sequences with functional domains, such as those from the Pfam database
apptainer exec \
    --bind $IPRSCAN_DATA_DIR:/opt/interproscan/data \
    --bind /data \
    --bind $SCRATCH:/temp \
    $IPRSCAN_IMG \
    /opt/interproscan/interproscan.sh \
    -appl pfam --disable-precalc -f TSV \
    --goterms --iprlookup --seqtype p \
    -i $PROTEIN_FILE -o $OUT_IPRSCAN_DIR

# Update GFF with InterProScan Results
$MAKER_BIN/ipr_update_gff $GFF_FILE_IN $OUT_IPRSCAN_DIR > $GFF_FILE_OUT