#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:05:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=dir_set_up
#SBATCH --output=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/dir_set_up_output_%j.o
#SBATCH --error=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/dir_set_up_error_%j.o

WORK_DIR=/data/users/aschmid/Euk-Genome-Organization-and-Annotation
ASSEMBLY_SOURCE=/data/users/aschmid/Genome-and-Transcriptome-Assembly/assembly/hifiasm/asm.bp.p_ctg.fa
ASSEMBLY_DEST=$WORK_DIR/assembly
TRINITY_ASSEMBLY_SOURCE=/data/users/aschmid/Genome-and-Transcriptome-Assembly/assembly/trinity.Trinity.fasta
TRINITY_ASSEMBLY_DEST=$WORK_DIR/assembly
HIFI_READS_SOURCE=/data/users/aschmid/Genome-and-Transcriptome-Assembly/reads/Ms-0/ERR11437313.fastq.gz
HIFI_READS_DEST=$WORK_DIR/reads/Ms-0
RNASEQ_READS_SOURCE=/data/users/aschmid/Genome-and-Transcriptome-Assembly/reads_filt/ERR754081_*.trimmed.fastq.gz
RNASEQ_READS_DEST=$WORK_DIR/reads/RNAseq_Sha

mkdir -p $ASSEMBLY_DEST
mkdir -p $TRINITY_ASSEMBLY_DEST
mkdir -p $HIFI_READS_DEST
mkdir -p $RNASEQ_READS_DEST

ln -s $ASSEMBLY_SOURCE $ASSEMBLY_DEST/assembly.fa
ln -s $HIFI_READS_SOURCE $HIFI_READS_DEST
ln -s $RNASEQ_READS_SOURCE $RNASEQ_READS_DEST
ln -s $TRINITY_ASSEMBLY_SOURCE $TRINITY_ASSEMBLY_DEST/trinity.fa