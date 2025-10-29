#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=1:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=TE_divergence
#SBATCH --output=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/TE_divergence_output_%j.o
#SBATCH --error=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/TE_divergence_error_%j.o

WORK_DIR=/data/users/aschmid/Euk-Genome-Organization-and-Annotation
OUT_DIR=$WORK_DIR/divergence
REPEATMASKER_FILE=$WORK_DIR/TE_annotation/EDTA/assembly.fa.mod.EDTA.anno/assembly.fa.mod.out
PERL_SCRIPT=$WORK_DIR/scripts/perlScripts/parseRM.pl

PERL_MODULE=BioPerl/1.7.8-GCCcore-10.3.0

mkdir -p $OUT_DIR

module add $PERL_MODULE
perl $PERL_SCRIPT -i $REPEATMASKER_FILE -l 50,1 -v

mv $REPEATMASKER_FILE.landscape* $OUT_DIR