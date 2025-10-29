#!/usr/bin/env bash

#SBATCH --time=7-0
#SBATCH --mem=120G
#SBATCH -p pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=50
#SBATCH --job-name=maker_gene_anno
#SBATCH --output=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/maker_gene_anno_output_%j.o
#SBATCH --error=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/maker_gene_anno_error_%j.o

COURSE_DIR="/data/courses/assembly-annotation-course/CDS_annotation"
WORK_DIR=/data/users/aschmid/Euk-Genome-Organization-and-Annotation
REPEATMASKER_DIR="/data/courses/assembly-annotation-course/CDS_annotation/softwares/RepeatMasker"
OUT_DIR=$WORK_DIR/gene_annotation
MAKER_IMG=/data/courses/assembly-annotation-course/CDS_annotation/containers/MAKER_3.01.03.sif

export PATH=$PATH:$REPEATMASKER_DIR

cd $OUT_DIR

module load OpenMPI/4.1.1-GCC-10.3.0
module load AUGUSTUS/3.4.0-foss-2021a

mpiexec --oversubscribe -n 50 apptainer exec \
--bind $SCRATCH:/TMP --bind $COURSE_DIR --bind $AUGUSTUS_CONFIG_PATH --bind $REPEATMASKER_DIR --bind /data \
$MAKER_IMG \
maker -mpi --ignore_nfs_tmp -TMP /TMP maker_opts.ctl maker_bopts.ctl maker_evm.ctl maker_exe.ctl