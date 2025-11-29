#!/usr/bin/env bash
# =============================================================================
# Script: 06_b_gene_anno_maker.slurm
# Purpose: Run MAKER gene annotation pipeline using MPI for parallel processing
#          This is the main gene annotation step that combines evidence from
#          multiple sources (protein homology, ESTs, ab initio predictions)
# =============================================================================

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

# Environment Setup
# Add RepeatMasker to PATH (required by MAKER)
export PATH=$PATH:$REPEATMASKER_DIR

cd $OUT_DIR

# Load required modules
module load OpenMPI/4.1.1-GCC-10.3.0
module load AUGUSTUS/3.4.0-foss-2021a

# Run MAKER using MPI for parallel processing
mpiexec --oversubscribe -n 50 apptainer exec \
--bind $SCRATCH:/TMP --bind $COURSE_DIR --bind $AUGUSTUS_CONFIG_PATH --bind $REPEATMASKER_DIR --bind /data \
$MAKER_IMG \
maker -mpi --ignore_nfs_tmp -TMP /TMP maker_opts.ctl maker_bopts.ctl maker_evm.ctl maker_exe.ctl