#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=32G
#SBATCH --time=1:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=maker_aed_filter_and_subset
#SBATCH --output=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/maker_aed_filter_and_subset_output_%j.o
#SBATCH --error=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/maker_aed_filter_and_subset_error_%j.o

COURSE_DIR="/data/courses/assembly-annotation-course/CDS_annotation"
WORK_DIR="/data/users/aschmid/Euk-Genome-Organization-and-Annotation"

MAKER_BIN="${COURSE_DIR}/softwares/Maker_v3.01.03/src/bin"
PERL_MODULE="BioPerl/1.7.8-GCCcore-10.3.0"
UCSC_MODULE="UCSC-Utils/448-foss-2021a"
MARIADB_MODULE="MariaDB/10.6.4-GCC-10.3.0"
R_MODULE=R-bundle-IBU/2023072800-foss-2021a-R-4.2.1

# Input/Output files
GFF_IN="${WORK_DIR}/gene_annotation/final/assembly.all.maker.noseq.renamed.iprscan.gff"
GFF_FILTERED_OUT="${WORK_DIR}/gene_annotation/final/assembly.all.maker.noseq.iprscan_quality_filtered.gff"
GFF_GENIC_ONLY="${WORK_DIR}/gene_annotation/final/filtered.genes.renamed.gff"
MRNA_ID_LIST="${WORK_DIR}/gene_annotation/final/list.mRNA.ids.txt"
TRANSCRIPT_FASTA="${WORK_DIR}/gene_annotation/final/assembly.all.maker.transcripts.renamed.fasta"
PROTEIN_FASTA="${WORK_DIR}/gene_annotation/final/assembly.all.maker.proteins.renamed.fasta"
TRANSCRIPT_FASTA_OUT="${WORK_DIR}/gene_annotation/final/assembly.all.maker.transcripts.renamed.filtered.fasta"
PROTEIN_FASTA_OUT="${WORK_DIR}/gene_annotation/final/assembly.all.maker.proteins.renamed.filtered.fasta"

module load $PERL_MODULE
module load $UCSC_MODULE
module load $MARIADB_MODULE
module load $R_MODULE 

# Filter by AED / InterProScan
perl "${MAKER_BIN}/quality_filter.pl" -s "${GFF_IN}" > "${GFF_FILTERED_OUT}"

# Keep only gene-related features
grep -P "\tgene\t|\tCDS\t|\texon\t|\tfive_prime_UTR\t|\tthree_prime_UTR\t|\tmRNA\t" "${GFF_FILTERED_OUT}" > "${GFF_GENIC_ONLY}"
  
# Check
cut -f3 "${GFF_GENIC_ONLY}" | sort | uniq -c

# Generate AED cumulative distribution file
perl "${MAKER_BIN}/AED_cdf_generator.pl" -b 0.025 "${GFF_GENIC_ONLY}" > "${GFF_GENIC_ONLY}.AED.txt"

# Make list of mRNA IDs (column 9: “ID=…”)
grep -P '\tmRNA\t' "${GFF_GENIC_ONLY}" \
    | awk '{print $9}' \
    | cut -d ';' -f1 \
    | sed 's/ID=//g' \
    > "${MRNA_ID_LIST}"

# Subset transcripts & proteins (faSomeRecords)
faSomeRecords "${TRANSCRIPT_FASTA}" "${MRNA_ID_LIST}" "${TRANSCRIPT_FASTA_OUT}"
faSomeRecords "${PROTEIN_FASTA}"    "${MRNA_ID_LIST}" "${PROTEIN_FASTA_OUT}"

Rscript $WORK_DIR/scripts/RAnalysis/07_iprscan_summary.R $GFF_FILTERED_OUT 

echo "Done."