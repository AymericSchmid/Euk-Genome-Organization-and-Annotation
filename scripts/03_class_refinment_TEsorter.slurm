#!/usr/bin/env bash
# =============================================================================
# Script: 03_class_refinment_TEsorter.slurm
# Purpose: Refine TE classification by extracting Copia and Gypsy sequences,
#          classifying them into clades using TEsorter, and appending clade
#          information to the TE annotation GFF file
# =============================================================================

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=1:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=class_refinment_TEsorter
#SBATCH --output=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/class_refinment_TEsorter_output_%j.o
#SBATCH --error=/data/users/aschmid/Euk-Genome-Organization-and-Annotation/output/class_refinment_TEsorter_error_%j.o

WORK_DIR=/data/users/aschmid/Euk-Genome-Organization-and-Annotation
TE_LIB=$WORK_DIR/TE_annotation/EDTA/assembly.fa.mod.EDTA.TElib.fa
TE_ANNOTATION=$WORK_DIR/TE_annotation/EDTA/assembly.fa.mod.EDTA.TEanno.gff3
TE_ANNOTATION_CLADES=$WORK_DIR/TE_annotation/TEanno.gff3
OUT_DIR=$WORK_DIR/TE_annotation
COPIA_SEQ_FILE=$OUT_DIR/Copia_sequences.fa
GYPSY_SEQ_FILE=$OUT_DIR/Gypsy_sequences.fa

# Expected TEsorter output file names
COPIA_CLS="${COPIA_SEQ_FILE##*/}.rexdb-plant.cls.tsv"
GYPSY_CLS="${GYPSY_SEQ_FILE##*/}.rexdb-plant.cls.tsv"

# Software modules and containers
SEQKIT_MODULE=SeqKit/2.6.1
TESORTER_IMG=/data/courses/assembly-annotation-course/CDS_annotation/containers/TEsorter_1.3.0.sif

mkdir -p $OUT_DIR

module load $SEQKIT_MODULE

# Extract all sequences containing "Copia" in their name or description
seqkit grep -r -p "Copia" $TE_LIB > $COPIA_SEQ_FILE

# Extract all sequences containing "Gypsy" in their name or description
seqkit grep -r -p "Gypsy" $TE_LIB > $GYPSY_SEQ_FILE

cd $OUT_DIR

# Classify Copia sequences into clades
apptainer exec --bind /data $TESORTER_IMG \
TEsorter $COPIA_SEQ_FILE -db rexdb-plant

# Classify Gypsy sequences into clades
apptainer exec --bind /data $TESORTER_IMG \
TEsorter $GYPSY_SEQ_FILE -db rexdb-plant

# Merge Classification Results
CLS_ALL="$OUT_DIR/te_clades.cls.tsv"
# Combine header from Copia file with data from both files (skip header from Gypsy)
{
  head -n 1 "$COPIA_CLS"
  tail -n +2 "$COPIA_CLS"
  tail -n +2 "$GYPSY_CLS"
} > "$CLS_ALL"

# Append Clade Information to TE Annotation GFF
# Use awk to read the classification table and append clade information
# to the TE annotation GFF file attributes column
awk -F'\t' -v OFS='\t' -v CLS="$CLS_ALL" '
BEGIN{
  # Read TEsorter classification table into memory
  while ((getline line < CLS) > 0) {
    if (line ~ /^[[:space:]]*$/) continue
    split(line, f, "\t")
    if (f[1] ~ /TE_/) {
      key = f[1]
      sub(/#.*/, "", key)          # TE_00000100_INT#LTR/Gypsy -> TE_00000100_INT
      ORD[key] = f[2]              # Order
      SF[key]  = f[3]              # Superfamily
      CLA[key] = (f[4] != "" ? f[4] : "")  # Clade (may be empty)
    }
  }
  close(CLS)
}
/^#/ { print; next }

{
  attrs=$9
  teid=""
  # Find the TE id in attributes (EDTA uses Name=..., sometimes Target=...)
  if (match(attrs, /Name=([^;]+)/, m))      teid = m[1]
  else if (match(attrs, /Target=([^;]+)/, m)) teid = m[1]

  add = ";clade="
  if (teid != "" && (teid in SF)) {
    # Append clade annotation, preserving existing attributes
    add = add CLA[teid]
  }
  attrs = attrs add
  $9 = attrs
  print
}
' "$TE_ANNOTATION" > "$TE_ANNOTATION_CLADES"