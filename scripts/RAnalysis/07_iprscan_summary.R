# =============================================================================
# Script: 07_iprscan_summary.R
# Purpose: Generate summary statistics for InterProScan annotations including
#          counts of annotated proteins, Pfam domains, and GO terms, plus
#          lists of the most common domains and GO terms.
# =============================================================================

# Load required libraries
library(tidyverse)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
gff_file <- args[1]  # Input: GFF file with InterProScan annotations

# Read GFF file
gff <- read_tsv(
  gff_file,
  comment = "#",
  col_names = FALSE
)

# Set column names
colnames(gff) <- c("chr", "src", "type", "start", "end", "score", "strand", "phase", "attr")

# Filter for InterProScan features only (Pfam, IPR, or GO terms)
gff_ipr <- gff %>% filter(str_detect(attr, "Pfam|IPR|GO:"))

# Extract parent IDs for InterProScan annotations
parents <- str_extract(gff_ipr$attr, "(?<=Parent=)[^;]+")

# Extract gene/protein IDs from GFF
gff$ID <- str_extract(gff$attr, "(?<=ID=)[^;]+")

# Extract Pfam domain annotations
pfam <- gff %>% filter(str_detect(attr, "Pfam"))
pfam_parents <- str_extract(pfam$attr, "(?<=Parent=)[^;]+")

# Extract GO term annotations
go <- gff %>% filter(str_detect(attr, "GO:"))
go_parents <- str_extract(go$attr, "(?<=Parent=)[^;]+")

# Print summary statistics
cat("Proteins with any InterProScan annotation:", n_distinct(parents), "\n")
cat("Proteins with Pfam domains:", n_distinct(pfam_parents), "\n")
cat("Proteins with GO terms:", n_distinct(go_parents), "\n")

# Calculate and print top 10 Pfam domains
pfam_counts <- pfam %>%
    mutate(domain = str_extract(attr, "Pfam:[^;]+")) %>%
    count(domain, sort = TRUE)
cat("Top 10 Pfam domains:\n")
print(head(pfam_counts, 10))

# Calculate and print top 10 GO terms
go_counts <- go %>%
    mutate(go_term = str_extract(attr, "GO:[0-9]+")) %>%
    count(go_term, sort = TRUE)
cat("Top 10 GO terms:\n")
print(head(go_counts, 10))