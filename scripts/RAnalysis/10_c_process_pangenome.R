# =============================================================================
# Script: 10_c_process_pangenome.R
# Purpose: Process pangenome data from GENESPACE to calculate core, accessory,
#          and species-specific genes. Generates frequency plots showing the
#          distribution of orthogroups and genes across genomes.
# =============================================================================

# Load required libraries
library(tidyverse)
library(readr)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]           # Working directory where pangenome_matrix.rds is located
focal_genome <- args[2] # Focal accession name for analysis
out_plot <- args[3]     # Output plot file path

# Load pangenome matrix from GENESPACE output
pangenome <- readRDS(file.path(wd, "pangenome_matrix.rds"))

# Identify genome columns (list-columns produced by query_pangenes())
genome_cols <- names(pangenome)[sapply(pangenome, is.list)]

# Convert to tibble for easier manipulation
pg <- as_tibble(pangenome)

# Function to clean gene lists: remove out-of-synteny genes (IDs ending with '*')
clean_gene_list <- function(v) {
  if (is.null(v) || length(v) == 0) {
    return(character(0))
  }
  v <- as.character(v)
  v <- v[!is.na(v)]
  v <- trimws(v)
  v <- v[!grepl("\\*$", v)]  # Drop genes ending with '*' (out-of-synteny)
  unique(v)
}

# Clean gene lists in all genome columns
pg <- pg %>%
  mutate(across(all_of(genome_cols), ~ lapply(.x, clean_gene_list)))

# Remove orthogroups that are now empty in ALL genomes after cleaning
pg <- pg %>%
  rowwise() %>%
  filter(
    sum(sapply(c_across(all_of(genome_cols)), length)) > 0
  ) %>%
  ungroup()

# Create presence/absence table (orthogroup-level)
# TRUE if the cell has a non-empty list of gene IDs
presence_tbl <- pg %>%
  transmute(
    pgID,
    across(
      all_of(genome_cols),
      ~ lengths(.x) > 0,
      .names = "{.col}"
    )
  )

# Calculate global orthogroup flags
n_genomes <- length(genome_cols)

pg_flags <- presence_tbl %>%
  mutate(
    # Count how many genomes have each orthogroup
    n_present = select(., all_of(genome_cols)) %>% rowSums(),

    # Core = present in ALL genomes
    is_core_all = (n_present == n_genomes),

    # Accessory = present in SOME genomes (not all)
    is_accessory = (n_present < n_genomes),

    # Is this orthogroup in our focal genome?
    present_in_focal = .data[[focal_genome]],

    # Focal-specific = ONLY in focal genome
    focal_specific = (present_in_focal & n_present == 1),

    # Lost in focal = present in other genomes but NOT focal
    lost_in_focal = (!present_in_focal & n_present > 0),

    # Almost-core lost in focal = in all genomes EXCEPT focal
    lost_core_without_focal = (!present_in_focal & n_present == (n_genomes - 1)),

    # Identify genes lost compared to TAIR10
    lost_vs_TAIR10 = (!present_in_focal) & .data[["TAIR10"]]
  ) %>%
  # Add simple category labels
  mutate(
    category = case_when(
      is_core_all ~ "core",
      n_present == 1 ~ "species_specific",
      TRUE ~ "accessory"
    )
  )

# Calculate gene counts per orthogroup × genome
# Convert list entries to integer counts (number of unique gene IDs)
count_genes <- function(gene_list) {
  if (is.null(gene_list) || length(gene_list) == 0) {
    return(0)
  }
  # Count unique genes
  length(unique(gene_list))
}

# Apply count function to all genome columns
gene_counts_matrix <- pg %>%
  select(pgID, all_of(genome_cols)) %>%
  mutate(across(
    all_of(genome_cols),
    ~ sapply(.x, count_genes)
  ))

# Calculate gene counts per genome and category
gene_counts_w_cat <- pg_flags %>%
  select(pgID, category) %>%
  left_join(gene_counts_matrix, by = "pgID")

# Calculate total genes per genome (sum gene counts by category)
gene_by_cat <- gene_counts_w_cat %>%
  pivot_longer(cols = all_of(genome_cols), names_to = "genome", values_to = "gene_count") %>%
  group_by(genome, category) %>%
  summarise(gene_count = sum(gene_count), .groups = "drop")

# Total genes per genome (sum across all categories)
gene_totals <- gene_by_cat %>%
  group_by(genome) %>%
  summarise(gene_total = sum(gene_count), .groups = "drop")

# Create per-genome summary table with categories
gene_counts_per_genome <- gene_by_cat %>%
  pivot_wider(names_from = category, values_from = gene_count, values_fill = 0) %>%
  rename(
    gene_core = core,
    gene_accessory = accessory,
    gene_specific = species_specific
  ) %>%
  left_join(gene_totals, by = "genome") %>%
  mutate(
    # Calculate percentages of genes by category
    percent_core = round(100 * gene_core / pmax(gene_total, 1), 2),
    percent_specific = round(100 * gene_specific / pmax(gene_total, 1), 2)
  )


# Create pangenome frequency plot: orthogroups and genes present in n genomes
# Count orthogroups by number of genomes they're present in
og_freq <- pg_flags %>%
  count(n_present, name = "count") %>%
  mutate(type = "Orthogroups")

# Create long table of (pgID, genome, gene) and attach orthogroup-level n_present
all_genes_with_presence <- pg %>%
  select(pgID, all_of(genome_cols)) %>%
  pivot_longer(cols = all_of(genome_cols), names_to = "genome", values_to = "genes_list") %>%
  # Drop empty/null lists
  filter(!sapply(genes_list, is.null) & sapply(genes_list, length) > 0) %>%
  unnest_longer(genes_list) %>%
  rename(gene = genes_list) %>%
  filter(!is.na(gene)) %>%
  mutate(gene = as.character(gene)) %>%
  distinct(pgID, genome, gene) %>%  # One row per distinct gene occurrence in a pgID/genome
  left_join(select(pg_flags, pgID, n_present), by = "pgID")

# Count genes by the n_present of their orthogroup
# This counts unique gene × orthogroup occurrences
gene_freq <- all_genes_with_presence %>%
  distinct(pgID, gene, n_present) %>%  # Unique gene in each orthogroup
  count(n_present, name = "count") %>%
  mutate(type = "Genes")

# Combine orthogroup and gene frequency data
freq_data <- bind_rows(og_freq, gene_freq)

# Create frequency plot
ggplot(freq_data, aes(x = n_present, y = count, fill = type)) +
  geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
  scale_x_continuous(
    breaks = 1:n_genomes,
    labels = 1:n_genomes
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = c("Genes" = "#0072B2", "Orthogroups" = "#D55E00")) +
  labs(
    x = "Number of genomes",
    y = "Count",
    fill = NULL,
    title = "Pangenome composition: distribution across genomes",
    subtitle = paste("Total genomes:", n_genomes)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# Save plot
ggsave(file.path(out_plot), width = 10, height = 6)

# Create and print summary table
freq_summary <- freq_data %>%
  pivot_wider(names_from = type, values_from = count, values_fill = 0) %>%
  mutate(
    label = case_when(
      n_present == n_genomes ~ "Core (all genomes)",
      n_present == 1 ~ "Species-specific (1 genome)",
      TRUE ~ paste0("Shared (", n_present, " genomes)")
    )
  ) %>%
  arrange(desc(n_present)) %>%
  select(n_present, label, Orthogroups, Genes)

print(freq_summary)
