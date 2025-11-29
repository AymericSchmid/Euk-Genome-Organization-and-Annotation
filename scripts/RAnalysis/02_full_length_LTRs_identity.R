# =============================================================================
# Script: 02_full_length_LTRs_identity.R
# Purpose: Analyze the identity of full-length LTR retrotransposons (LTR-RTs)
#          annotated in a genome. Reads GFF file with LTR-RT annotations,
#          extracts LTR identity values, merges with TEsorter classification data,
#          and generates plots showing the distribution of LTR identities per
#          clade within the Copia and Gypsy superfamilies.
# =============================================================================

# Load required libraries
library(tidyverse)
library(data.table)
library(cowplot)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
gff_file <- args[1]  # GFF file with LTR-RT annotations
cls_file <- args[2]  # TEsorter classification TSV file

# Read and preprocess LTR-RT annotation GFF file
message("Reading GFF: ", gff_file)
anno <- read.table(gff_file, sep = "\t", header = FALSE)

# Remove subfeatures (Terminal repeats, TSDs) to keep only top-level TE annotations
exclude_feats <- c("long_terminal_repeat", "repeat_region", "target_site_duplication")
anno <- anno %>% filter(!V3 %in% exclude_feats)

# Extract Name and ltr_identity from attributes column (V9) using regex
anno <- anno %>%
  rowwise() %>%
  mutate(
    # Extract Name=... from attributes (V9)
    Name = str_extract(V9, "(?<=Name=)[^;]+"),
    # Extract ltr_identity=... value
    Identity = as.numeric(str_extract(V9, "(?<=ltr_identity=)[^;]+")),
    # Compute element length as end - start
    length = as.numeric(V5) - as.numeric(V4)
  ) %>%
  # Keep only columns needed for downstream analysis
  select(V1, V4, V5, V3, Name, Identity, length)

# Read TEsorter classification table
message("Reading classification: ", cls_file)
cls <- fread(cls_file, sep = "\t", header = TRUE)
setnames(cls, 1, "TE")

# Process TEsorter classification data
# TEsorter encodes internal domain classification as TEName_INT#Classification
# Split on '#' to separate name and classification, keep only internal-domain
# matches (Name ends with _INT), and strip _INT suffix
cls <- cls %>%
  separate(TE, into = c("Name", "Classification"), sep = "#", fill = "right") %>%
  filter(str_detect(Name, "_INT")) %>%
  mutate(Name = str_remove(Name, "_INT$"))

# Merge annotation with classification table
# Use left join to keep all annotated TEs even if they have no classification match
anno_cls <- merge(anno, cls, by = "Name", all.x = TRUE)

# Quality control: check data within plotting range
binwidth <- 0.005  # Bin width for histogram
xlims <- c(0.80, 1.00)  # X-axis limits for identity values
outside <- anno_cls %>%
  dplyr::filter(!is.na(Identity)) %>%
  summarise(
    below = sum(Identity < xlims[1]),
    above = sum(Identity > xlims[2]),
    total = dplyr::n()
  )
message(sprintf("Identities outside plotting range: <%.2f = %d, >%.2f = %d (%.1f%% of %d).",
                xlims[1], outside$below, xlims[2], outside$above,
                100*(outside$below+outside$above)/outside$total, outside$total))

# Calculate clade statistics for ordering facets and showing sample sizes
clade_stats <- anno_cls %>%
  filter(Superfamily %in% c("Copia", "Gypsy"), !is.na(Identity), !is.na(Clade)) %>%
  group_by(Superfamily, Clade) %>%
  summarise(
    n = sum(between(Identity, xlims[1], xlims[2])), 
    med_id = median(Identity[between(Identity, xlims[1], xlims[2])], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(clade_lab = paste0(Clade, " (n=", n, ")"))

# Attach clade labels and order facets by median identity
anno_cls <- anno_cls %>%
  left_join(clade_stats, by = c("Superfamily", "Clade")) %>%
  group_by(Superfamily) %>%
  mutate(clade_lab = forcats::fct_reorder(clade_lab, med_id, .desc = TRUE)) %>%
  ungroup()

# Quick summary: counts per Superfamily and Clade
# (may be NA if classification is missing)
message("Counts per Superfamily")
print(table(anno_cls$Superfamily, useNA = "ifany"))
message("Counts per Clade")
print(table(anno_cls$Clade, useNA = "ifany"))


# Calculate global y-axis maximum for consistent scaling across plots
# Bin identity values into consistent breaks and find maximum count
global_ymax <- anno_cls %>%
  filter(Superfamily %in% c("Copia", "Gypsy"), !is.na(Identity)) %>%
  count(Superfamily, Clade, Identity = cut(Identity, seq(xlims[1], xlims[2], by = binwidth))) %>%
  pull(n) %>%
  max(na.rm = TRUE)

# Define plotting function for one superfamily
plot_by_clade <- function(df, sf, ymax) {
  df %>%
    dplyr::filter(
      Superfamily == sf, 
      !is.na(Identity), 
      !is.na(Clade), 
      Identity >= xlims[1], 
      Identity <= xlims[2]
    ) %>%
    ggplot(aes(x = Identity)) +
    geom_histogram(
      binwidth = binwidth, 
      boundary = xlims[1], 
      closed = "right", 
      color = "black",
      fill = ifelse(sf == "Copia", "#E69F00", "#009E73")
    ) +
    # Add vertical lines at 95% and 98% identity thresholds
    geom_vline(xintercept = c(0.95, 0.98), linetype = "dashed", linewidth = 0.3) +
    facet_wrap(~clade_lab, ncol = 1, scales = "fixed") +
    scale_x_continuous(limits = xlims, breaks = seq(xlims[1], xlims[2], 0.05)) +
    scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
    theme_cowplot() +
    theme(
      strip.background = element_rect(fill = "#f0f0f0"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(title = sf, x = "Identity", y = "Count")
}

# Generate plots for Copia and Gypsy superfamilies
p_copia <- plot_by_clade(anno_cls, "Copia", global_ymax)
p_gypsy <- plot_by_clade(anno_cls, "Gypsy", global_ymax)

# Combine plots side-by-side using cowplot
combined <- plot_grid(p_copia, p_gypsy, ncol = 2, rel_widths = c(1, 1))

# Save plots in both PNG and PDF formats
ggsave("02_LTR_Copia_Gypsy_cladelevel.png", combined, width = 12, height = 10, dpi = 300)
ggsave("02_LTR_Copia_Gypsy_cladelevel.pdf", combined, width = 12, height = 10)