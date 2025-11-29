# =============================================================================
# Script: 05_plot_div.R
# Purpose: Generate TE divergence landscape plot showing the distribution of
#          transposable element sequences by divergence percentage from consensus.
#          Plots sequence amount (Mbp) rather than counts to account for
#          fragmented TEs due to nested insertions and deletions.
# =============================================================================

# Load required libraries
library(reshape2)
library(tidyverse)
library(data.table)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
data <- args[1]  # Input file: divergence landscape table from parseRM.pl

# Read divergence landscape data
rep_table <- fread(data, header = FALSE, sep = "\t")

# Set column names: Rname, Rclass, Rfam, and divergence bins (1-50%)
colnames(rep_table) <- c("Rname", "Rclass", "Rfam", 1:50)

# Filter out unknown families
rep_table <- rep_table %>% filter(Rfam != "unknown")

# Create combined family label (Class/Family)
rep_table$fam <- paste(rep_table$Rclass, rep_table$Rfam, sep = "/")

# Reshape data from wide to long format
rep_table.m <- melt(rep_table)

# Remove peak at divergence = 1% (library sequences are copies in the genome,
# which inflate this low divergence peak)
rep_table.m <- rep_table.m[-c(which(rep_table.m$variable == 1)), ]

# Set factor levels for proper ordering in plot
# Order: LTR/Copia, LTR/Gypsy, DNA transposons, Helitron, MITEs
rep_table.m$fam <- factor(rep_table.m$fam, levels = c(
  "LTR/Copia", "LTR/Gypsy", 
  "DNA/DTA", "DNA/DTC", "DNA/DTH", "DNA/DTM", "DNA/DTT", "DNA/Helitron",
  "MITE/DTA", "MITE/DTC", "MITE/DTH", "MITE/DTM"
))

# Convert divergence from percentage to decimal (variable is percent divergence)
rep_table.m$distance <- as.numeric(rep_table.m$variable) / 100

# Remove Helitrons as EDTA cannot annotate them properly
# Reference: https://github.com/oushujun/EDTA/wiki/Making-sense-of-EDTA-usage-and-outputs---Q&A
rep_table.m <- rep_table.m %>% filter(fam != "DNA/Helitron")

# Disable implicit Rplots.pdf device
pdf(NULL)

# Create divergence landscape plot
# Weight by sequence amount (Mbp) rather than counts to account for fragmented TEs
ggplot(rep_table.m, aes(fill = fam, x = distance, weight = value / 1000000)) +
  geom_bar() +
  cowplot::theme_cowplot() +
  scale_fill_brewer(palette = "Paired") +
  xlab("Distance") +
  ylab("Sequence (Mbp)") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 1, size = 9, hjust = 1), 
    plot.title = element_text(hjust = 0.5)
  )

# Save plot
ggsave(filename = "05_TE_divergence_landscape.pdf", width = 10, height = 5, useDingbats = FALSE)

# Note: Why plot in Mbp instead of counts?
# Consider a scenario with many small TE fragments due to nested insertions and deletions.
# Using counts would over-represent these fragments, while using Mbp (sequence amount)
# provides a more accurate representation of the actual TE content in the genome.
