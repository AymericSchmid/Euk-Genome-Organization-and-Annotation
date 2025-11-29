# =============================================================================
# Script: 07_genes_AED_distribution.R
# Purpose: Visualize the distribution of Annotation Edit Distance (AED) scores
#          for gene annotations. Creates both histogram and cumulative distribution
#          plots to assess annotation quality (lower AED = better quality).
# =============================================================================

# Load required libraries
library(tidyverse)
library(cowplot)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
aed_file <- args[1]   # Input: AED cumulative distribution file
plot_file <- args[2]  # Output: PDF filename for plots

# Read AED cumulative distribution file
aed <- read.table(aed_file, header = TRUE)
colnames(aed) <- c("AED", "cumulative_fraction")

# Recover individual AED bin frequencies from cumulative distribution
# freq = fraction of genes in each AED bin
aed <- aed %>%
    arrange(AED) %>%
    mutate(freq = c(cumulative_fraction[1], diff(cumulative_fraction)))

# Create histogram of AED distribution
p1 <- ggplot(aed, aes(x = AED, y = freq)) +
    geom_col(fill = "#1b9e77", color = "black", width = 0.02) +
    theme_cowplot() +
    labs(
        title = "AED Score Distribution (Histogram)",
        x = "Annotation Edit Distance (AED)",
        y = "Proportion of genes"
    )

# Create cumulative distribution curve
p2 <- ggplot(aed, aes(x = AED, y = cumulative_fraction)) +
    geom_line(color = "#d95f02", size = 1) +
    geom_point(size = 1.5) +
    theme_cowplot() +
    labs(
        title = "Cumulative AED Distribution",
        x = "Annotation Edit Distance (AED)",
        y = "Cumulative proportion of genes"
    )

# Combine plots side-by-side
combined <- plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1))

# Save combined plot
ggsave(plot_file, combined, width = 12, height = 6)