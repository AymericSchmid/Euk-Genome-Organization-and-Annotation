library(tidyverse)
library(cowplot)

args <- commandArgs(trailingOnly = TRUE)

aed_file <- args[1]
plot_file <- args[2]

# Read AED cumulative distribution file
aed <- read.table(aed_file, header = TRUE)
colnames(aed) <- c("AED", "cumulative_fraction")

# Recover individual AED values
aed <- aed %>%
    arrange(AED) %>%
    mutate(freq = c(cumulative_fraction[1], diff(cumulative_fraction))) # freq = fraction of genes in that AED bin

# AED Histogram
p1 <- ggplot(aed, aes(x = AED, y = freq)) +
    geom_col(fill = "#1b9e77", color = "black", width = 0.02) +
    theme_cowplot() +
    labs(title = "AED Score Distribution (Histogram)",
        x = "Annotation Edit Distance (AED)",
        y = "Proportion of genes")

# AED Cumulative curve
p2 <- ggplot(aed, aes(x = AED, y = cumulative_fraction)) +
    geom_line(color = "#d95f02", size = 1) +
    geom_point(size = 1.5) +
    theme_cowplot() +
    labs(title = "Cumulative AED Distribution",
        x = "Annotation Edit Distance (AED)",
        y = "Cumulative proportion of genes")

combined <- plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1))
ggsave(plot_file, combined, width = 12, height = 6)