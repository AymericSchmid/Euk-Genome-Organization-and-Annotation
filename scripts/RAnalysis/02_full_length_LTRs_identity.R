# In this script, we analyze the identity of full-length LTR retrotransposons (LTR-RTs)
# annotated in a genome. We read a GFF file with LTR-RT annotations,
# extract the LTR identity values, merge with classification data from TEsorter,
# and generate plots showing the distribution of LTR identities per clade within
# the Copia and Gypsy superfamilies.
library(tidyverse)
library(data.table)
library(cowplot)

#-------------------------------------------------
# Input files (edit paths if needed)
#-------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
gff_file <- args[1]
cls_file <- args[2]

#-------------------------------------------------
# Read and preprocess input data
#-------------------------------------------------
message("Reading GFF: ", gff_file)
anno <- read.table(gff_file, sep = "\t", header = FALSE)

# Remove subfeatures (Terminal repeats, TSDs) so we keep top-level TE annotations
exclude_feats <- c("long_terminal_repeat", "repeat_region", "target_site_duplication")
anno <- anno %>% filter(!V3 %in% exclude_feats)

# Extract Name and ltr_identity from the ninth column (attributes). This uses regex.

anno <- anno %>%
  rowwise() %>%
  mutate(
    # extract Name=... from attributes (V9)
    Name = str_extract(V9, "(?<=Name=)[^;]+"),
    # extract ltr_identity=... 
    Identity = as.numeric(str_extract(V9, "(?<=ltr_identity=)[^;]+")),
    # compute length as end - start
    length = as.numeric(V5) - as.numeric(V4)
  ) %>%
  # keep only the columns used downstream
  select(V1, V4, V5, V3, Name, Identity, length)

message("Reading classification: ", cls_file)
# Read classification table (TE name in first column). If your file doesn't have a header, set header=FALSE.
cls <- fread(cls_file, sep = "\t", header = TRUE)
setnames(cls, 1, "TE")

# TEsorter outputs encode the internal domain classification as TEName_INT#Classification. We split on '#',
# then keep only rows that correspond to internal-domain matches (Name ends with _INT), and strip _INT.
cls <- cls %>%
  separate(TE, into = c("Name", "Classification"), sep = "#", fill = "right") %>%
  filter(str_detect(Name, "_INT")) %>%
  mutate(Name = str_remove(Name, "_INT$"))

## Merge annotation with classification table
# Use a left join so all annotated TEs are kept even if they have no classification match
anno_cls <- merge(anno, cls, by = "Name", all.x = TRUE)

# Small QC on x-limits
binwidth <- 0.005
xlims <- c(0.80, 1.00)
outside <- anno_cls %>%
  dplyr::filter(!is.na(Identity)) %>%
  summarise(below = sum(Identity < xlims[1]),
            above = sum(Identity > xlims[2]),
            total = dplyr::n())
message(sprintf("Identities outside plotting range: <%.2f = %d, >%.2f = %d (%.1f%% of %d).",
                xlims[1], outside$below, xlims[2], outside$above,
                100*(outside$below+outside$above)/outside$total, outside$total))

# clade stats to order facets and show n
inrange <- between(anno_cls$Identity, xlims[1], xlims[2])
clade_stats <- anno_cls %>%
  filter(Superfamily %in% c("Copia","Gypsy"), !is.na(Identity), !is.na(Clade)) %>%
  group_by(Superfamily, Clade) %>%
  summarise(n = sum(between(Identity, xlims[1], xlims[2])), 
            med_id = median(Identity[between(Identity, xlims[1], xlims[2])], na.rm=TRUE),
            .groups="drop") %>%
  mutate(clade_lab = paste0(Clade, " (n=", n, ")"))


# attach the label + facet order back to the main table
anno_cls <- anno_cls %>%
  left_join(clade_stats, by = c("Superfamily","Clade")) %>%
  group_by(Superfamily) %>%
  mutate(clade_lab = forcats::fct_reorder(clade_lab, med_id, .desc=TRUE)) %>%
  ungroup()


# Quick checks: how many per Superfamily/Clade (may be NA if classification missing)
message("Counts per Superfamily")
print(table(anno_cls$Superfamily, useNA = "ifany"))
message("Counts per Clade")
print(table(anno_cls$Clade, useNA = "ifany"))


#-------------------------------------------------
# Plot setup
#-------------------------------------------------
# change global_ymax to per-superfamily ymax
global_ymax <- anno_cls %>%
  filter(Superfamily %in% c("Copia", "Gypsy"), !is.na(Identity)) %>%
  # bin Identity into consistent breaks and count occurrences
  count(Superfamily, Clade, Identity = cut(Identity, seq(xlims[1], xlims[2], by = binwidth))) %>%
  pull(n) %>%
  max(na.rm = TRUE)


#-------------------------------------------------
# Plot function for one superfamily
#-------------------------------------------------
plot_by_clade <- function(df, sf, ymax) {
  df %>%
    dplyr::filter(Superfamily == sf, !is.na(Identity), !is.na(Clade), Identity >= xlims[1], Identity <= xlims[2]) %>%
    ggplot(aes(x = Identity)) +
    geom_histogram(binwidth = binwidth, boundary = xlims[1], closed = "right", color = "black",
                   fill = ifelse(sf == "Copia", "#E69F00", "#009E73")) +
    geom_vline(xintercept = c(0.95, 0.98), linetype = "dashed", linewidth = 0.3) +
    facet_wrap(~clade_lab, ncol = 1, scales = "fixed") +   # use the new label
    scale_x_continuous(limits = xlims, breaks = seq(xlims[1], xlims[2], 0.05)) +
    scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
    theme_cowplot() +
    theme(strip.background = element_rect(fill = "#f0f0f0"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5)) +
    labs(title = sf, x = "Identity", y = "Count")
}

#-------------------------------------------------
# Generate Copia and Gypsy plots
#-------------------------------------------------
p_copia <- plot_by_clade(anno_cls, "Copia", global_ymax)
p_gypsy <- plot_by_clade(anno_cls, "Gypsy", global_ymax)

# Combine with cowplot side-by-side 
combined <- plot_grid(p_copia, p_gypsy, ncol = 2, rel_widths = c(1, 1))

ggsave("02_LTR_Copia_Gypsy_cladelevel.png", combined, width = 12, height = 10, dpi = 300)
ggsave("02_LTR_Copia_Gypsy_cladelevel.pdf", combined, width = 12, height = 10)