library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

gff_file <- args[1]

gff <- read_tsv(
  gff_file,
  comment = "#",
  col_names = FALSE
)

colnames(gff) <- c("chr","src","type","start","end","score","strand","phase","attr")

# InterProScan features only
gff_ipr <- gff %>% filter(str_detect(attr, "Pfam|IPR|GO:"))

# Extract ID each annotation belongs to
parents <- str_extract(gff_ipr$attr, "(?<=Parent=)[^;]+")

# Extract IDs
gff$ID <- str_extract(gff$attr, "(?<=ID=)[^;]+")

# Pfam domains (IPR entries that map to Pfam)
pfam <- gff %>% filter(str_detect(attr, "Pfam"))
pfam_parents <- str_extract(pfam$attr, "(?<=Parent=)[^;]+")

# GO terms
go <- gff %>% filter(str_detect(attr, "GO:"))
go_parents   <- str_extract(go$attr,   "(?<=Parent=)[^;]+")

# Summary counts
cat("Proteins with any InterProScan annotation:", n_distinct(parents), "\n")
cat("Proteins with Pfam domains:", n_distinct(pfam_parents), "\n")
cat("Proteins with GO terms:", n_distinct(go_parents), "\n")

# Top 10 Pfam domains
pfam_counts <- pfam %>%
    mutate(domain = str_extract(attr, "Pfam:[^;]+")) %>%
    count(domain, sort = TRUE)
cat("Top 10 Pfam domains:\n")
print(head(pfam_counts, 10))

# Top 10 GO terms
go_counts <- go %>%
    mutate(go_term = str_extract(attr, "GO:[0-9]+")) %>%
    count(go_term, sort = TRUE)
cat("Top 10 GO terms:\n")
print(head(go_counts, 10))