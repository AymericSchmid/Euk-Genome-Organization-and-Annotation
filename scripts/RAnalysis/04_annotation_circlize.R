# Load the circlize package
library(circlize)
library(tidyverse)
library(ComplexHeatmap)

# Load the TE annotation GFF3 file
args = commandArgs(trailingOnly=TRUE)

assembly_fai_file <- args[1]
te_gff_file <- args[2]
te_gff_data <- read.table(te_gff_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
gene_gff_file <- ""
gene_gff_data <- NULL
out_pdf_file <- ""

if (length(args) == 3) {
    out_pdf_file <- args[3]
} else if (length(args) == 4) {
    gene_gff_file <- args[3]
    out_pdf_file <- args[4]
    gene_gff_data <- read.table(gene_gff_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
} else {
    stop("Insufficient arguments provided. Please provide at least assembly_fai_file, te_gff_file, and out_pdf_file.")
}

# Check the superfamilies present in the GFF3 file, and their counts
superfam_counts <- table(te_gff_data$V3)
print(superfam_counts)
superfams_sorted = rev(names(sort(superfam_counts)))

# custom ideogram data
## To make the ideogram data, you need to know the lengths of the scaffolds.
## There is an index file that has the lengths of the scaffolds, the `.fai` file.
## To generate this file you need to run the following command in bash:
## samtools faidx assembly.fasta
## This will generate a file named assembly.fasta.fai
## You can then read this file in R and prepare the custom ideogram data

custom_ideogram <- read.table(assembly_fai_file, header = FALSE, stringsAsFactors = FALSE)
custom_ideogram$chr <- custom_ideogram$V1
custom_ideogram$start <- 1
custom_ideogram$end <- custom_ideogram$V2
custom_ideogram <- custom_ideogram[, c("chr", "start", "end")]
custom_ideogram <- custom_ideogram[order(custom_ideogram$end, decreasing = T), ]
sum(custom_ideogram$end[1:20])

# Select only the first 20 longest scaffolds, You can reduce this number if you have longer chromosome scale scaffolds
custom_ideogram <- custom_ideogram[1:20, ]

# Function to filter GFF3 data based on Superfamily (You need one track per Superfamily)
filter_superfamily <- function(te_gff_data, superfamily, custom_ideogram) {
    filtered_data <- te_gff_data[te_gff_data$V3 == superfamily, ] %>%
        as.data.frame() %>%
        mutate(chrom = V1, start = V4, end = V5, strand = V6) %>%
        select(chrom, start, end, strand) %>%
        filter(chrom %in% custom_ideogram$chr)
    return(filtered_data)
}

filter_clade <- function(te_gff_data, clade, custom_ideogram) {
    filtered_data <- te_gff_data %>%
        mutate(clade_extracted = str_match(V9, "(?<=;clade=)[^;]+")[,1]) %>%
        filter(!is.na(clade_extracted), clade_extracted %in% clade, V1 %in% custom_ideogram$chr) %>%
        transmute(chrom = V1, start = as.integer(V4), end   = as.integer(V5))
    return(filtered_data)
}

filter_genes <- function(gene_gff_data, custom_ideogram) {
    filtered_data <- gene_gff_data[gene_gff_data$V3 == "gene", ] %>%
        filter(V1 %in% custom_ideogram$chr) %>%
        transmute(chrom = V1, start = as.integer(V4), end = as.integer(V5))
    return(filtered_data)
}

pdf(out_pdf_file, width = 10, height = 10)
gaps <- c(rep(1, length(custom_ideogram$chr) - 1), 5) # Add a gap between scaffolds, more gap for the last scaffold
circos.par(start.degree = 90, gap.after = 1, track.margin = c(0, 0), gap.degree = gaps)
# Initialize the circos plot with the custom ideogram
circos.genomicInitialize(custom_ideogram)

# Plot te density
circos.genomicDensity(filter_superfamily(te_gff_data, "Gypsy_LTR_retrotransposon", custom_ideogram), count_by = "number", col = "darkgreen", track.height = 0.07, window.size = 1e5)
circos.genomicDensity(filter_superfamily(te_gff_data, "Copia_LTR_retrotransposon", custom_ideogram), count_by = "number", col = "darkred", track.height = 0.07, window.size = 1e5)
# Also plotting the top 2 most abundant superfamilies
circos.genomicDensity(filter_superfamily(te_gff_data, superfams_sorted[1], custom_ideogram), count_by = "number", col = "darkblue", track.height = 0.07, window.size = 1e5)
circos.genomicDensity(filter_superfamily(te_gff_data, superfams_sorted[2], custom_ideogram), count_by = "number", col = "darkorange", track.height = 0.07, window.size = 1e5)
# Plot Athila and CRM clades
circos.genomicDensity(filter_clade(te_gff_data, "Athila", custom_ideogram), count_by = "number", col = "red", track.height = 0.07, window.size = 1e5)
circos.genomicDensity(filter_clade(te_gff_data, "CRM", custom_ideogram), count_by = "number", col = "blue", track.height = 0.07, window.size = 1e5)
# Plot gene density if gene annotation is provided
if (!is.null(gene_gff_data)) {
    circos.genomicDensity(filter_genes(gene_gff_data, custom_ideogram), count_by = "number", col = "black", track.height = 0.07, window.size = 1e5)
}
circos.clear()

legend_labels <- c("Gypsy_LTR_retrotransposon", "Copia_LTR_retrotransposon", superfams_sorted[1], superfams_sorted[2], "Athila", "CRM")
legend_colors <- c("darkgreen", "darkred", "darkblue", "darkorange", "red", "blue")
if (!is.null(gene_gff_data)) {
    legend_labels <- c(legend_labels, "Genes")
    legend_colors <- c(legend_colors, "black")
}

lgd <- Legend(
    title = "Superfamily", at = legend_labels,
    legend_gp = gpar(fill = legend_colors)
)
draw(lgd, x = unit(8, "cm"), y = unit(10, "cm"), just = c("center"))

dev.off()

# Plot the distribution of Athila and CRM clades (known centromeric TEs in Brassicaceae).
# You need to run the TEsorter on TElib to get the clades classification from the TE library