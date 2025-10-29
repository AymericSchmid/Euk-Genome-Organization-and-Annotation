# Load the circlize package
library(circlize)
library(tidyverse)
library(ComplexHeatmap)

# Load the TE annotation GFF3 file
args = commandArgs(trailingOnly=TRUE)
gff_file <- args[1]
gff_data <- read.table(gff_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
assembly_fai_file <- args[2]

# Check the superfamilies present in the GFF3 file, and their counts
superfam_counts <- table(gff_data$V3)
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
filter_superfamily <- function(gff_data, superfamily, custom_ideogram) {
    filtered_data <- gff_data[gff_data$V3 == superfamily, ] %>%
        as.data.frame() %>%
        mutate(chrom = V1, start = V4, end = V5, strand = V6) %>%
        select(chrom, start, end, strand) %>%
        filter(chrom %in% custom_ideogram$chr)
    return(filtered_data)
}

filter_clade <- function(gff_data, clade, custom_ideogram) {
    print(gff_data$V9)
    filtered_data <- gff_data[gff_data$V9 == clade, ] %>%
        as.data.frame() %>%
        mutate(chrom = V1, start = V4, end = V5, strand = V6) %>%
        select(chrom, start, end, strand) %>%
        filter(chrom %in% custom_ideogram$chr)
    return(filtered_data)
}

pdf("04_TE_density.pdf", width = 10, height = 10)
gaps <- c(rep(1, length(custom_ideogram$chr) - 1), 5) # Add a gap between scaffolds, more gap for the last scaffold
circos.par(start.degree = 90, gap.after = 1, track.margin = c(0, 0), gap.degree = gaps)
# Initialize the circos plot with the custom ideogram
circos.genomicInitialize(custom_ideogram)

# Plot te density
circos.genomicDensity(filter_superfamily(gff_data, "Gypsy_LTR_retrotransposon", custom_ideogram), count_by = "number", col = "darkgreen", track.height = 0.07, window.size = 1e5)
circos.genomicDensity(filter_superfamily(gff_data, "Copia_LTR_retrotransposon", custom_ideogram), count_by = "number", col = "darkred", track.height = 0.07, window.size = 1e5)
# Also plotting the top 2 most abundant superfamilies
circos.genomicDensity(filter_superfamily(gff_data, superfams_sorted[1], custom_ideogram), count_by = "number", col = "darkblue", track.height = 0.07, window.size = 1e5)
circos.genomicDensity(filter_superfamily(gff_data, superfams_sorted[2], custom_ideogram), count_by = "number", col = "darkorange", track.height = 0.07, window.size = 1e5)
circos.genomicDensity(filter_clade(gff_data, "Athila", custom_ideogram), count_by = "number", col = "red", track.height = 0.07, window.size = 1e5)
circos.genomicDensity(filter_clade(gff_data, "CRM", custom_ideogram), count_by = "number", col = "green", track.height = 0.07, window.size = 1e5)
circos.clear()

lgd <- Legend(
    title = "Superfamily", at = c("Gypsy_LTR_retrotransposon", "Copia_LTR_retrotransposon", superfams_sorted[1], superfams_sorted[2], "Athila", "CRM"),
    legend_gp = gpar(fill = c("darkgreen", "darkred", "darkblue", "darkorange", "red", "green"))
)
draw(lgd, x = unit(8, "cm"), y = unit(10, "cm"), just = c("center"))

dev.off()

# Plot the distribution of Athila and CRM clades (known centromeric TEs in Brassicaceae).
# You need to run the TEsorter on TElib to get the clades classification from the TE library