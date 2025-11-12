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

# ---- counts & ordering ----
superfam_counts  <- table(te_gff_data$V3)
print(superfam_counts)
superfams_sorted <- rev(names(sort(superfam_counts)))

# avoid duplicating Copia/Gypsy as "top2"
base_sf <- c("Gypsy_LTR_retrotransposon","Copia_LTR_retrotransposon", "helitron")
top2    <- setdiff(superfams_sorted, base_sf)[1:2]

# ---- ideogram (top 14 scaffolds) ----
custom_ideogram <- read.table(assembly_fai_file, header = FALSE, stringsAsFactors = FALSE)
custom_ideogram <- custom_ideogram[, c(1,2)]
colnames(custom_ideogram) <- c("chr","len")
custom_ideogram <- custom_ideogram |>
  arrange(desc(len)) |>
  mutate(start = 1, end = len) |>
  select(chr, start, end) |>
  slice(1:14)

# ---- helpers ----
filter_superfamily <- function(te_gff_data, superfamily, custom_ideogram) {
  te_gff_data[te_gff_data$V3 == superfamily, ] |>
    as.data.frame() |>
    mutate(chrom = V1, start = as.integer(V4), end = as.integer(V5)) |>
    select(chrom, start, end) |>
    filter(chrom %in% custom_ideogram$chr)
}

filter_clade <- function(te_gff_data, clade, custom_ideogram) {
  te_gff_data |>
    mutate(clade_extracted = stringr::str_match(V9, "(?<=;clade=)[^;]+")[,1]) |>
    filter(!is.na(clade_extracted), clade_extracted %in% clade, V1 %in% custom_ideogram$chr) |>
    transmute(chrom = V1, start = as.integer(V4), end = as.integer(V5))
}

filter_genes <- function(gene_gff_data, custom_ideogram) {
  gene_gff_data[gene_gff_data$V3 == "gene", ] |>
    filter(V1 %in% custom_ideogram$chr) |>
    transmute(chrom = V1, start = as.integer(V4), end = as.integer(V5))
}

ref_sector <- custom_ideogram$chr[1] 
track_height <- 0.1

okabe <- c(
  black = "#000000", orange = "#E69F00", sky = "#56B4E9", blue = "#0072B2",
  vermilion = "#D55E00", purple = "#CC79A7", green = "#009E73",
  yellow = "#F0E442", grey = "#999999"
)

col_gypsy  <- "#009E73"     
col_copia  <- "#E69F00"  
col_top1   <- "#F0E442"       
col_top2   <- "#B5B5B5"
col_athila <- "#CC79A7"    
col_crm    <- "#56B4E9"
col_gene   <- "#333333"

density_track <- function(df, col, ws = 1e5, h = 0.085,
                          label = NULL, ref_sector = NULL,
                          offset_mm = 0.9, cex_lab = 0.6,
                          bending = TRUE) {
  if (nrow(df) == 0) return(invisible(NULL))
  circos.genomicDensity(
    df, count_by = "number", window.size = ws,
    col = col, border = TRUE, lwd = 0.5, baseline = "bottom",
    track.height = h
  )
  if (!is.null(label) && !is.null(ref_sector)) {
    # place a single label on the chosen sector, inside the ring
    ti   <- get.current.track.index()
    xlim <- get.cell.meta.data("xlim", sector.index = ref_sector, track.index = ti)
    ylim <- get.cell.meta.data("ylim", sector.index = ref_sector, track.index = ti)
    xmid <- mean(xlim)
    ytop <- ylim[2]
    y <- ytop - convert_y(offset_mm, "mm",
                          sector.index = ref_sector, track.index = ti)
    circos.text(xmid, y, labels = label,
                sector.index = ref_sector, track.index = ti,
                facing = if (bending) "bending.inside" else "inside",
                niceFacing = TRUE, cex = cex_lab, adj = c(0.5, 1))
  }
}

pdf(out_pdf_file, width = 10, height = 10)
# ---- plotting ----

gaps <- c(rep(1, nrow(custom_ideogram)))
circos.clear()
circos.par(
  start.degree = 90,
  gap.after = gaps,
  cell.padding = c(0.002, 0.002, 0.002, 0.002),
  track.margin = c(0.002, 0.002),
  points.overflow.warning = FALSE
)
circos.genomicInitialize(custom_ideogram)

# choose a reference scaffold to place ring labels (the longest)
ref_sector <- custom_ideogram$chr[1]

# tune window size a bit by genome scale (larger genome -> larger window)
ws <- if (max(custom_ideogram$end) > 2e7) 2e5 else 1e5

# ---- rings (outer -> inner) ----
density_track(filter_superfamily(te_gff_data, "Gypsy_LTR_retrotransposon", custom_ideogram),
              col_gypsy, ws, h = 0.09, label = "Gypsy",  ref_sector = ref_sector)

density_track(filter_superfamily(te_gff_data, "Copia_LTR_retrotransposon", custom_ideogram),
              col_copia, ws, h = 0.09, label = "Copia",  ref_sector = ref_sector)

# top1 / top2 (exclude Gypsy/Copia already handled)
if (length(top2) >= 1) {
  density_track(filter_superfamily(te_gff_data, top2[1], custom_ideogram),
                col_top1, ws, h = 0.09, label = top2[1], ref_sector = ref_sector)
}
if (length(top2) >= 2) {
  density_track(filter_superfamily(te_gff_data, top2[2], custom_ideogram),
                col_top2, ws, h = 0.09, label = top2[2], ref_sector = ref_sector)
}

density_track(filter_clade(te_gff_data, "Athila", custom_ideogram),
              col_athila, ws, h = 0.09, label = "Athila", ref_sector = ref_sector)

density_track(filter_clade(te_gff_data, "CRM", custom_ideogram),
              col_crm, ws, h = 0.09, label = "CRM", ref_sector = ref_sector)

if (!is.null(gene_gff_data)) {
  density_track(filter_genes(gene_gff_data, custom_ideogram),
                col_gene, ws, h = 0.09, label = "Genes", ref_sector = ref_sector)
}

# caption
grid::grid.text("Tracks (outer -> inner): Gypsy, Copia, Top1, Top2, Athila, CRM, Genes",
                x = grid::unit(0.5, "npc"), y = grid::unit(0.02, "npc"),
                just = c("center", "bottom"),
                gp = grid::gpar(cex = 0.85, fontface = "bold"))

# ---- grouped legend (optional but recommended) ----
lg_superfam <- ComplexHeatmap::Legend(title = "Superfamilies",
                                      labels = c("Gypsy", "Copia"),
                                      legend_gp = grid::gpar(fill = c(col_gypsy, col_copia))
)
lg_topfeat <- ComplexHeatmap::Legend(title = "Top feature types",
                                     labels = c(top2[1], top2[2]),
                                     legend_gp = grid::gpar(fill = c(col_top1, col_top2))
)
lg_clades <- ComplexHeatmap::Legend(title = "Clades",
                                    labels = c("Athila", "CRM"),
                                    legend_gp = grid::gpar(fill = c(col_athila, col_crm))
)
lg_genes <- ComplexHeatmap::Legend(title = "Genes",
                                   labels = "Gene density",
                                   legend_gp = grid::gpar(fill = col_gene)
)
lgd <- ComplexHeatmap::packLegend(lg_superfam, lg_topfeat, lg_clades, lg_genes,
                                  column_gap = grid::unit(3, "mm"))
grid::grid.text("Tracks", x = grid::unit(0.98, "npc"), y = grid::unit(0.98, "npc"),
                just = c("right", "top"), gp = grid::gpar(fontface = "bold"))
ComplexHeatmap::draw(lgd, x = grid::unit(0.98, "npc"), y = grid::unit(0.94, "npc"),
                     just = c("right", "top"))

# finish
circos.clear()


dev.off()

# Plot the distribution of Athila and CRM clades (known centromeric TEs in Brassicaceae).
# You need to run the TEsorter on TElib to get the clades classification from the TE library