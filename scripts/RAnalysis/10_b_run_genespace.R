# =============================================================================
# Script: 10_b_run_genespace.R
# Purpose: Run GENESPACE for comparative genomics analysis to identify syntenic
#          orthologs across multiple genomes and generate pangenome data.
#          Saves pangenome matrix for downstream analysis.
# =============================================================================

# Load required libraries
library(GENESPACE)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]           # Working directory with prepared BED and FASTA files
mcscanx_path <- args[2]  # Path to MCScanX executable for synteny detection

# Initialize GENESPACE parameters
gpar <- init_genespace(wd = wd, path2mcscanx = mcscanx_path, nCores = 1)

# Run GENESPACE analysis
out <- run_genespace(gpar)

# Query pangenome data
# bed = NULL: use all genes
# refGenome = "TAIR10": use TAIR10 as reference
# transform = TRUE: transform coordinates
# showArrayMem/showNSOrtho: display memory and non-syntenic ortholog information
pangenome <- query_pangenes(
  out, 
  bed = NULL, 
  refGenome = "TAIR10", 
  transform = TRUE, 
  showArrayMem = TRUE, 
  showNSOrtho = TRUE, 
  maxMem2Show = Inf
)

# Save pangenome object as RDS file for downstream analysis
saveRDS(pangenome, file = file.path(wd, "pangenome_matrix.rds"))

# Note: In downstream scripts, load the pangenome matrix with:
# pangenome <- readRDS(file.path(wd, "pangenome_matrix.rds"))
# Then use it for analyses such as calculating core, accessory, and species-specific genes