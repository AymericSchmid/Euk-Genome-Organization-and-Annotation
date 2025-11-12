library(GENESPACE)

args <- commandArgs(trailingOnly = TRUE)

# get the folder where the genespace workingDirectory is located
wd <- args[1]
mcscanx_path <- args[2]  # path to MCScanX executable

gpar <- init_genespace(wd = wd, path2mcscanx = mcscanx_path, nCores = 1)

# run genespace
out <- run_genespace(gpar)
pangenome <- query_pangenes(out, bed = NULL, refGenome = "TAIR10", transform = TRUE, showArrayMem = TRUE, showNSOrtho = TRUE, maxMem2Show = Inf)

# save pangenome object as rds
saveRDS(pangenome, file = file.path(wd, "pangenome_matrix.rds"))

# in your next script, you can load the pangenome matrix with:
# pangenome <- readRDS(file.path(wd, "pangenome_matrix.rds"))
# and then use it for downstream analyses, e.g., calculating core, accessory and specific genes