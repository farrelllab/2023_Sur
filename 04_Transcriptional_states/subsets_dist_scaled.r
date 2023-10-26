#This script is to find relationship between the distances calculated on gene expression space using an union of variable genes
##Calculate distances on normalized gene expression - using data slot

# To find: Relationship between distances for individual tissues calculated on two different PCAs

# Get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("Must provide at least one argument.")
sample = as.character(args[1])

# Load libraries
library(Seurat)
library(Matrix)
library(gridExtra)
library(tidyverse)

# Set base path
base.path <- "/data/CSD/zfext/results/04f-SubsetsV6"
old.base.path <- "/data/CSD/zfext/results/04c-SubsetsV3"


##Assign paths
dist.path.data <- "/data/CSD/zfext/results/04f-SubsetsV6/dist_subsets_data/"
dist.path.scaled <- "/data/CSD/zfext/results/04f-SubsetsV6/dist_subsets_scaled/"
mat.path.data <- "/data/CSD/zfext/LTA/results/13-EPS/mat_subsets_data/"
mat.path.scaled <- "/data/CSD/zfext/LTA/results/13-EPS/mat_subsets_scaled/"
#plot.path <- "/data/CSD/zfext/LTA/results/13-EPS/scatter_plots/"
#conv.path <- "/data/CSD/zfext/LTA/results/13-EPS/conversion/"
cells.path <- paste0(base.path, "/cells_subsets_v6/")
plot.path <- "/data/CSD/zfext/LTA/results/13-EPS/plots/"

##Load mama
message(paste0(Sys.time(), ": Loading mama"))
mama <- readRDS("/data/CSD/zfext/LTA/results/02-Clustering/merged_timepoints/obj/merged_mama_ds_integrated_2021-01-11_seurat_UMAP_recalc.rds")

##Load all variable genes
message(paste0(Sys.time(), ": List all variable genes"))
message(paste0(Sys.time(), ": Get the curated variable gene list"))
var.genes.all <- scan("/data/CSD/zfext/LTA/results/13-EPS/var_genes_all_tissues_filtered.txt", what = "character")

##Now calculate distance
##Get cells for individual tissue-specific subsets
##cells.list <- scan(paste0(cells.path, sample, "_cell_ids.txt"), what = "character")

##Load seurat object
obj <- readRDS(paste0(base.path, "/obj_subsets_v6/", sample, "_seurat.rds"))

##Get cells
cells.list <- WhichCells(obj)

# Check that all genes and cells are present in mama
missing.genes <- setdiff(var.genes.all, rownames(mama@assays$RNA@scale.data))
missing.cells <- setdiff(cells.list, colnames(mama@assays$RNA@scale.data))
genes.use <- intersect(var.genes.all, rownames(mama@assays$RNA@scale.data))
cells.use <- intersect(cells.list, colnames(mama@assays$RNA@scale.data))
message("Missing ", length(missing.genes), " genes: ", paste(missing.genes, sep = ", "))
message("Missing ", length(missing.cells), " cells: ", paste(missing.cells, sep = ", "))


##Calculate a distance matrix for all pairs of cells in this object. 
message(paste0(Sys.time(), ": Calculating distance matrix for ",  sample))
dist.mat <- dist(t(mama@assays$RNA@scale.data[genes.use, cells.use]))

#Save distance matrix
message(paste0(Sys.time(), ": Saving distance matrix for ",  sample))
saveRDS(dist.mat, paste0(dist.path.scaled, sample, "_mama_dist_scaled.rds"))

message(paste0(Sys.time(), ": Recovering RAM"))
rm(mama)
gc()

##Convert to matrix
#message(paste0("Converting distance matrix to matrix for",  sample))
#mat <- as.matrix(dist.mat)
#message(paste0("Saving converted matrix for",  sample))
#saveRDS(mat, paste0(mat.path.scaled, sample, "_mama_dist_scaled.rds"))

message(paste0(Sys.time(), ": Convert distance matrix to regular matrix for ", sample))
mat <- as.matrix(dist.mat)
message(paste0(Sys.time(), ": Saving distance matrix."))
saveRDS(mat, file = paste0(dist.path.scaled, sample, "_dist_mat_scaled.rds"))

# Plot and save distances between cells in this data to determine a reasonable downstream EPS.
message(Sys.time(), ": Histogram of cell-cell distances")
pdf(file = paste0(plot.path, "dist_all_scaled_", sample, ".pdf"), width = 8, height = 8)
hist(mat, breaks = 100, main = "All cell-cell distances", xlab = "Distance")
dev.off()

message(Sys.time(), ": Calculating KNN.")
nn <- apply(mat[,], 1, function(x) {
  y <- sort(x, decreasing = F)
  return(y[c(2,11,21,51)])
})
saveRDS(nn, file = paste0(dist.path.scaled, sample, "_nn_scaled.rds"))

message(Sys.time(), ": Plotting dists to NNs")
pdf(file = paste0(plot.path, "dist_nn_scaled_", sample, ".pdf"), width = 8, height = 8)
hist(nn[1,], breaks = 100, main = "Distance to 1st NN", xlab = "Distance")
hist(nn[2,], breaks = 100, main = "Distance to 10th NN", xlab = "Distance")
hist(nn[3,], breaks = 100, main = "Distance to 20th NN", xlab = "Distance")
hist(nn[4,], breaks = 100, main = "Distance to 50th NN", xlab = "Distance")
dev.off()
