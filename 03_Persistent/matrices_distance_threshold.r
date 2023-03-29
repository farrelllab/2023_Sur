#This script is to calculate sparse matrices showing neighbors for each tissue object

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
base.path <- "/data/CSD/zfext/LTA/results/13-EPS"
old.base.path <- "/data/CSD/zfext/results/04f-SubsetsV6"

# Set relative paths
old.mama.obj.path <- paste0("/data/CSD/zfext/results/02-Clustering/merged_timepoints/obj/merged_ALLhpf__igraph-weighted.rds")
new.mama.obj.path <- paste0("/data/CSD/zfext/LTA/results/02-Clustering/merged_timepoints/obj/merged_mama_ds_integrated_2021-01-11_seurat_UMAP_recalc.rds")

#base.path <- "~/Documents/zfext-biowulf/02-Clustering/"
##Assign paths
save.path <- paste0(base.path, "/knn_neighbors_data/")
plot.path <- paste0(base.path, "/plots_mat_data/")
dist.path <- paste0(old.base.path, "/dist_subsets_data/")
#seurat.path <- paste0(old.base.path, "/obj_subsets_v6/")

dist.path.data <- "/data/CSD/zfext/results/04f-SubsetsV6/dist_subsets_data/"
dist.path.scaled <- "/data/CSD/zfext/results/04f-SubsetsV6/dist_subsets_scaled/"
mat.path.data <- "/data/CSD/zfext/LTA/results/13-EPS/mat_subsets_data/"
mat.path.scaled <- "/data/CSD/zfext/LTA/results/13-EPS/mat_subsets_scaled/"
cells.path <- paste0(base.path, "/cells_subsets_v6/")


#message(paste0(Sys.time(), ": Creating subset for", sample))
#Load mama
#message(paste0(Sys.time(), ": Loading mama"))
#mama <- readRDS("/data/CSD/zfext/LTA/results/02-Clustering/merged_timepoints/obj/merged_mama_ds_integrated_2021-01-11_seurat_UMAP_recalc.rds")
#mama <- readRDS("/data/CSD/zfext/LTA/results/02-Clustering/merged_timepoints/obj/merged_mama_ds_slim_seurat_25Jan2022.rds")


##Load mama distance matrices 
##Load matrices
message(paste0(Sys.time(), ": Loading distance matrices for ", sample))
mat <- readRDS(paste0(dist.path.data, sample, "_dist_mat_data.rds"))

#Long term and short term categorization across tissues at Epsilon values of both 34 and 35 were quite similar at this stage, so we generated matrices for both these eps values.
message(paste0(Sys.time(), ": Iterate over many eps values for ", sample))
for(eps in c(34, 35)){

##Find distances that is below eps - those cells would be considered transcriptionally similar
message(paste0(Sys.time(), "Find elements in the matrix that is below the eps value"))
matches <- which(mat < eps, arr.ind = TRUE)
print(matches)

##keeping the same dimensions as the distance matrix, create a binary matrix to represent the most similar cells for each cell. 
message(paste0(Sys.time(), "Create a binary matrix of the same dimensions as the distance matrix"))
binary_matrix <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
colnames(binary_matrix) <- colnames(mat)
rownames(binary_matrix) <- rownames(mat)
binary_matrix[matches] <- 1

message(paste0(Sys.time(), "Save matrix"))
saveRDS(binary_matrix, paste0(save.path, sample, "_neighbor_matrix_", eps, "_binary.rds"))

#nearest <- sapply(rownames(mat), function(cell) mean(abs(obj@meta.data[setdiff(colnames(mat)[which(mat[cell,] <= eps)], cell), "stage.nice"] - obj@meta.data[cell, "stage.nice"])))
#message(paste0(Sys.time(), "Calculate neighbors as a list"))
#nn <- sapply(rownames(mat), function(cell) setdiff(colnames(mat)[which(mat[cell, ] <= eps)], cell))
#neighbors <- list(nearest, nn)
#names(neighbors) <- c("mean.stage.diff", "neighbors")

#message(paste0(Sys.time(), ": Saving"))
#saveRDS(nearest, file = paste0(save.path, sample, "_mean_stage_diff_eps_", eps, ".rds"))
#saveRDS(nn, file = paste0(save.path, sample, "_neighbors_eps_", eps, ".rds"))
}