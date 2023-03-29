#This script is to find relationship between the distances calculated on gene expression space using an union of variable genes calculated on each tissue specific dataset
##Calculate distances on normalized gene expression - using data slot

# To find: distances between cells of individual tissue subsets in gene expression space

# Get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("Must provide at least one argument.")
sample = as.character(args[1])

# Load libraries
library(Seurat)
library(Matrix)
library(gridExtra)

# Set base path
base.path <- "/data/CSD/zfext/results/04f-SubsetsV6"
old.base.path <- "/data/CSD/zfext/results/04c-SubsetsV3"


##Assign paths
dist.path.data <- "/data/CSD/zfext/results/04f-SubsetsV6/dist_subsets_data/"
dist.path.scaled <- "/data/CSD/zfext/results/04f-SubsetsV6/dist_subsets_scaled/"
mat.path.data <- "/data/CSD/zfext/LTA/results/13-EPS/mat_subsets_data/"
mat.path.scaled <- "/data/CSD/zfext/LTA/results/13-EPS/mat_subsets_scaled/"
cells.path <- paste0(base.path, "/cells_subsets_v6/")
plot.path <- "/data/CSD/zfext/LTA/results/13-EPS/plots/"

##Load global dataset
message(paste0(Sys.time(), ": Loading global dataset"))
mama <- readRDS("/data/CSD/zfext/LTA/results/02-Clustering/merged_timepoints/obj/merged_mama_ds_integrated_2021-01-11_seurat_UMAP_recalc.rds")

##Load all variable genes
message(paste0(Sys.time(), ": List all variable genes"))
paths <- paste0("/data/CSD/zfext/results/04f-SubsetsV6/var_subsets_v6/")
var.genes <- lapply(list.files(paths, pattern = "_vargenes.txt", full.names = T), scan, what = "character")

message(paste0(Sys.time(), ": Take a union of all variable genes"))
var.genes.all <- unlist(unique(var.genes))

##Load tissue specific seurat object
message(paste0(Sys.time(), ": Load seurat object"))
obj <- readRDS(paste0(base.path, "/obj_subsets_v6/", sample, "_seurat.rds"))

##Now calculate distance
##Get cells for individual tissue-specific subsets
#cells.list <- scan(paste0(cells.path, sample, "_cell_ids.txt"), what = "character")

message(paste0(Sys.time(), ": get the list of cells"))
cells.list <- WhichCells(obj)

# Check that all genes and cells are present in mama
missing.genes <- setdiff(var.genes.all, rownames(mama@assays$RNA@data))
missing.cells <- setdiff(cells.list, colnames(mama@assays$RNA@data))
genes.use <- intersect(var.genes.all, rownames(mama@assays$RNA@data))
cells.use <- intersect(cells.list, colnames(mama@assays$RNA@data))
message("Missing ", length(missing.genes), " genes: ", paste(missing.genes, sep = ", "))
message("Missing ", length(missing.cells), " cells: ", paste(missing.cells, sep = ", "))

message(paste0(Sys.time(), ": Calculate distance matrix"))
dist.mat <- dist(t(as.matrix(mama@assays$RNA@data[genes.use, cells.use])))

#Save distance matrix
message(paste0(Sys.time(), ": Saving distance matrix for ",  sample))
saveRDS(dist.mat, paste0(dist.path.data, sample, "_mama_dist_data.rds"))

message(paste0(Sys.time(), ": Recovering RAM"))
rm(mama)
gc()

##Convert to matrix
#message(paste0("Converting distance matrix to matrix for",  sample))
#mat <- as.matrix(dist.mat)
#message(paste0("Saving converted matrix for",  sample))
#saveRDS(mat, paste0(mat.path.data, sample, "_mama_dist_data.rds"))
#message(paste0(Sys.time(), ": Loading distance matrix for ", sample))
#mat <- readRDS(file = paste0(dist.path, sample, "_mama_dist_data.rds"))

##Convert the distance object into a matrix
message(paste0(Sys.time(), ": Convert distance matrix to regular matrix for ", sample))
mat <- as.matrix(dist.mat)
message(paste0(Sys.time(), ": Saving distance matrix."))
saveRDS(mat, file = paste0(dist.path.data, sample, "_dist_mat_data.rds"))

# Plot and save distances between cells in this data to determine a reasonable downstream EPS (i.e., constant distance representing transcriptional similarity).
message(Sys.time(), ": Histogram of cell-cell distances")
pdf(file = paste0(plot.path, "dist_all_data_", sample, ".pdf"), width = 8, height = 8)
hist(mat, breaks = 100, main = "All cell-cell distances", xlab = "Distance")
dev.off()

message(Sys.time(), ": Calculating KNN.")
nn <- apply(mat[,], 1, function(x) {
  y <- sort(x, decreasing = F)
  return(y[c(2,11,21,51)])
})
saveRDS(nn, file = paste0(dist.path.data, sample, "_nn_data.rds"))

message(Sys.time(), ": Plotting dists to NNs")
pdf(file = paste0(plot.path, "dist_nn_data_", sample, ".pdf"), width = 8, height = 8)
hist(nn[1,], breaks = 100, main = "Distance to 1st NN", xlab = "Distance")
hist(nn[2,], breaks = 100, main = "Distance to 10th NN", xlab = "Distance")
hist(nn[3,], breaks = 100, main = "Distance to 20th NN", xlab = "Distance")
hist(nn[4,], breaks = 100, main = "Distance to 50th NN", xlab = "Distance")
dev.off()
