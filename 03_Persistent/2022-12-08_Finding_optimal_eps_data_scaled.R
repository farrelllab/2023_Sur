library(Seurat)
library(Matrix)

#Load mama
mama <- readRDS("~/Box/zfext/02-Clustering/mama+ds_integrated/merged_mama_ds_slim_seurat_25Jan2022.rds")

# Load cell annotations to introduce for downsample purposes
cell.annot.path <- "~/Box/zfext/annotations_celltype_curated_newMama/celltype_annotations_df.tsv"
message(Sys.time(), ": Loading ", cell.annot.path)
cell.annot <- read.table(cell.annot.path, sep = " ", header = T, stringsAsFactors = F)
cell.annot$stage.group <- plyr::mapvalues(
  from = c("  3-4", "  5-6", "  7-9", " 10-12", " 14-21", " 24-34", " 36-46", 
           " 48-58", " 60-70", " 72-82", " 84-94", " 96-106", " 108-118", 
           "120", "5-6", "7-9", "10-12", "14-21", "24-34", "36-46", "48-58", 
           "60-70", "72-82", "84-94", "96-106", "108-118"),
  to =   c("  3-4", "  5-6", "  7-9", " 10-12", " 14-21", " 24-34", " 36-46", 
           " 48-58", " 60-70", " 72-82", " 84-94", " 96-106", " 108-118", 
           "120", "  5-6", "  7-9", " 10-12", " 14-21", " 24-34", " 36-46", " 48-58", 
           " 60-70", " 72-82", " 84-94", " 96-106", " 108-118"),
  x = cell.annot$stage.group
)
cell.annot$tissue <- unlist(lapply(strsplit(x = cell.annot$clust, split = "\\."), function(x) x[1]))
cell.annot$clust.sg <- paste0(cell.annot$clust, "_", gsub(" ", "", cell.annot$stage.group))
cell.annot$clust.stage <- paste0(cell.annot$clust, "_", gsub(" ", "", cell.annot$stage.nice))
cell.annot$tissue.sg <- paste0(cell.annot$tissue, "_", gsub(" ", "", cell.annot$stage.group))
cell.annot$tissue.stage <- paste0(cell.annot$tissue, "_", gsub(" ", "", cell.annot$stage.nice))

mama@meta.data <- cbind(
  mama@meta.data[,setdiff(colnames(mama@meta.data), colnames(cell.annot))],
  cell.annot[rownames(mama@meta.data),]
)

##Load distance matrix
dist.path.data <- "~/Box/zfext/global_analysis_curated/transient_vs_long-term/2022-11-30/dist_mat_data/"
dist.path.scaled <- "~/Box/zfext/global_analysis_curated/transient_vs_long-term/2022-11-30/dist_mat_scaled/"

##Read distance matrix
sample <- "pronephros"
dist.mat.data <- readRDS(paste0(dist.path.data, sample, "_dist_mat_data.rds"))
dist.mat.scaled <- readRDS(paste0(dist.path.scaled, sample, "_dist_mat_scaled.rds"))

##Use the precomputed distance matrix to find neighbors to query cells
##This function uses a cell id to find nearest neighbors within a certain user-defined radius distance
#' Calculate nearest neighbors within a specified epsilon
#' @param cell (a single cell id or several cell ids) 
#' @param eps (Character) a certain threshold radius
find_eps_neighbors <- function(cell, eps){
  neighbors <- setdiff(colnames(mat)[which(mat[cell, ] <= eps)], cell)
  stage  <- abs(mama@meta.data[neighbors, "stage.nice"] - mama@meta.data[cell, "stage.nice"])
  stage.mean <- mean(stage)
  nearest <- list(neighbors, stage)
  return(nearest)
}

##As an example calculate the nearest neighbor mean stage difference for one mamaect (use endoderm as example)
sample <- "fin"
dist.mat.data <- readRDS(paste0(dist.path.data, sample, "_dist_mat_data.rds"))
mat <- dist.mat.data

for(eps in c("20", "25", "30", "35", "40", "45", "50", "55", "60")){
nearest <- sapply(rownames(mat), function(cell) mean(abs(mama@meta.data[setdiff(colnames(mat)[which(mat[cell,] <= eps)], cell), "stage.nice"] - mama@meta.data[cell, "stage.nice"])))
nn <- sapply(rownames(mat), function(cell) setdiff(colnames(mat)[which(mat[cell, ] <= eps)], cell))
#neighbors <- list(nearest, nn)
#names(neighbors) <- c("mean.stage.diff", "neighbors")

saveRDS(nearest, paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/2022-11-30/knn_neighbors_data/", sample, "_data_stage_diff_eps_", eps, ".rds"))
saveRDS(nn, paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/2022-11-30/knn_neighbors_data/", sample, "_data_neighbors_eps_", eps, ".rds"))
#saveRDS(neighbors, paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/2022-11-30/knn_neighbors_data/", sample, "_data_neighbors_stage.diff.rds"))
}

rm(list = c("dist.mat.data", "dist.mat.scaled"))
shh <- gc()

##I did calculate stage differences for 
d <- dist(t(GetAssayData(mama, slot = "data")))

library(RcppArmadillo)
library(RcppXPtrUtils)
install.packages("parallelDist")
library(parallelDist)
install.packages("distances")
library(distances)
set.seed(123)

dist <- distances(t(mama@assays$RNA@data))
dist <- parDist(t(mama@assays$RNA@data), method = "euclidean")


  
