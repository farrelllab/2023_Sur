#!/bin/Rscript

# mamasmooth_fuzzyc_cluster.R
# Jeff Farrell - 04-21-2022
# Can fuzzy clustering work on all mama?

# Arguments:
#       1: c (Number of clusters to generate)
#       2: m (Fuzziness parameter for Mfuzz: 1 is hard clustering, >1 allows multiple membership)

library(Seurat)
library(Mfuzz)
library(Matrix)
#library(heatmaply)

source("/data/CSD/zfext/LTA/scripts/11-Modules/knn_smooth.R")
source("/data/CSD/zfext/LTA/scripts/11-Modules/2022-04-13 Gene Module Functions.R")

## Get arguments -----------------------------

message(paste0(Sys.time(), ": Starting."))
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) stop("Please provide 2 arguments. (See R script for details.)")

n <- as.numeric(args[1])
m <- as.numeric(args[2])

# Load data --------------

# Load mama with smoothed data
mama.path <- paste0("/data/CSD/zfext/LTA/results/12-KNNSmooth/merged_mama_ds_slim_seurat_25Jan2022_smoothedk5.rds")
message(Sys.time(), ": Loading", mama.path)
mama <- readRDS(mama.path)

# Load cell annotations to introduce for downsample purposes
cell.annot <- read.table("/data/CSD/zfext/LTA/scripts/12-KNNSmooth/mama_celltype_annotations_df.tsv", sep = " ", header = T, stringsAsFactors = F)
cell.annot$stage.group <- plyr::mapvalues(
  from = c("  3-4", "  5-6", "  7-9", " 10-12", " 14-21", " 24-34", " 36-46", 
           " 48-58", " 60-70", " 72-82", " 84-94", " 96-106", " 108-118", 
           "120", "5-6", "7-9", "10-12", "14-21", "24-34", "36-46", "48-58", 
           "60-70", "72-82", "84-94", "96-106", "108-118", "4-6", " 120", " 7-9"),
  to =   c("  3-4", "  5-6", "  7-9", " 10-12", " 14-21", " 24-34", " 36-46", 
           " 48-58", " 60-70", " 72-82", " 84-94", " 96-106", " 108-118", 
           "120", "  5-6", "  7-9", " 10-12", " 14-21", " 24-34", " 36-46", " 48-58", 
           " 60-70", " 72-82", " 84-94", " 96-106", " 108-118", "  5-6", "120", "  7-9"),
  x = cell.annot$stage.group
)
cell.annot$clust.sg <- paste0(cell.annot$clust, "_", gsub(" ", "", cell.annot$stage.group))
mama@meta.data <- cbind(
  mama@meta.data[,setdiff(colnames(mama@meta.data), colnames(cell.annot))],
  cell.annot[rownames(mama@meta.data),]
)

# Calculate modules ------------------

# Genes to operate on: all variable genes used in the object, plus the best markers of each cluster
# This ends up being too many genes, unfortunately, so going to have to choose your faves.
message(Sys.time(), ": Determining genes to use.")
n.cells <- ncol(mama)
row.chunks <- chunk(1:nrow(mama), 50)
prop.exp <- unlist(lapply(row.chunks, function(chunk) {
  return(rowSums(mama@assays$RNA@data[chunk,] > 0.1) / n.cells)
}))
names(prop.exp) <- rownames(mama)
# Going to consider genes expressed in at least 0.1% of cells and fewer than 75%
prop.exp <- prop.exp * 100
genes.use <- intersect(names(which(prop.exp >= 0.1 & prop.exp < 75)), rownames(mama))
message(Sys.time(), ": Going to use ", length(genes.use), " of ", nrow(mama), " genes in mama.")

# Going to downsample to 50% of cells at random? Because the data is too big!
# Add PCA back
#mama2 <- readRDS("/data/CSD/zfext/LTA/results/02-Clustering/merged_timepoints/obj/merged_mama_ds_slim.rds")
#mama@reductions$pca <- mama2@reductions$pca
#rm(mama2); gc()
# Create KNN
#mama <- FindNeighbors(mama, reduction = "pca", dims = 1:100, k.param = 10, compute.SNN = T)
#rownames(mama@neighbors$RNA.nn@nn.idx) <- mama@neighbors$RNA.nn@cell.names

#half.cells <- seq(1, ncol(mama), 2)

# Going to downsample:
# Min 20% of each stage-cluster
# Min 50 cells per cluster
message(Sys.time(), ": Determining cells to use")

clust.n <- table(cell.annot$clust.sg)
clust.n.mod <- unlist(lapply(clust.n, function(n) {
  if (n < 50) return(n)
  return(max(ceiling(n * .20), 50))
}))
cells.use <- unlist(lapply(names(clust.n.mod), function(clust) {
  n <- clust.n.mod[[clust]]
  cells <- rownames(cell.annot)[cell.annot$clust.sg == clust]
  cells.do <- cells[ceiling(1:n * (length(cells) / n))]
  return(cells.do)
}))
cells.use <- intersect(cells.use, colnames(mama))

message(Sys.time(), ": Going to use ", length(cells.use), " of ", ncol(mama), " cells in mama.")

exp.data <- ExpressionSet(as.matrix(
  mama@assays$RNA@data[intersect(genes.use, rownames(mama)), intersect(cells.use, colnames(mama))]
))
exp.data.std <- Mfuzz::standardise(exp.data)

# Generate c-means fuzzy clustering
message(Sys.time(), ": Beginning clustering.")
cl <- mfuzz(exp.data.std, c=n, m=m, verbose = T, iter.max = 50)

message(Sys.time(), ": Saving clustering as a backup.")      
saveRDS(cl, file = paste0("/data/CSD/zfext/LTA/results/11-Modules/cl_mama_c", n, "_m", m, ".rds"))

# Filter clusters
# Add clustering result to Seurat
mama <- fuzzy.clust.to.seurat(mama, cl = cl, min.a = 0.2, min.gene.per.mod = 5, module.combine.thresh = 0.95, do.plot = F, add.to.meta = T)

# Do plots!
file.prefix <- paste0("c", n, "_m", m)
plot.path <- paste0("/data/CSD/zfext/LTA/results/11-Modules/fc_mama_plots/", file.prefix, "/")

dir.create(plot.path, recursive = T)
write.csv2(mama@reductions$fc@feature.loadings, file=paste0(plot.path, file.prefix, "_feature_loadings.csv"))
saveRDS(mama@reductions$fc, file = paste0(plot.path, file.prefix, "_mamaFCreduc.rds"))
plot.all.modules.mamaonly(obj = obj, mama = mama, plot.path = plot.path)






