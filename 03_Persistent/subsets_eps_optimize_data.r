#This script is to optimize a radius or distance in the global dataset that would be considered an universal constant ('epsilon') across all tissues.
# To find: timings of differentiation of various cell types across tissues, presence of long-term cycling states in tissues
#at various eps radius values and check where the neighbors lie

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

#Load slim version of global dataset
message(paste0(Sys.time(), ": Loading mama"))
#mama <- readRDS("/data/CSD/zfext/LTA/results/02-Clustering/merged_timepoints/obj/merged_mama_ds_integrated_2021-01-11_seurat_UMAP_recalc.rds")
mama <- readRDS("/data/CSD/zfext/LTA/results/02-Clustering/merged_timepoints/obj/merged_mama_ds_slim_seurat_25Jan2022.rds")

##Add annotations
cell.annot <- read.table("/data/CSD/zfext/LTA/results/02-Clustering/merged_timepoints/meta.data/celltype_annotations_df.tsv", sep = " ", header = T, stringsAsFactors = F)
cell.annot$tissue <- unlist(lapply(strsplit(x = cell.annot$clust, split = "\\."), function(x) x[1]))
cell.annot$clust.sg <- paste0(cell.annot$clust, "_", gsub(" ", "", cell.annot$stage.group))

mama@meta.data <- cbind(
  mama@meta.data[,setdiff(colnames(mama@meta.data), colnames(cell.annot))],
  cell.annot[rownames(mama@meta.data),]
)

##Recovering RAM
obj <- mama
rm(list = c("mama"))
shh <- gc()


##Load global distance matrices 
##Load matrices
message(paste0(Sys.time(), ": Loading distance matrices for ", sample))
mat <- readRDS(paste0(dist.path.data, sample, "_dist_mat_data.rds"))

#Iterate over many epsilon values to find optimal transcriptomic similarity 
message(paste0(Sys.time(), ": Iterate over many eps values for ", sample))
for(eps in c(20, 25, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 50, 55, 60)){
nearest <- sapply(rownames(mat), function(cell) mean(abs(obj@meta.data[setdiff(colnames(mat)[which(mat[cell,] <= eps)], cell), "stage.nice"] - obj@meta.data[cell, "stage.nice"])))
nn <- sapply(rownames(mat), function(cell) setdiff(colnames(mat)[which(mat[cell, ] <= eps)], cell))
neighbors <- list(nearest, nn)
names(neighbors) <- c("mean.stage.diff", "neighbors")

message(paste0(Sys.time(), ": Saving"))
saveRDS(nearest, file = paste0(save.path, sample, "_mean_stage_diff_eps_", eps, ".rds"))
saveRDS(nn, file = paste0(save.path, sample, "_neighbors_eps_", eps, ".rds"))

##Long-term states would be all those states where the stage differences are atleast more than 3 days
message(paste0(Sys.time(), ": Get long-term and short-term states at 24h threshold"))
nearest.long.term.24 <- names(nearest)[which(nearest >= 24)]
##Short term states would be states that last less than 8 hours
nearest.short.term.24 <- names(nearest)[which(nearest < 24)]

message(paste0(Sys.time(), ": Get long-term and short-term states at 36h threshold"))
nearest.long.term.36 <- names(nearest)[which(nearest >= 36)]
##Short term states would be states that last less than 8 hours
nearest.short.term.36 <- names(nearest)[which(nearest < 36)]

message(paste0(Sys.time(), ": Get long-term and short-term states at 48h threshold"))
nearest.long.term.48 <- names(nearest)[which(nearest >= 48)]
##Short term states would be states that last less than 8 hours
nearest.short.term.48 <- names(nearest)[which(nearest < 48)]

#Calculate proportion of long-term states per cluster
#message(paste0(Sys.time(), ": get proportion of long and short-term states per clusters"))
#long.term.clust.prop <- table(obj@meta.data[nearest.long.term, "clust"])/table(obj@meta.data$clust)
#short.term.clust.prop <- table(obj@meta.data[nearest.short.term, "clust"])/table(obj@meta.data$clust)

#message(paste0(Sys.time(), ": Convert proportions into a dataframe"))
#long.prop <- as.data.frame(long.term.clust.prop)
#colnames(long.prop) <- c("celltype", "proportion")

##Create a slot in the object and save these cells under a name
message(paste0(Sys.time(), ": Create a slot in object for these states for 24h"))
obj@meta.data$cell.states.24 <- NA
obj@meta.data[nearest.long.term.24, "cell.states.24"] <- "long-term"
obj@meta.data[nearest.short.term.24, "cell.states.24"] <- "short-term"

##Create a slot in the object and save these cells under a name
message(paste0(Sys.time(), ": Create a slot in object for these states for 36h"))
obj@meta.data$cell.states.36 <- NA
obj@meta.data[nearest.long.term.36, "cell.states.36"] <- "long-term"
obj@meta.data[nearest.short.term.36, "cell.states.36"] <- "short-term"

##Create a slot in the object and save these cells under a name
message(paste0(Sys.time(), ": Create a slot in object for these states for 36h"))
obj@meta.data$cell.states.48 <- NA
obj@meta.data[nearest.long.term.48, "cell.states.48"] <- "long-term"
obj@meta.data[nearest.short.term.48, "cell.states.48"] <- "short-term"


##Plotting
nst <- as.data.frame(nearest)
colnames(nst) <- c("stg.diff")
message(paste0(Sys.time(), ": Plotting graphs for ", sample))
pdf(file = paste0(plot.path, sample, "_eps_data_", eps, "_cell_states.pdf"), width = 28, height = 36)
p <- nst %>% ggplot(aes(x = stg.diff)) + geom_histogram(binwidth = 3, fill = "#69b3a2", color="#e9ecef", alpha=0.9) + ggtitle("Mean_stage_difference")
y <- DimPlot(obj, group.by = "cell.states.24") + ggplot2::ggtitle(paste0("states above 24 hr"))
z <- DimPlot(obj, group.by = "cell.states.36") + ggplot2::ggtitle(paste0("states above 36 hr"))
h <- DimPlot(obj, group.by = "cell.states.48") + ggplot2::ggtitle(paste0("states above 48 hr")) 
gridExtra::grid.arrange(grobs = list(p, y, z, h), ncol = 2)
dev.off()
}


