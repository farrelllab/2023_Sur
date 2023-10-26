## This code is for plotting the UMAP and markers plots for the newly cropped tissue specific seurat objects

# Get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("Must provide at least one argument.")
sample = as.character(args[1])

library(Seurat)
library(Matrix)
library(gridExtra)

# Set base path
base.path <- "/data/CSD/zfext/results/04f-SubsetsV6"
old.base.path <- "/data/CSD/zfext/results/04d-SubsetsV4"

# Set relative paths
#old.mama.obj.path <- paste0("/data/CSD/zfext/results/02-Clustering/merged_timepoints/obj/merged_ALLhpf__igraph-weighted.rds")
#new.mama.obj.path <- paste0("/data/CSD/zfext/LTA/results/02-Clustering/merged_timepoints/obj/merged_ALLhpf_201028.rds")

#base.path <- "~/Documents/zfext-biowulf/02-Clustering/"

counts.path <- paste0(base.path, "/count_subsets_v6/")
obj.path <- paste0(base.path, "/obj_subsets_v6/")
var.path <- paste0(base.path, "/var_subsets_v6/")
plot.path <- paste0(base.path, "/plots_subsets_v6/")
meta.path <- paste0(base.path, "/meta_subsets_v6/")
seurat.path <- paste0(base.path, "/obj_subsets_v6/")


# Load data
message(paste0(Sys.time(), ": Loading data for: ", sample))
obj <- readRDS(file=paste0(obj.path, sample, "_seurat.rds"))

# Find markers
message(paste0(Sys.time(), ": Finding markers"))
Idents(obj) <- 'RNA_snn_res.2'
posMarkers.wilcox <- sapply(levels(Idents(obj)), function(i) FindMarkers(obj,i,assay.type='RNA',test.use="wilcox",only.pos=T, min.cells.group = 1),simplify=F)
posMarkers.roc <- sapply(levels(Idents(obj)), function(i) FindMarkers(obj,i,test.use="roc",only.pos=T),simplify=F)
markers <- list(wilcox = posMarkers.wilcox, roc = posMarkers.roc)
saveRDS(markers, file=paste0(seurat.path, sample, "_markers.rds"))

# Save object
message(paste0(Sys.time(), ": Saving"))
saveRDS(obj, file=paste0(seurat.path, sample, "_seurat.rds"))

# Do plots
message(paste0(Sys.time(), ": Plotting basic"))
pdf(file=paste0(plot.path, sample, "_clusters.pdf"), width = 8, height = 8)
DimPlot(obj, group.by="stage.nice") + ggplot2::ggtitle(paste0(sample, ": stage.nice"))
DimPlot(obj, group.by = "stage.group", split.by="stage.group", ncol=3) + ggplot2::ggtitle(paste0(sample)) + ggplot2::theme(legend.position = "none")
DimPlot(obj, group.by="RNA_snn_res.2", label = T) + ggplot2::theme(legend.position = "none") + ggplot2::ggtitle(paste0(sample, ": res 2"))
DimPlot(obj, group.by="RNA_snn_res.3", label = T) + ggplot2::theme(legend.position = "none") + ggplot2::ggtitle(paste0(sample, ": res 3"))
dev.off()

# Do cluster plots
message(paste0(Sys.time(), ": Plotting markers"))
pdf(file=paste0(plot.path, sample, "_markers.pdf"), width = 16, height = 16)
for (i in levels(Idents(obj))) {
  genes.plot <- unique(c(rbind(rownames(posMarkers.wilcox[[i]])[1:8], rownames(posMarkers.roc[[i]])[1:8])))[1:8]
  b <- lapply(genes.plot, FeaturePlot, object = obj)
  b[[9]] <- DimPlot(obj, group.by = "RNA_snn_res.2", cells.highlight = WhichCells(obj, idents = i)) + ggplot2::ggtitle(paste0(sample, ": ", i)) + ggplot2::theme(legend.position = "none")
  gridExtra::grid.arrange(grobs=b)
}
dev.off()

