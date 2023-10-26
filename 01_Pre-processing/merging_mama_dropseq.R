##This script is to combine the 2018 DropSeq dataset (3.3-12hpf) remapped to the Lawson transcriptome with the newly generated dataset (14-120 hpf) using the 10X genomics platform to create the final composite dataset (3.3-120 hpf). This was done using the two Seurat objects separately created for the 2018 dataset and the newly generated 10X genomics dataset 

##Load libraries
library(Seurat)
library(URD)
#Loading required package: ggplot2
#Loading required package: Matrix
library(ggplot2)
library(Matrix)

##Load the newly generated 10X dataset seurat object mapped to the Lawson transcriptome
mama <- readRDS("/data/CSD/zfext/LTA/results/02-Clustering/merged_timepoints/obj/merged_ALLhpf_210122_keep.rds")

##load  the 2018 DropSeq dataset seurat object mapped to the Lawson transcriptome
obj.seurat <- readRDS("/data/CSD/2018-Dropseq/obj/zf2018_seurat_LTA.rds")

##Load the 2018 DropSeq dataset URD object
obj.urd <- readRDS("/data/CSD/2018-Dropseq/obj/urd_zf_LTA.rds")

##To merge the two seurat objects, the merge command was used. The merge command merges the raw count matrices of the two seurat objects and creates a new Seurat object with a resultant combined raw counts matrix. No additional cell ids were added as the cell ids for these datasets were already different.                   
obj <- merge(obj.seurat, mama)

# Mitochondrial and ribosomal genes: 
mito.genes <- grep("^mt-", rownames(obj@assays$RNA@counts), value=T)
ribo.genes <- intersect(rownames(obj@assays$RNA@counts), c("rpl18a", "rps16", "rplp2l", "rps13", "rps17", "rpl34", 
                                                           "rpl13", "rplp0", "rpl36a", "rpl12", "rpl7a", 
                                                           "rpl19", "rps2", "rps15a", "rpl3", "rpl27", "rpl23", "rps11", 
                                                           "rps27a", "rpl5b", "rplp2", "rps26l", "rps10", 
                                                           "rpl5a", "rps28", "rps8a", "rpl7", "rpl37", "rpl24", "rpl9", 
                                                           "rps3a", "rps6", "rpl8", "rpl31", "rpl18", "rps27.2", "rps19", 
                                                           "rps9", "rpl28", "rps7", "rpl7l1", "rps29", "rpl6", 
                                                           "rps8b", "rpl10a", "rpl13a", 
                                                           "rpl39", "rpl26", "rps24", "rps3", "rpl4", 
                                                           "rpl35a", "rpl38", "rplp1", "rps27.1", "rpl15", "rps18", "rpl30", 
                                                           "rpl11", "rpl14", "rps5", "rps21", "rpl10", "rps26", 
                                                           "rps12", "rpl37.1", "rpl35", "rpl17", "rpl23a", "rps14", "rpl29", 
                                                           "rps15", "rpl22", "rps23", "rps25", "rpl21", 
                                                           "rpl22l1", "rpl36", "rpl32", "rps27l"))
obj[["percent.mt"]] <- PercentageFeatureSet(obj, features = mito.genes)
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, features = ribo.genes)

message(paste0(Sys.time(), ": Saving results"))
saveRDS(obj, file="/data/CSD/zfext/results/merged_mama_DS_integrated_12012021_seurat.rds")

# Normalize and regress
message(paste0(Sys.time(), ": Normalizing Data"))
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

message(paste0(Sys.time(), ": Saving results"))
saveRDS(obj, file="/data/CSD/zfext/results/merged_mama_DS_integrated_12012021_seurat.rds")

message(paste0(Sys.time(), ": Regressing out mito/ribo"))
obj <- ScaleData(obj, vars.to.regress = c("percent.mt", "percent.ribo"))

message(paste0(Sys.time(), ": Saving results"))
saveRDS(obj, file="/data/CSD/zfext/results/merged_mama_DS_integrated_12012021_seurat.rds")

# Variable Feature Selection
message(paste0(Sys.time(), ": Variable Feature Selection"))
obj <- FindVariableFeatures(obj, nfeatures = 2000)
# Save object and variable genes for downstream DR + clustering?
message(paste0(Sys.time(), ": Saving results"))
write(obj@assays$RNA@var.features, file=paste0("/data/CSD/zfext/results/merged_mama_DS_integrated_12012021_vargenes.txt"))
saveRDS(obj, file="/data/CSD/zfext/results/merged_mama_DS_integrated_12012021_seurat.rds")

# Dimensionality Reduction: PCA
message(paste0(Sys.time(), ": PCA"))
obj <- RunPCA(obj, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
obj <- JackStraw(obj, assay="RNA", dims=100)
obj <- ScoreJackStraw(obj, dims = 1:100)

if (any(obj@reductions$pca@jackstraw@overall.p.values[,2] > 1e-3)) {
  dims.use <- min(which(obj@reductions$pca@jackstraw@overall.p.values[,2] > 1e-3)) - 1
} else {
  dims.use <- 30
}

# Plot variable feature selection and PCA selection
pdf(file= paste0("/data/CSD/zfext/results/merged_mama_DS_integrated_12012021_var.pdf"), width=8, height=8)
plot(VariableFeaturePlot(obj) + ggplot2::ggtitle(sample))
ElbowPlot(obj, ndims=100) + ggplot2::ggtitle(sample) + ggplot2::geom_vline(xintercept=dims.use+0.5, color='red')
dev.off()


# Save object - temporarily in case of failure
message(paste0(Sys.time(), ": Saving - temp"))
saveRDS(obj, file="/data/CSD/zfext/results/merged_mama_DS_integrated_12012021_seurat.rds")

# Dimensionality Reduction: UMAP
message(paste0(Sys.time(), ": UMAP"))
#nn.use <- floor(sqrt(ncol(obj@assays$RNA@counts))/10)*10
nn.use <- 50
obj <- RunUMAP(obj, dims=1:dims.use, n.neighbors=nn.use)

# Save object - temporarily in case of failure
message(paste0(Sys.time(), ": Saving - temp"))
saveRDS(obj, file="/data/CSD/zfext/results/merged_mama_DS_integrated_12012021_seurat.rds")

# Clustering
message(paste0(Sys.time(), ": Clustering"))
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:dims.use, nn.eps = 0, k.param = nn.use)

if(ncol(obj) < 35000) {
  cat('--> Leiden algorithm, using defaults (method="matrix", weights=NULL) ...\n')
  obj <- FindClusters(obj, algorithm=4, resolution = c(0.75,1,2,3), n.start=50, random.seed=17)
} else {
  cat('--> Leiden algorithm; using method="igraph" with weights=T ...\n')
  obj <- FindClusters(obj, algorithm=4, resolution = c(0.75,1,2,3), n.start=50, random.seed=17, method='igraph', weights=T)
}
Idents(obj) <- 'RNA_snn_res.2'

# Save object
message(paste0(Sys.time(), ": Saving"))
saveRDS(obj, file="/data/CSD/zfext/results/merged_mama_DS_integrated_12012021_seurat.rds")

# Save the metadata table
write.table(obj@meta.data, file=paste0("/data/CSD/zfext/results/merged_mama_DS_integrated_12012021_meta.txt"))

message(paste0(Sys.time(), ": Fixing all stage information"))
obj@meta.data$stage.nice <- NA
for (i in obj@meta.data$timepoint) {
  obj@meta.data[which(obj@meta.data$timepoint == i), "stage.nice"] <- as.numeric(sapply(strsplit(i, split='hpf', fixed=TRUE), function(x) (x[1])))
}

for (i in obj@meta.data$hpf) {
  obj@meta.data[which(obj@meta.data$hpf == i), "stage.nice"] <- i
}

##Group stages into smaller bins to create a slot called stage.groups
message(paste0(Sys.time(), ": Creating a slot for stage groups"))
obj@meta.data[obj@meta.data$stage.nice == 5.3, "stage.nice"] <- 5
obj@meta.data[obj@meta.data$stage.nice == 4.7, "stage.nice"] <- 4
obj@meta.data[obj@meta.data$stage.nice == 4.3, "stage.nice"] <- 4
obj@meta.data[obj@meta.data$stage.nice == 3.3, "stage.nice"] <- 3

obj@meta.data$stage.group <- plyr::mapvalues(x = obj@meta.data$stage.nice,
                                             from = c(seq(3, 12, 1), 14, 16, 18, 21, seq(24, 120, 2)),
                                             to = c(rep(c("3-4", "5-6"), each = 2), rep(c("7-9", "10-12"), each = 3),
                                                    rep("14-21", 4), 
                                                    rep(c("24-34", "36-46", "48-58", "60-70", "72-82", "84-94", "96-106", "108-118"), each = 6),
                                                    "120"
                                             )
)

# Save the metadata table again
write.table(obj@meta.data, file=paste0("/data/CSD/zfext/results/merged_mama_DS_integrated_12012021_meta.txt"))

# Save object
message(paste0(Sys.time(), ": Saving"))
saveRDS(obj, file="/data/CSD/zfext/results/merged_mama_DS_integrated_12012021_seurat.rds")

# Do plots
message(paste0(Sys.time(), ": Plotting basic"))
pdf(file=paste0(plot.path, sample, "_clusters.pdf"), width = 8, height = 8)
DimPlot(obj, group.by="stage") + ggplot2::ggtitle(paste0(sample, ": stage"))
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





