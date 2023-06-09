#Subsetting the new cropped subsets from new mama object 
##This object has all curated annotations at the level of broader tissues and more granular cell types

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

# Set relative paths
old.mama.obj.path <- paste0("/data/CSD/zfext/results/02-Clustering/merged_timepoints/obj/merged_ALLhpf__igraph-weighted.rds")
new.mama.obj.path <- paste0("/data/CSD/zfext/LTA/results/02-Clustering/merged_timepoints/obj/merged_mama_ds_integrated_2021-01-11_seurat_UMAP_recalc.rds")

#base.path <- "~/Documents/zfext-biowulf/02-Clustering/"

##Set the other paths to save counts, variable genes, plots, meta.data and the tissue specific seurat objects
counts.path <- paste0(base.path, "/count_subsets_v6/")
obj.path <- paste0(base.path, "/mama_crops_v6/")
var.path <- paste0(base.path, "/var_subsets_v6/")
plot.path <- paste0(base.path, "/plots_subsets_v6/")
meta.path <- paste0(base.path, "/meta_subsets_v6/")
seurat.path <- paste0(base.path, "/obj_subsets_v6/")
cells.path <- paste0(base.path, "/cells_subsets_v6/")
save.path <- paste0(base.path, "/cells_from_object/")

##First we need to create a seurat object for each tissue (here labeled as sample)
message(paste0(Sys.time(), ": Creating subset for", sample))

#Load the global dataset (merged 2018 Dropseq and 10X datasets)
message(paste0(Sys.time(), ": Loading mama"))
mama <- readRDS("/data/CSD/zfext/LTA/results/02-Clustering/merged_timepoints/obj/merged_mama_ds_integrated_2021-01-11_seurat_UMAP_recalc.rds")

#Loading tissue specific clustering object where the cells are already sorted into various tissue categories
message(paste0(Sys.time(), ": Loading tissue clustering object"))
tissue.annot <- readRDS("/data/CSD/zfext/LTA/results/02-Clustering/merged_timepoints/meta.data/tissue_annot_2.rds")

##Add the tissue cluster information to the metadata of the global data seurat object
message(paste0(Sys.time(), ": Adding the tissue clustering to mama"))
mama@meta.data$tissue.annot <- tissue.annot

#message(paste0(Sys.time(), ": Retrieve cell lists for ", sample))
cells.list <- rownames(mama@meta.data)[which(mama@meta.data$tissue.annot == sample)]

message(paste0(Sys.time(), ": Writing cells for", sample))
write(cells.list, paste0(cells.path, sample, "_cell_ids.txt"))
cells.length <- length(cells.list)

message(paste0(Sys.time(), ": ", cells.length, " cells are present in" , sample))

message(paste0(Sys.time(), ": Align cells with that in mama"))
cell.length.init <- length(cells.list)
cells.keep <- intersect(cells.list, colnames(mama@assays$RNA@counts))
cell.length.intersect <- length(cells.keep)
message(paste0(Sys.time(), ": ", cell.length.intersect, " of ", cell.length.init, " cells in the list were present in this object."))


# Get relevant counts and metadata
message(paste0(Sys.time(), ": Getting relevant counts and metadata"))
subset.counts <- mama@assays$RNA@counts[, cells.keep]
subset.meta <- mama@meta.data[cells.keep, ]

message(paste0(Sys.time(), ": Create Seurat object"))
obj <- CreateSeuratObject(counts = subset.counts, meta.data = subset.meta, min.features = 0, min.cells = 2)

rm(list = c("mama"))
shh <- gc()

# Write out sample count matrices in MTX format needed for gene module analysis
message(paste0(Sys.time(), ": Savings counts as MTX"))
save.path <- paste0(counts.path, sample, ".mtx")
Matrix::writeMM(obj = obj@assays$RNA@counts, file = save.path)
write(rownames(obj@assays$RNA@counts), file = paste0(save.path, ".rn"))
write(colnames(obj@assays$RNA@counts), file = paste0(save.path, ".cn"))

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

# Normalize and regress
message(paste0(Sys.time(), ": Normalizing Data"))
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

##Regress out mito and ribo
message(paste0(Sys.time(), ": Regressing out mito/ribo"))
obj <- ScaleData(obj, vars.to.regress = c("percent.mt", "percent.ribo"))

# Variable Feature Selection
message(paste0(Sys.time(), ": Variable Feature Selection"))	
obj <- FindVariableFeatures(obj, nfeatures = 2000)

# Save object and variable genes for downstream DR + clustering?
message(paste0(Sys.time(), ": Saving results"))
write(obj@assays$RNA@var.features, file=paste0(var.path, sample, "_vargenes.txt"))
saveRDS(obj, file=paste0(seurat.path, sample, "_seurat.rds"))

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

message(paste0(Sys.time(), ": Dimensions used equal to :", dims.use))

# Plot variable feature selection and PCA selection
pdf(file=paste0(var.path, sample, "_var.pdf"), width=8, height=8)
plot(VariableFeaturePlot(obj) + ggplot2::ggtitle(sample))
ElbowPlot(obj, ndims=100) + ggplot2::ggtitle(sample) + ggplot2::geom_vline(xintercept=dims.use+0.5, color='red')
dev.off()


# Save object - temporarily in case of failure
message(paste0(Sys.time(), ": Saving - temp"))
saveRDS(obj, file=paste0(seurat.path, sample, "_seurat.rds"))

# Dimensionality Reduction: UMAP
message(paste0(Sys.time(), ": UMAP"))
#nn.use <- floor(sqrt(ncol(obj@assays$RNA@counts))/10)*10
nn.use <- 50
obj <- RunUMAP(obj, dims=1:dims.use, n.neighbors=nn.use)

# Save object - temporarily in case of failure
message(paste0(Sys.time(), ": Saving - temp"))
saveRDS(obj, file=paste0(seurat.path, sample, "_seurat.rds"))

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
saveRDS(obj, file=paste0(seurat.path, sample, "_seurat.rds"))


# Save the metadata table
write.table(obj@meta.data, file=paste0(meta.path, sample, "_meta.txt"))

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

##Generating the UMAP plot where cells are categorized by their stage (hpf) - Figure 1B
##Define the colors for the stages
stage.colors.new <- c(
  colorRampPalette(c("#f4be1d", "#efb918", "#bd8700", "#b37d00", "#a97300", "#9f6900", "#955f00", "#8b5500"))(8), # 3, 4, 5, 6, 7, 8, 9, 10
  colorRampPalette(c("#FFC0CB", "#FFB6C1", "#FF69B4", "#DA70D6", "#BA55D3", "#C71585", "#FF1493"))(7), # 11, 12, 14, 16, 18, 21, 24
  rep(RColorBrewer::brewer.pal(6, "Oranges"), each = 2), # 26 / 28, 30 / 32, 34 / 36, 38 / 40, 42 / 44, 46
  rep(RColorBrewer::brewer.pal(6, "Greens"), each = 2), # 48, 50 / 52, 54 / 56, 58 / 60, 62 / 64, 66 / 68, 70,  72
  rep(RColorBrewer::brewer.pal(6, "Blues"), each = 2),  # 74 / 76, 78 / 80, 82 / 84, 86 / 88, 90 / 92, 94, 96
  rep(RColorBrewer::brewer.pal(7, "Purples"), each = 2) #  98 / 100, 102 / 104, 106 / 108, 110 / 112, 114 / 116, 118 / 120
)

dpi <- 300
png("~/Desktop/mama_stage.png", width = 5*dpi, height = 5*dpi)
DimPlot(mama, group.by = "stage.nice", cols = stage.colors.new, raster = F)
dev.off()


##Plot the broader tissue catgories on the UMAP plot - Figure 1C and Fig S1B
# Load cell annotations to introduce for downsample purposes
cell.annot.path <- "~/Box/zfext/annotations_celltype_curated_newMama/celltype_annotations_df.tsv"
message(Sys.time(), ": Loading ", cell.annot.path)
cell.annot <- read.table(cell.annot.path, sep = " ", header = T, stringsAsFactors = F)
cell.annot$stage.group <- plyr::mapvalues(
  from = c("3-4", "5-6", "7-9", " 10-12", " 14-21", " 24-34", " 36-46", 
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


## CREATE ADDITIONAL META DATA IN MAMA --------------------

# Create any additional needed meta.data
mama@meta.data$hpf <- as.numeric(mama@meta.data$stage.nice)

##Replace all underscores in the tissue.subset slot
mama@meta.data$tissue.annot <- plyr::mapvalues(
  from = c("blastula", "periderm", "axial_mesoderm", "gastrula", "neural", "hematopoietic", "muscle_superset", "endoderm", "cephalic_mesoderm",
           "PGCs", "pronephros", "unknown", "mesenchyme", "basal_epidermis", "ionocytes_mucous-secreting", "taste_olfactory", "glial_cells", "eye",
           "otic", "pigment-cells", "fin", "non-skeletal_muscle"),
  to =   c("blastula", "periderm", "axial", "gastrula", "neural", "hematopoietic", "muscle", "endoderm", "cephalic",
           "PGCs", "pronephros", "unknown", "mesenchyme", "epidermis", "ionocytes", "taste", "glial", "eye",
           "otic", "pigment", "fin", "mural"),
  x = mama@meta.data$tissue.subsets
)

# Create stage.group / cluster
mama@meta.data$tissue.group <- URD::cutstring(delimiter = "-", field = 1, x = 
                                                gsub(" ", "", 
                                                     apply(mama@meta.data[,c("tissue.annot", "stage.group")], 1, paste0, collapse = "_")
                                                )
)

mama@meta.data$tissue.hpf <- gsub(" ", "", apply(mama@meta.data[,c("tissue.annot", "hpf")], 1, paste0, collapse = "_"))
mama@meta.data$tissue_stage.diff <- paste0(mama@meta.data$clust, "_", gsub(" ", "", mama@meta.data$stage.diff))

##Now load the cell-embeddings calculated on mama
all.emb <- mama.full.red@cell.embeddings
all.emb <- as.data.frame(all.emb)
all.emb$tissue.clust <- mama@meta.data$tissue.group

# Load cluster annotations 
is.empty <- function(x) {
  is.na(x) | is.null(x) | (sapply(trimws(x), nchar) == 0)
}

annot.confirmed <- readxl::read_xlsx("~/Box/zfext/02-Clustering/mama_tissue_clustering/annot/mama_annotations_full_extended_jeffedit_newtissues_curated_clustering.xlsx")
annot.confirmed <- annot.confirmed[!is.empty(annot.confirmed$clust),3:ncol(annot.confirmed)]
rownames(annot.confirmed) <- annot.confirmed$clust
annot.confirmed$confirmed <- T

##Add annotations to mama metadata
mama@meta.data$identity.super <- NA
for(i in annot.confirmed$clust){
  mama@meta.data[which(mama@meta.data$clust == i), "identity.super"] <- annot.confirmed[which(annot.confirmed$clust == i), "tissue"]
}

mama@meta.data[rownames(mama@meta.data)[which(mama@meta.data$tissue.subsets == "cephalic")], "identity.super"] <- "cephalic"

##Plot UMAP with tissue categories - shown in Figure 1C
dpi <- 300
png("mama_curated_clustering_new_clusters.png", width = 5*dpi, height = 5*dpi)
DimPlot(mama, group.by = "identity.super", raster = F, cols = cluster_colors, na.value = "#F5F5F5") + NoAxes() + NoLegend()
dev.off()


#Split the stage plots into larger groups and plot them separately - Figure S1A
##As the stage plot was too hard to visualize, here I split the stage UMAP plot into smaller pieces according to the colors below
stage.colors.new <- c(
  colorRampPalette(c("#f4be1d", "#efb918", "#bd8700", "#b37d00", "#a97300", "#9f6900", "#955f00", "#8b5500"))(8), # 3, 4, 5, 6, 7, 8, 9, 10
  colorRampPalette(c("#FFC0CB", "#FFB6C1", "#FF69B4", "#DA70D6", "#BA55D3", "#C71585", "#FF1493"))(7), # 11, 12, 14, 16, 18, 21, 24
  rep(RColorBrewer::brewer.pal(6, "Oranges"), each = 2), # 26 / 28, 30 / 32, 34 / 36, 38 / 40, 42 / 44, 46
  rep(RColorBrewer::brewer.pal(6, "Greens"), each = 2), # 48, 50 / 52, 54 / 56, 58 / 60, 62 / 64, 66 / 68, 70
  rep(RColorBrewer::brewer.pal(6, "Blues"), each = 2),  # 72, 74 / 76, 78 / 80, 82 / 84, 86 / 88, 90 / 92, 94
  rep(RColorBrewer::brewer.pal(7, "Purples"), each = 2) # 96, 98 / 100, 102 / 104, 106 / 108, 110 / 112, 114 / 116, 118 / 120
)

colors.stg.3.10 <- setNames(colorRampPalette(c("#f4be1d", "#efb918", "#bd8700", "#b37d00", "#a97300", "#9f6900", "#955f00", "#8b5500"))(8), unique(mama@meta.data$stage.nice)[1:8])
colors.stg.12.21 <- setNames(colorRampPalette(c("#FFC0CB", "#FFB6C1", "#FF69B4", "#DA70D6", "#BA55D3", "#C71585", "#FF1493"))(7), unique(mama@meta.data$stage.nice)[9:15])
colors.stg.24.46 <- setNames(rep(RColorBrewer::brewer.pal(6, "Oranges"), each = 2), unique(mama@meta.data$stage.nice)[16:27])
colors.stg.48.70 <- setNames(rep(RColorBrewer::brewer.pal(6, "Greens"), each = 2), unique(mama@meta.data$stage.nice)[28:39])
colors.stg.72.94 <- setNames(rep(RColorBrewer::brewer.pal(6, "Blues"), each = 2), unique(mama@meta.data$stage.nice)[40:51])
colors.stg.96.120 <- setNames(rep(RColorBrewer::brewer.pal(7, "Purples"), each = 2), unique(mama@meta.data$stage.nice)[52:63])

dpi <- 300
png("~/Box/zfext/results/04c-SubsetsV3/figures_for_manuscript/Figure_1_mama_cell_atlas/mama_stage_3-10hpf.png", width = 5*dpi, height = 5*dpi)
DimPlot(mama, group.by = "stage.nice", cols = colors.stg.3.10, raster = F, na.value = "#f5f7f6") + NoLegend() + NoAxes()
dev.off()

png("~/Box/zfext/results/04c-SubsetsV3/figures_for_manuscript/Figure_1_mama_cell_atlas/mama_stage_12-21hpf.png", width = 5*dpi, height = 5*dpi)
DimPlot(mama, group.by = "stage.nice", cols = colors.stg.12.21, raster = F, na.value = "#f5f7f6") + NoLegend() + NoAxes()
dev.off()

png("~/Box/zfext/results/04c-SubsetsV3/figures_for_manuscript/Figure_1_mama_cell_atlas/mama_stage_24-46hpf.png", width = 5*dpi, height = 5*dpi)
DimPlot(mama, group.by = "stage.nice", cols = colors.stg.24.46, raster = F, na.value = "#f5f7f6") + NoLegend() + NoAxes()
dev.off()

png("~/Box/zfext/results/04c-SubsetsV3/figures_for_manuscript/Figure_1_mama_cell_atlas/mama_stage_48-70hpf.png", width = 5*dpi, height = 5*dpi)
DimPlot(mama, group.by = "stage.nice", cols = colors.stg.48.70, raster = F, na.value = "#f5f7f6") + NoLegend() + NoAxes()
dev.off()

png("~/Box/zfext/results/04c-SubsetsV3/figures_for_manuscript/Figure_1_mama_cell_atlas/mama_stage_72-94hpf.png", width = 5*dpi, height = 5*dpi)
DimPlot(mama, group.by = "stage.nice", cols = colors.stg.72.94, raster = F, na.value = "#f5f7f6") + NoLegend() + NoAxes()
dev.off()

png("~/Box/zfext/results/04c-SubsetsV3/figures_for_manuscript/Figure_1_mama_cell_atlas/mama_stage_96-120hpf.png", width = 5*dpi, height = 5*dpi)
DimPlot(mama, group.by = "stage.nice", cols = colors.stg.96.120, raster = F, na.value = "#f5f7f6") + NoLegend() + NoAxes()
dev.off()






