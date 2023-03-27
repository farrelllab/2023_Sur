##This code is to understand what cell states persist for longer periods and what cell states persist for shorter periods. 
library(Seurat)
library(URD)
library(dbscan)
library(RColorBrewer)

##Load mama
mama <- readRDS("~/Box/zfext/02-Clustering/mama+ds_integrated/merged_mama_ds_slim_seurat_25Jan2022.rds")

#Load the reductions slot for mama
mama.reductions <- readRDS("~/Desktop/merged_mama_ds_reductions_2022_01-18.rds")
mama@reductions <- mama.reductions

# Load cell annotations to introduce for downsample purposes
cell.annot.path <- "~/Box/zfext/annotations_celltype_curated_newMama/celltype_annotations_df_full.tsv"
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


## CREATE ADDITIONAL META DATA IN MAMA --------------------

# Create any additional needed meta.data
mama@meta.data$hpf <- as.numeric(mama@meta.data$stage.nice)

##Replace all underscores in the tissue.subset slot
mama@meta.data$tissue.annot <- plyr::mapvalues(
  from = c("blastula", "periderm", "axial_mesoderm", "gastrula", "neural", "hematopoietic", "muscle_superset", "endoderm", "cephalic_mesoderm",
           "PGCs", "pronephros", "cephalic", "mesenchyme", "basal_epidermis", "ionocytes_mucous-secreting", "taste_olfactory", "glial_cells", "eye",
           "otic", "pigment-cells", "fin", "non-skeletal_muscle"),
  to =   c("blastula", "periderm", "axial", "gastrula", "neural", "hematopoietic", "muscle", "endoderm", "cephalic",
           "PGCs", "pronephros", "cephalic", "mesenchyme", "epidermis", "ionocytes", "taste", "glial", "eye",
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


##Load the various stage.difference objects calculated on the cluster and plot on a UMAP how long individual cells persist for
##Ok, now generate a script such that we can sort every single cell based on their stage difference
##Load the "nearest" objects for each tissue and save the values as a slot in mama
sample <- "axial"
eps <- "40"
nearest.axial <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.axial <- as.data.frame(nearest.axial)
colnames(nearest.axial) <- "stage.diff"

sample <- "epidermis"
nearest.basal <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.basal <- as.data.frame(nearest.basal)
colnames(nearest.basal) <- "stage.diff"

sample <- "blastomeres"
nearest.blastula <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.blastula <- as.data.frame(nearest.blastula)
colnames(nearest.blastula) <- "stage.diff"

sample <- "endoderm"
nearest.endo <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.endo <- as.data.frame(nearest.endo)
colnames(nearest.endo) <- "stage.diff"

sample <- "eye"
nearest.eye <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.eye <- as.data.frame(nearest.eye)
colnames(nearest.eye) <- "stage.diff"

sample <- "fin"
nearest.fin <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.fin <- as.data.frame(nearest.fin)
colnames(nearest.fin) <- "stage.diff"

sample <- "glial"
nearest.glia <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.glia <- as.data.frame(nearest.glia)
colnames(nearest.glia) <- "stage.diff"

sample <- "hematopoietic"
nearest.hema <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.hema <- as.data.frame(nearest.hema)
colnames(nearest.hema) <- "stage.diff"

sample <- "ionocytes"
nearest.iono <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.iono <- as.data.frame(nearest.iono)
colnames(nearest.iono) <- "stage.diff"

sample <- "mesenchyme"
nearest.mese <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.mese <- as.data.frame(nearest.mese)
colnames(nearest.mese) <- "stage.diff"

sample <- "mural"
nearest.mural <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.mural <- as.data.frame(nearest.mural)
colnames(nearest.mural) <- "stage.diff"

sample <- "muscle"
nearest.muscle <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.muscle <- as.data.frame(nearest.muscle)
colnames(nearest.muscle) <- "stage.diff"

sample <- "neural"
nearest.neural <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.neural <- as.data.frame(nearest.neural)
colnames(nearest.neural) <- "stage.diff"

sample <- "otic"
nearest.otic <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.otic <- as.data.frame(nearest.otic)
colnames(nearest.otic) <- "stage.diff"

sample <- "periderm"
nearest.periderm <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.periderm <- as.data.frame(nearest.periderm)
colnames(nearest.periderm) <- "stage.diff"

sample <- "pgc"
nearest.pgc <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.pgc <- as.data.frame(nearest.pgc)
colnames(nearest.pgc) <- "stage.diff"

sample <- "pigment"
nearest.pigment <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.pigment <- as.data.frame(nearest.pigment)
colnames(nearest.pigment) <- "stage.diff"

sample <- "pronephros"
nearest.pron <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.pron <- as.data.frame(nearest.pron)
colnames(nearest.pron) <- "stage.diff"

sample <- "taste"
nearest.taste <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.taste <- as.data.frame(nearest.taste)
colnames(nearest.taste) <- "stage.diff"

nearest.total <- rbind(nearest.axial, nearest.basal, nearest.blastula, nearest.endo, nearest.eye, nearest.fin, nearest.gast, nearest.glia, nearest.hema,
                       nearest.iono, nearest.mese, nearest.mural, nearest.muscle, nearest.neural, nearest.otic, nearest.periderm, nearest.pgc, 
                       nearest.pigment, nearest.pron, nearest.taste)

nearest.total <- rbind(nearest.axial, nearest.basal, nearest.blastula, nearest.endo, nearest.eye, nearest.fin, nearest.glia, nearest.hema,
                       nearest.iono, nearest.mese, nearest.mural, nearest.muscle, nearest.neural, nearest.otic, nearest.periderm, nearest.pgc, 
                       nearest.pigment, nearest.pron, nearest.taste)

cells.all <- WhichCells(mama)
cells.rest <- setdiff(rownames(mama@meta.data), rownames(nearest.total))
cells.rest <- rownames(mama@meta.data)[mama@meta.data$tissue.subsets == "cephalic"]

obj <- subset(mama, cells = cells.rest)
message(paste0(Sys.time(), ": Normalizing Data"))
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
message(paste0(Sys.time(), ": Regressing out mito/ribo"))
obj <- ScaleData(obj, vars.to.regress = c("percent.mt", "percent.ribo"))
# Variable Feature Selection
message(paste0(Sys.time(), ": Variable Feature Selection"))
obj <- FindVariableFeatures(obj, nfeatures = 2000)
# Save object and variable genes for downstream DR + clustering?
message(paste0(Sys.time(), ": Saving results"))
write(obj@assays$RNA@var.features, file=paste0(var.path, sample, "_vargenes.txt"))
seurat.path <- "~/Box/zfext/annotations_celltype_curated_newMama/cephalic/"
saveRDS(obj, file=paste0(seurat.path, sample, "_seurat.rds"))


# Dimensionality Reduction: PCA
message(paste0(Sys.time(), ": PCA"))
obj <- RunPCA(obj, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
obj <- JackStraw(obj, assay="RNA", dims=100)
obj <- ScoreJackStraw(obj, dims = 1:100)

dims.use <- 30

# Plot variable feature selection and PCA selection
pdf(file=paste0(var.path, sample, "_var.pdf"), width=8, height=8)
plot(VariableFeaturePlot(obj) + ggplot2::ggtitle(sample))
ElbowPlot(obj, ndims=100) + ggplot2::ggtitle(sample) + ggplot2::geom_vline(xintercept=dims.use+0.5, color='red')
dev.off()

nn.use <- 50
obj <- RunUMAP(obj, dims=1:18, n.neighbors=nn.use)

obj <- FindNeighbors(obj, reduction = "pca", dims = 1:dims.use, nn.eps = 0, k.param = nn.use)
obj <- FindClusters(obj, algorithm=4, resolution = c(2,3), n.start=50, random.seed=17)

dist.mat  <- dist(t(obj@assays$RNA@data))
mat <- as.matrix(dist.mat)

save.path <- "~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/"
plot.path <- "~/Box/zfext/global_analysis_curated/transient_vs_long-term/plots_mat_data_eps_compare/"
sample <- "cephalic"
for(eps in c(40, 40, 42, 43, 44, 45)){
  nearest <- sapply(rownames(mat), function(cell) mean(abs(obj@meta.data[setdiff(colnames(mat)[which(mat[cell,] <= eps)], cell), "stage.nice"] - obj@meta.data[cell, "stage.nice"])))
  #nn <- sapply(rownames(mat), function(cell) setdiff(colnames(mat)[which(mat[cell, ] <= eps)], cell))
  #neighbors <- list(nearest, nn)
  #names(neighbors) <- c("mean.stage.diff", "neighbors")
  
  message(paste0(Sys.time(), ": Saving"))
  saveRDS(nearest, file = paste0(save.path, sample, "_mean_stage_diff_eps_", eps, ".rds"))
  #saveRDS(nn, file = paste0(save.path, sample, "_neighbors_eps_", eps, ".rds"))
  
  sample <- "cephalic"
  eps <- 40
  nearest.cephalic <- readRDS(file = paste0(save.path, sample, "_mean_stage_diff_eps_", eps, ".rds"))
  
  ##Long-term states would be all those states where the stage differences are atleast more than 2 days
  message(paste0(Sys.time(), ": Get long-term and short-term states"))
  long.term.cephalic <- names(nearest.cephalic)[which(nearest.cephalic >= 48)]
  ##Short term states would be states that last less than 8 hours
  short.term.cephalic <- names(nearest.cephalic)[which(nearest.cephalic < 48)]
  
  #Calculate proportion of long-term states per cluster
  message(paste0(Sys.time(), ": get proportion of long and short-term states per clusters"))
  long.term.clust.prop <- table(obj@meta.data[long.term.cephalic, "RNA_snn_res.2"])/table(obj@meta.data$RNA_snn_res.2)
  short.term.clust.prop <- table(obj@meta.data[short.term.cephalic, "RNA_snn_res.2"])/table(obj@meta.data$RNA_snn_res.2)
  
  message(paste0(Sys.time(), ": Convert proportions into a dataframe"))
  long.prop <- as.data.frame(long.term.clust.prop)
  colnames(long.prop) <- c("celltype", "proportion")
  
  ##Create a slot in the object and save these cells under a name
  message(paste0(Sys.time(), ": Create a slot in object for these states"))
  obj@meta.data$cell.states <- NA
  obj@meta.data[long.term.cephalic, "cell.states"] <- "long-term"
  obj@meta.data[short.term.cephalic, "cell.states"] <- "short-term"
  
  nst <- as.data.frame(nearest)
  colnames(nst) <- c("stg.diff")
  message(paste0(Sys.time(), ": Plotting graphs for ", sample))
  pdf(file = paste0(plot.path, sample, "_eps_", eps, "_cell_states.pdf"), width = 28, height = 36)
  p <- nst %>% ggplot(aes(x = stg.diff)) + geom_histogram(binwidth = 3, fill = "#69b3a2", color="#e9ecef", alpha=0.9) + ggtitle("Mean_stage_difference")
  x <- DimPlot(obj, group.by = "stage.group")
  y <- long.prop %>% ggplot(aes(x = celltype, y = proportion)) + geom_bar(stat = "identity")
  z <- DimPlot(obj, group.by = "cell.states")
  gridExtra::grid.arrange(grobs = list(x, p, z, y), ncol = 2)
  dev.off()
}

nearest.cephalic <- readRDS(paste0(save.path, "cephalic_mean_stage_diff_eps_40.rds"))
nearest.cephalic <- as.data.frame(nearest.cephalic)
colnames(nearest.cephalic) <- "stage.diff"
nearest.all <- rbind(nearest.total, nearest.cephalic)

cells.out <- rownames(nearest.all)[which(rownames(nearest.all) %in% rownames(mama@meta.data))]
nearest.clean <- subset(nearest.all, cells = cells.extra)
nearest.clean <- as.data.frame(nearest.all[(rownames(nearest.all) %in% cells.out), ])
rownames(nearest.clean) <- cells.out
colnames(nearest.clean) <- "stage.diff"
cells.extra <- setdiff(rownames(nearest.all), rownames(mama@meta.data))

mama@meta.data <- cbind(mama@meta.data, nearest.clean[rownames(mama@meta.data), ])
colnames(mama@meta.data)[77] <- c("stage.diff")
##Get rid of all NaN values
mama@meta.data[which(mama@meta.data$stage.diff == "NaN"), "stage.diff"] <- 0
saveRDS(mama@meta.data$stage.diff, "~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_stage.diff_data_eps_40_with_neural.rds")
stage.diff <- readRDS("~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_stage.diff_data_eps_40_with_neural.rds")
mama@meta.data$stage.diff <- stage.diff

DimPlot(mama, group.by = "stage.diff", raster = F) + NoAxes() + NoLegend() + scale_fill_continuous(type = "gradient")
library(RColorBrewer)
pdf("~/Desktop/stage.diff.pdf")
DimPlot(mama, group.by = "stage.diff", raster = F) + NoAxes() + NoLegend() + scale_fill_continuous(type = "gradient")
dev.off()

##Categorize the different cells into how long they persist for using this eps setting
##Define cells that persist for different lengths of time
cells.under.24 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff < 24)]
cells.24.48 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff >= 24 & mama@meta.data$stage.diff < 48)]
cells.48.72 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff >= 48 & mama@meta.data$stage.diff < 72)]
#cells.36.48 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff >= 36 & mama@meta.data$stage.diff < 48)]
#cells.48.60 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff >= 48 & mama@meta.data$stage.diff < 60)]
#cells.60.72 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff >= 60 & mama@meta.data$stage.diff < 72)]
#cells.72.84 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff >= 72 & mama@meta.data$stage.diff < 84)]
cells.72 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff >= 72)]

##Now assign a slot in mama and plot these cells
mama@meta.data$long_vs_short <- NA
mama@meta.data[cells.under.24, "long_vs_short"] <- "<24hr"
mama@meta.data[cells.24.48, "long_vs_short"] <- "24-48hr"
mama@meta.data[cells.48.72, "long_vs_short"] <- "48-72hr"
mama@meta.data[cells.72, "long_vs_short"] <- ">72hr"

colors.stg <- setNames(c("#E9967A", "#1E90FF", "#00FF7F", "#FF1493"), c("<24hr", "24-48hr", "48-72hr",  ">72hr"))
dpi <- 300
png("~/Box/zfext/global_analysis_curated/transient_vs_long-term/figures/mama_transient_vs_persistent_states_data_eps_40_v3.png", height = 5*dpi, width = 5*dpi)
DimPlot(mama, group.by = "long_vs_short", cols = colors.stg, raster = F, na.value = "#f5f7f6") + NoAxes() + NoLegend()
dev.off()

saveRDS(mama@meta.data$long_vs_short, "~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_long_vs_short_stage_diff_data_eps_40.rds")
long_vs_short <- readRDS("~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_long_vs_short_stage_diff_data_eps_40.rds")
mama@meta.data$long_vs_short <- long_vs_short

##Load a dataframe with relevant information
df <- readRDS("~/Box/zfext/portaltest_static_celltype_curated/data/cyc/mama_tissue_cc_cell_class.rds")

##Assign cells which are dividing and not dividing
mama@meta.data <- cbind(mama@meta.data, df$cc.phase)
colnames(mama@meta.data)[68] <- "cc.phase"
DimPlot(mama, group.by = "cc.phase")
cells.non.cycling <- rownames(mama@meta.data)[which(is.na(mama@meta.data$cc.phase))]
cells.non.cycling <- rownames(mama@meta.data)[which(mama@meta.data$cc.phase == "non-cycling")]
cells.cycling <- rownames(mama@meta.data)[which(mama@meta.data$cc.phase == "G1/S-phase")]
cells.also.cycling <- rownames(mama@meta.data)[which(mama@meta.data$cc.phase == "G2M-phase")]
cells.non.cycling.total <- unlist(unique(list(cells.non.cycling, cells.also.non.cycling)))
mama@meta.data[cells.non.cycling.total, "cc.phase"] <- "non-cycling"
cells.cycling <- setdiff(rownames(mama@meta.data), cells.non.cycling)
cells.cycling.good <- setdiff(cells.cycling, cells.to.exclude)
cells.cycling.total <- unlist(unique(list(cells.cycling, cells.also.cycling)))
mama@meta.data$cc.status <- NA
mama@meta.data[cells.non.cycling, "cc.status"] <- "non-cycling"
mama@meta.data[cells.cycling.total, "cc.status"] <- "cycling"

color.stg.cc <- setNames(c("#4728E6", "#F2D9D0", "#3B5A57"), c("G1/S-phase", "G2M-phase", "non-cycling"))

png("~/Box/zfext/global_analysis_curated/transient_vs_long-term/figures/mama_cell_cycle_score_v2.png", width = 5*dpi, height = 5*dpi)
DimPlot(mama, group.by = "cc.phase", raster = F, cols = color.stg.cc) + NoAxes() + NoLegend()
dev.off()

##Load crude annotations and add that as a column to mama
curated.clusters <- readRDS("~/Box/zfext/02-Clustering/mama+ds_integrated/tissue_curated_clusters_mama_figure.rds")
mama@meta.data$curated.clusters <- curated.clusters
DimPlot(mama, group.by = "curated.clusters")

cells.in.clusters <- readRDS("~/Box/zfext/02-Clustering/mama+ds_integrated/mama_cluster_group_names.rds")
mama@meta.data$cells.in.clusters <- cells.in.clusters
DimPlot(mama, group.by = "cells.in.clusters")

df <- cbind(df, mama@meta.data$curated.clusters)
df <- cbind(df, mama@meta.data$cc.status)
colnames(df) <- c("stage.nice", "stage.group", "tissue.annot", "tissue.group", "tissue.stage", "clust", "clust.sg", "clust.stage", "s.score", "g2m.score", "cc.phase", "cells.class", "stage.diff", "long_vs_short", "annot", "cc.status")

##Divide into groups that includes cycling and non-cycling
mama@meta.data$stage.diff.cc <- NA
mama@meta.data$stage.diff.cc <- paste0(mama@meta.data$long_vs_short, "_", mama@meta.data$cc.status)
color.stg.cc <- setNames(c("#CCCCFF", "#E9967A", "#FF8000", "#FD1212", "#86B0DA", "#08519c", "#40E0D0", "#3B5A57", "#B2FF66", "#009900", "#CCCCFF", "#4728E6", "#FFCCFF", "#CC00CC", "#696969", "#404040"), c("<24hr_cycling", "<24hr_non-cycling", "24-48hr_cycling", "24-48hr_non-cycling", "48-72hr_cycling", "48-72hr_non-cycling", ">72hr_cycling", ">72hr_non-cycling"))

saveRDS(mama@meta.data$stage.diff.cc, "~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_stage_diff_cc_status_data_eps_40.rds")

##Just plot the cycling states
color.stg.cycling <- setNames(c("#E9967A", "#1E90FF", "#00FF7F", "#FF1493"), c("<24hr_cycling", "24-48hr_cycling", "48-72hr_cycling", ">72hr_cycling"))

color.stg.non.cycling <- setNames(c("#E9967A", "#1E90FF", "#00FF7F", "#FF1493"), c("<24hr_non-cycling", "24-48hr_non-cycling", "48-72hr_non-cycling", ">72hr_non-cycling"))

dpi <- 300
png("~/Box/zfext/global_analysis_curated/transient_vs_long-term/figures/mama_transient_vs_persistent_states_cc.state_eps_40_v3.png", height = 5*dpi, width = 5*dpi)
DimPlot(mama, group.by = "stage.diff.cc", cols = color.stg.cc, raster = F, na.value = "#f5f7f6") + NoAxes()
dev.off()

dpi <- 300
png("~/Box/zfext/global_analysis_curated/transient_vs_long-term/figures/mama_transient_vs_persistent_states_only_cycling_eps_35_v3.png", height = 5*dpi, width = 5*dpi)
DimPlot(mama, group.by = "stage.diff.cc", cols = color.stg.cycling, raster = F, na.value = "#f5f7f6") + NoAxes() + NoLegend()
dev.off()

dpi <- 300
png("~/Box/zfext/global_analysis_curated/transient_vs_long-term/figures/mama_transient_vs_persistent_states_only_non-cycling_eps_35_v2.png", height = 5*dpi, width = 5*dpi)
DimPlot(mama, group.by = "stage.diff.cc", cols = color.stg.non.cycling, raster = F, na.value = "#f5f7f6") + NoAxes() + NoLegend()
dev.off()

##Check what clusters belong to these categories
cells.24h.non.cycling <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff.cc == "<24hr_non-cycling")]
table(mama@meta.data[cells.24h.non.cycling, "clust"])


########### PLOTTING ##################
##Create timeline plots for how long various stages are persistent
##Consider plotting separate timeline plots for cycling and non-cycling populations
##Subset a dataframe from mama metadata
##Load in the precalculated cell cycle score on mama
cc.score <- readRDS("~/Box/zfext/02-Clustering/mama+ds_integrated/merged_mama_cellCycle_score.rds")
mama@meta.data <- cbind(mama@meta.data, cc.score)
dat <- mama@meta.data[,c("stage.nice", "stage.group", "tissue.stage", "tissue.subsets", "clust", "clust.sg", "s.score", "g2m.score", "cc.phase", "long_vs_short", "curated.clusters", "cc.status", "stage.diff.cc")]
p <- DimPlot(mama, group.by = "stage.diff.cc", cols = color.stg.cc, raster = F) + NoAxes()
color.stg.cc <- unique(ggplot_build(p)$data[[1]]$colour)

# Load cluster annotations 
is.empty <- function(x) {
  is.na(x) | is.null(x) | (sapply(trimws(x), nchar) == 0)
}
annot.confirmed <- readxl::read_xlsx("~/Box/zfext/02-Clustering/mama_tissue_clustering/annot/mama_annotations_full_extended_jeffedit_newtissues_curated_clustering.xlsx")
annot.confirmed <- annot.confirmed[!is.empty(annot.confirmed$clust),3:ncol(annot.confirmed)]
rownames(annot.confirmed) <- annot.confirmed$
annot.confirmed$confirmed <- T


##Merge the annotations to the dataframe
colnames(annot.confirmed$UI) <- "clust"
dat <- merge(dat, annot.confirmed[, c("clust", "identity.super")], by = "clust")

dat$color <- NA
dat[which(dat$stage.diff.cc == "<12hr_cycling"), "color"] <- "#F2D9D0"
dat[which(dat$stage.diff.cc == "<12hr_non-cycling"), "color"] <- "#E9967A"
dat[which(dat$stage.diff.cc == "12-24hr_cycling"), "color"] <- "#800000"
dat[which(dat$stage.diff.cc == "12-24hr_non-cycling"), "color"] <- "#FD1212"
dat[which(dat$stage.diff.cc == "24-36hr_cycling"), "color"] <- "#86B0DA"
dat[which(dat$stage.diff.cc == "24-36hr_non-cycling"), "color"] <- "#08519c"
dat[which(dat$stage.diff.cc == "36-48hr_cycling"), "color"] <- "#40E0D0"
dat[which(dat$stage.diff.cc == "36-48hr_non-cycling"), "color"] <- "#3B5A57"
dat[which(dat$stage.diff.cc == "48-60hr_cycling"), "color"] <- "#4DEC53"
dat[which(dat$stage.diff.cc == "48-60hr_non-cycling"), "color"] <- "#054707"
dat[which(dat$stage.diff.cc == "60-72hr_cycling"), "color"] <- "#5b4e9f"
dat[which(dat$stage.diff.cc == "60-72hr_non-cycling"), "color"] <- "#4728E6"
dat[which(dat$stage.diff.cc == "72-84hr_cycling"), "color"] <- "#94699F"
dat[which(dat$stage.diff.cc == "72-84hr_non-cycling"), "color"] <- "#C023E8"
dat[which(dat$stage.diff.cc == ">84hr_cycling"), "color"] <- "#2F4F4F"
dat[which(dat$stage.diff.cc == ">84hr_non-cycling"), "color"] <- "#696969"

pdf("~/Box/zfext/global_analysis_curated/transient_vs_long-term/figures/mama_cell_states_segment_plot_cc.pdf", height = 30, width = 26)
ggplot(dat, aes(x=0, xend=stage.diff, y=reorder(identity.super, stage.diff), yend=reorder(identity.super, stage.diff), color = stage.diff.cc)) + geom_text(aes(x=stage.diff, label = NA, hjust=-0.3), family="sans") +
  geom_segment(size=4, aes(colour=color)) +
  scale_colour_identity() +
  theme(axis.text=element_text(size=10,family="sans"),
        axis.title=element_text(size=13,face="bold",family="sans"),
        strip.text.y = element_text(size=12,family="sans"), 
        plot.title=element_text(size=14,face="bold",family="sans")) +      
  facet_grid(tissue.annot~., scales = "free") + labs(x = "hpf",
                                                     y = "cell-types",
                                                     color = "stage.diff.cc") +
  scale_fill_manual(name = NULL,
                    values = c(
                      "<12hr_cycling" =   "#F2D9D0",
                      "<12hr_non-cycling" = "#E9967A",
                      "12-24hr_cycling" = "#800000",
                      "12-24hr_non-cycling" = "#FD1212",
                      "24-36hr_cycling" = "#86B0DA",
                      "24-36hr_non-cycling" = "#08519c",
                      "36-48hr_cycling" = "#40E0D0" ,
                      "36-48hr_non-cycling" = "#3B5A57" ,
                      "48-60hr_cycling" = "#4DEC53",
                      "48-60hr_non-cycling" = "#054707",
                      "60-72hr_cycling" = "#5b4e9f",
                      "60-72hr_non-cycling" = "#4728E6",
                      "72-84hr_cycling" = "#94699F",
                      "72-84hr_non-cycling" = "#C023E8"
                    ),
                    breaks = unique(dat$stage.diff.cc)
  ) +
  theme(axis.text.x = element_text(color = "gray60", size = 10)) +
  theme(legend.position = "top") +
  theme_bw()

dev.off()

##Calculate percentage of cells in these different categories
table(mama@meta.data$stage.diff.cc)
table(mama@meta.data[rownames(mama@meta.data)[which(mama@meta.data$stage.diff.cc == "24-36hr_non-cycling")], "identity.super"])

mama@meta.data$identity.super <- NA
for(i in annot.confirmed$clust){
  mama@meta.data[which(mama@meta.data$clust == i), "identity.super"] <- annot.confirmed[which(annot.confirmed$clust == i), "identity.super"]
}

df <- as.data.frame(table(mama@meta.data$stage.diff.cc))
colnames(df) <- c("category", "cell_no.")

# Barplot
ggplot(df, aes(x=category, y=cell_no.)) + 
  geom_bar(stat = "identity")


##Do a plot only for dividing cells
##Subset dataframe to only include cycling cells
dt <- dat[which(dat$cc.status == "cycling"), ]

pdf("~/Box/zfext/global_analysis_curated/transient_vs_long-term/figures/mama_cell_states_segment_plot_cycling.pdf", height = 30, width = 26)
ggplot(dt, aes(x=0, xend=stage.diff, y=reorder(identity.super, stage.diff), yend=reorder(identity.super, stage.diff), color = stage.diff.cc)) + geom_text(aes(x=stage.diff, label = NA, hjust=-0.3), family="sans") +
  geom_segment(size=4, aes(colour=color)) +
  scale_colour_identity() +
  theme(axis.text=element_text(size=10,family="sans"),
        axis.title=element_text(size=13,face="bold",family="sans"),
        strip.text.y = element_text(size=12,family="sans"), 
        plot.title=element_text(size=14,face="bold",family="sans")) +      
  facet_grid(tissue.annot~., scales = "free") + labs(x = "hpf",
                                                     y = "cell-types",
                                                     color = "stage.diff.cc") +
  scale_fill_manual(name = NULL,
                    values = c(
                      "<12hr_cycling" =   "#F2D9D0",
                      "<12hr_non-cycling" = "#E9967A",
                      "12-24hr_cycling" = "#800000",
                      "12-24hr_non-cycling" = "#FD1212",
                      "24-36hr_cycling" = "#86B0DA",
                      "24-36hr_non-cycling" = "#08519c",
                      "36-48hr_cycling" = "#40E0D0" ,
                      "36-48hr_non-cycling" = "#3B5A57" ,
                      "48-60hr_cycling" = "#4DEC53",
                      "48-60hr_non-cycling" = "#054707",
                      "60-72hr_cycling" = "#5b4e9f",
                      "60-72hr_non-cycling" = "#4728E6",
                      "72-84hr_cycling" = "#94699F",
                      "72-84hr_non-cycling" = "#C023E8"
                    ),
                    breaks = unique(dat$stage.diff.cc)
  ) +
  theme(axis.text.x = element_text(color = "gray60", size = 10)) +
  theme(legend.position = "top") +
  theme_bw()

dev.off()

##Ok so we have a dataframe with the stage difference information of every cluster at various stages
##Make a segment plot - hpf (stage) on the x-axis and tissue names on the Y-axis

##Ok subset dataset to one tissue at a time which contains long-term proliferative states
df <- dat[which(dat$tissue.annot == "hematopoietic"), ]
df.cyc <- df[which(df$cc.phase != "non-cycling"), ]

dx <- df.cyc[, c("clust", "long_vs_short", "stage.diff.cc", "identity.super")]
dx <- unique(dx)
dx$stage.start <- NA
dx$stage.stop <- NA

for(i in dx$clust){
      for(j in dx$stage.diff.cc){
      mt <- df.cyc[which(df.cyc$clust == i & df.cyc$stage.diff.cc == j), ]
      dx[which(dx$clust == i & dx$stage.diff.cc == j), "stage.start"] <- min(mt$stage.nice)
      }
  }

for(i in dx$clust){
  for(j in dx$stage.diff.cc){
    mt <- df.cyc[which(df.cyc$clust == i & df.cyc$stage.diff.cc == j), ]
    dx[which(dx$clust == i & dx$stage.diff.cc == j), "stage.stop"] <- max(mt$stage.nice)
  }
}

dm <- dx[which(dx$clust == "hema.6"), ]

library(ggplot2)
pdf("hematopoietic_cell_states_time_lengths.pdf", width = 28, height = 36)
ggplot(data=dx, aes(color=stage.diff.cc))+
  geom_segment(aes(x=stage.start, xend=stage.stop, y=stage.diff.cc, yend = stage.diff.cc),lwd=4)+
  facet_grid(clust~.)+xlab("HPF")
dev.off()


##Muscle_superset
df <- dat[which(dat$tissue.annot == "muscle"), ]
df.cyc <- df[which(df$cc.phase != "non-cycling"), ]

dx <- df.cyc[, c("clust", "long_vs_short", "stage.diff.cc", "identity.super")]
dx <- unique(dx)
dx$stage.start <- NA
dx$stage.stop <- NA

for(i in dx$clust){
  for(j in dx$stage.diff.cc){
    mt <- df.cyc[which(df.cyc$clust == i & df.cyc$stage.diff.cc == j), ]
    dx[which(dx$clust == i & dx$stage.diff.cc == j), "stage.start"] <- min(mt$stage.nice)
  }
}

for(i in dx$clust){
  for(j in dx$stage.diff.cc){
    mt <- df.cyc[which(df.cyc$clust == i & df.cyc$stage.diff.cc == j), ]
    dx[which(dx$clust == i & dx$stage.diff.cc == j), "stage.stop"] <- max(mt$stage.nice)
  }
}

dm <- dx[which(dx$clust == "hema.6"), ]

library(ggplot2)
pdf("muscle_cell_states_time_lengths.pdf", width = 28, height = 36)
ggplot(data=dx, aes(color=stage.diff.cc))+
  geom_segment(aes(x=stage.start, xend=stage.stop, y=stage.diff.cc, yend = stage.diff.cc),lwd=4)+
  facet_grid(clust~.)+xlab("HPF")
dev.off()


##Muscle_superset
df <- dat
df.nocyc <- df[which(df$cc.phase == "non-cycling"), ]

dx <- df.nocyc[, c("clust", "long_vs_short", "stage.diff.cc", "identity.super")]
dx <- unique(dx)
dx$stage.start <- NA
dx$stage.stop <- NA

for(i in dx$clust){
  for(j in dx$stage.diff.cc){
    mt <- df.nocyc[which(df.nocyc$clust == i & df.nocyc$stage.diff.cc == j), ]
    dx[which(dx$clust == i & dx$stage.diff.cc == j), "stage.start"] <- min(mt$stage.nice)
  }
}

for(i in dx$clust){
  for(j in dx$stage.diff.cc){
    mt <- df.nocyc[which(df.nocyc$clust == i & df.nocyc$stage.diff.cc == j), ]
    dx[which(dx$clust == i & dx$stage.diff.cc == j), "stage.stop"] <- max(mt$stage.nice)
  }
}

dm <- dx[which(dx$clust == "hema.6"), ]

library(ggplot2)
pdf("muscle_cell_states_time_lengths.pdf", width = 28, height = 36)
ggplot(data=dx, aes(color=stage.diff.cc))+
  geom_segment(aes(x=stage.start, xend=stage.stop, y=stage.diff.cc, yend = stage.diff.cc),lwd=4)+
  facet_grid(clust~.)+xlab("HPF")
dev.off()


##Check all long-term states in the eye
dat <- as.data.frame(table(mama@meta.data[long.cycling.cells, "clust"]))
colnames(dat) <- c("clust", "long_cyc")
dat$tissue <- unlist(lapply(strsplit(x = as.character(dat$clust), split = "\\."), function(x) x[1]))
colnames(dat) <- c("clust", "long_cyc", "tissue")

dm <- as.data.frame(table(mama@meta.data[long.non_cycling.cells, "clust"]))
colnames(dm) <- c("clust", "long_non-cyc")
dm$tissue <- unlist(lapply(strsplit(x = as.character(dm$clust), split = "\\."), function(x) x[1]))

dx <- merge(dat, dm)
dx$long.ratio <- dx$long_cyc/dx$`long_non-cyc`
dx[which(dx$long.ratio > 1), "clust"]

cells.48.72 <- rownames(mama@meta.data)[which(mama@meta.data$long_vs_short == "48-72hr")]
cells.above.72 <- rownames(mama@meta.data)[which(mama@meta.data$long_vs_short == ">72hr")]
cells.above.48.total <- unlist(unique(list(cells.48.72, cells.above.72)))

t1 <- table(mama@meta.data[cells.above.48.total, "stage.nice"])

sum(t1[as.numeric(names(t1)) >= 48])
sum(t1[as.numeric(names(t1)) >= 72])
sum(t1[as.numeric(names(t1)) >= 14 & as.numeric(names(t1)) <= 28])

t2 <- table(mama@meta.data[long.cycling.cells, "stage.nice"])
sum(t2[as.numeric(names(t2)) >= 48])
sum(t2[as.numeric(names(t2)) >= 72])
sum(t2[as.numeric(names(t2)) >= 14 & as.numeric(names(t2)) <= 28])

t3 <- table(mama@meta.data[cells.non.cycling, "stage.nice"])
sum(t3[as.numeric(names(t3)) >= 72])

t4 <- table(mama@meta.data$stage.group)


library(dplyr)
data <- master.df  %>%
  group_by(stage, tissue) %>%
  summarise(n = sum(value)) %>%
  mutate(percentage = n / sum(n))

##Periderm and erythroblasts
cells.48.72.above <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff.cc == "48-72hr_cycling")]
cells.72.above <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff.cc == ">72hr_cycling")]
cells.long <- unlist(unique(list(cells.48.72.above, cells.72.above)))
cells.late <- rownames(mama@meta.data)[which(mama@meta.data$stage.group == "120")]

cells.in.late <- cells.long[cells.long %in% cells.late]

Idents(mama) <- mama@meta.data$clust




