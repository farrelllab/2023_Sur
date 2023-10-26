##Cluster the merged dropseq and 10X dataset to generate the Figures 1B, 1C

##Load libraries
library(Seurat)
library(URD)
library(Matrix)

##Load mama
mama <- readRDS("~/Box/zfext/02-Clustering/mama+ds_integrated/merged_mama_ds_slim_seurat_25Jan2022.rds")

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

##Plot the stage plots as in Figure 1B

##Set colors for stages
stage.colors.new <- c(
  colorRampPalette(c("#f4be1d", "#efb918", "#bd8700", "#b37d00", "#a97300", "#9f6900", "#955f00", "#8b5500"))(8), # 3, 4, 5, 6, 7, 8, 9, 10
  colorRampPalette(c("#FFC0CB", "#FFB6C1", "#FF69B4", "#DA70D6", "#BA55D3", "#C71585", "#FF1493"))(7), # 11, 12, 14, 16, 18, 21, 24
  rep(RColorBrewer::brewer.pal(6, "Oranges"), each = 2), # 26 / 28, 30 / 32, 34 / 36, 38 / 40, 42 / 44, 46
  rep(RColorBrewer::brewer.pal(6, "Greens"), each = 2), # 48, 50 / 52, 54 / 56, 58 / 60, 62 / 64, 66 / 68, 70,  72
  rep(RColorBrewer::brewer.pal(6, "Blues"), each = 2),  # 74 / 76, 78 / 80, 82 / 84, 86 / 88, 90 / 92, 94, 96
  rep(RColorBrewer::brewer.pal(7, "Purples"), each = 2) #  98 / 100, 102 / 104, 106 / 108, 110 / 112, 114 / 116, 118 / 120
)

##Plot the UMAP as in Figure 1B
dpi <- 300
png("~/Desktop/mama_stage.png", width = 5*dpi, height = 5*dpi)
DimPlot(mama, group.by = "stage.nice", cols = stage.colors.new, raster = F)
dev.off()

##Plot the UMAP as in Figure 1C
dpi <- 300
png("~/Desktop/mama_tissue_subsets.png", width = 5*dpi, height = 5*dpi)
DimPlot(mama, group.by = "tissue_subsets", cols = stage.colors.new, raster = F)
dev.off()

##Color the UMAPs by individual groups of stages as shown in Figure S1A
colors.stg.3.10 <- setNames(c("#708090", "#1A85FF", "#0000FF", "#2b8cbe", "#17A589", "#00FF00", "#f4be1d", "#8b5500"), unique(mama@meta.data$stage.nice)[1:8])
colors.stg.12.21 <- setNames(c("#708090", "#1A85FF", "#0000FF", "#2b8cbe", "#17A589", "#00FF00", "#f4be1d", "#8b5500"), unique(mama@meta.data$stage.nice)[9:15])
colors.stg.24.46 <- setNames(c("#708090", "#1A85FF", "#0000FF", "#2b8cbe", "#17A589", "#00FF00", "#f4be1d", "#8b5500", "#FF7F00", "#FF1493", "#FF00FF", "#8B008B"), unique(mama@meta.data$stage.nice)[16:27])
colors.stg.48.70 <- setNames(c("#708090", "#1A85FF", "#0000FF", "#2b8cbe", "#17A589", "#00FF00", "#f4be1d", "#8b5500", "#FF7F00", "#FF1493", "#FF00FF", "#8B008B"), unique(mama@meta.data$stage.nice)[28:39])
colors.stg.72.94 <- setNames(c("#708090", "#1A85FF", "#0000FF", "#2b8cbe", "#17A589", "#00FF00", "#f4be1d", "#8b5500", "#FF7F00", "#FF1493", "#FF00FF", "#8B008B", "#F08080", "#FF0000"), unique(mama@meta.data$stage.nice)[40:51])
colors.stg.96.120 <- setNames(c("#708090", "#1A85FF", "#0000FF", "#2b8cbe", "#17A589", "#00FF00", "#f4be1d", "#8b5500", "#FF7F00", "#FF1493", "#FF00FF", "#8B008B", "#F08080", "#FF0000"), unique(mama@meta.data$stage.nice)[52:63])

dpi <- 300
png("~/Box/zfext/results/04c-SubsetsV3/figures_for_manuscript/Figure_1_mama_cell_atlas/mama_stage_3-10hpf_v2.png", width = 5*dpi, height = 5*dpi)
DimPlot(mama, cols = colors.stg.3.10, raster = F, na.value = "#f5f7f6", order = 3:10) + NoLegend() + NoAxes()
dev.off()

png("~/Box/zfext/results/04c-SubsetsV3/figures_for_manuscript/Figure_1_mama_cell_atlas/mama_stage_12-21hpf_v2.png", width = 5*dpi, height = 5*dpi)
DimPlot(mama, cols = colors.stg.12.21, raster = F, na.value = "#f5f7f6", order = 12:21) + NoLegend() + NoAxes()
dev.off()

png("~/Box/zfext/results/04c-SubsetsV3/figures_for_manuscript/Figure_1_mama_cell_atlas/mama_stage_24-46hpf_v2.png", width = 5*dpi, height = 5*dpi)
DimPlot(mama, cols = colors.stg.24.46, raster = F, na.value = "#f5f7f6", order = 24:46) + NoLegend() + NoAxes()
dev.off()

png("~/Box/zfext/results/04c-SubsetsV3/figures_for_manuscript/Figure_1_mama_cell_atlas/mama_stage_48-70hpf_v2.png", width = 5*dpi, height = 5*dpi)
DimPlot(mama, cols = colors.stg.48.70, raster = F, na.value = "#f5f7f6", order = 48:70) + NoLegend() + NoAxes()
dev.off()

png("~/Box/zfext/results/04c-SubsetsV3/figures_for_manuscript/Figure_1_mama_cell_atlas/mama_stage_72-94hpf_v2.png", width = 5*dpi, height = 5*dpi)
DimPlot(mama, cols = colors.stg.72.94, raster = F, na.value = "#f5f7f6", order = 72:94) + NoLegend() + NoAxes()
dev.off()

png("~/Box/zfext/results/04c-SubsetsV3/figures_for_manuscript/Figure_1_mama_cell_atlas/mama_stage_96-120hpf_v2.png", width = 5*dpi, height = 5*dpi)
DimPlot(mama, cols = colors.stg.96.120, raster = F, na.value = "#f5f7f6", order = 96:120) + NoLegend() + NoAxes()
dev.off()

##Use ggplot2 to plot all cells in mama and color the relevant stages by their colors
library(ggplot2)

cells.3.10 <- rownames(mama@meta.data)[which(mama@meta.data$stage.nice %in% c(3:10))]
cells.12.21 <- rownames(mama@meta.data)[which(mama@meta.data$stage.nice %in% c(12:21))]
cells.24.46 <- rownames(mama@meta.data)[which(mama@meta.data$stage.nice %in% c(24:46))]
cells.48.70 <- rownames(mama@meta.data)[which(mama@meta.data$stage.nice %in% c(48:70))]
cells.72.94 <- rownames(mama@meta.data)[which(mama@meta.data$stage.nice %in% c(72:94))]
cells.96.120 <- rownames(mama@meta.data)[which(mama@meta.data$stage.nice %in% c(96:120))]

# Assuming umap_df is your dataframe and cell_group is the column defining groups
umap_df <- as.data.frame(mama@reductions$umap@cell.embeddings)
umap_df$stage.nice <- mama@meta.data$stage.nice
#for(i in rownames(umap_df)){
#  umap_df[i, "stage.nice"] <- mama@meta.data[i, "stage.nice"]
#}

# Plot all cells in gray
p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(color = "#C0C0C0", size = 1)

# Overlay specific groups with different colors
overlay_colors <- c("#708090", "#1A85FF", "#0000FF", "#2b8cbe", "#17A589", "#00FF00", "#f4be1d", "#8b5500", "#FF7F00", "#FF1493", "#FF00FF", "#8B008B")  # Define colors for different groups
groups_to_overlay <- c(unique(umap_df$stage.nice)[52:63]) # Define groups to overlay

for (i in seq_along(groups_to_overlay)) {
  group <- groups_to_overlay[i]
  color <- overlay_colors[i]
  
  overlay_points <- umap_df[umap_df$stage.nice == group, ]
  p <- p + geom_point(data = overlay_points, aes(x = UMAP_1, y = UMAP_2), color = color, size = 1.5) + theme_minimal() +  # Choose a minimal theme
    theme(panel.grid = element_blank(),  # Remove grid lines
          panel.background = element_rect(fill = "white"))  # Set white background
}

# Print the plot
dpi <- 300
png("~/Box/zfext/results/04c-SubsetsV3/figures_for_manuscript/Figure_1_mama_cell_atlas/mama_stage_96-120hpf_v4.png", width = 5*dpi, height = 5*dpi)
p
dev.off()


##Clustering the global data to show the 19 major tissues as in Figure S1H and in Daniocell
dpi <- 300
png("~/Desktop/mama_tissue_categories.png", width = 5*dpi, height = 5*dpi)
DimPlot(mama, group.by = "tissue_categories", cols = stage.colors.new, raster = F)
dev.off()




