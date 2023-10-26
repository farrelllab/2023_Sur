library(Seurat)
library(URD)

#Grabbing smooth muscle cells from different sources:
#Smooth muscle cell markers: acta2/tagln, myl9a, lmod1b, cnn1b, desmb etc. 
#Let's see where these genes are expressed?

##Read in new mural cell object
obj.smc <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/mural/obj_seurat/mural_seurat.rds")
markers <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/mural/obj_seurat/mural_markers.rds")

##Plot  stage information on the non-skeletal muscle atlas UMAP - Figure S4A
stage.colors.new <- c(
  colorRampPalette(c("#BA55D3", "#C71585", "#FF1493"))(7), # 14, 16, 18, 21`
  rep(RColorBrewer::brewer.pal(6, "Oranges"), each = 2), # 24, 26 / 28, 30 / 32, 34 / 36, 38 / 40, 42 / 44, 46
  rep(RColorBrewer::brewer.pal(6, "Greens"), each = 2), # 48, 50 / 52, 54 / 56, 58 / 60, 62 / 64, 66 / 68, 70
  rep(RColorBrewer::brewer.pal(6, "Blues"), each = 2),  # 72, 74 / 76, 78 / 80, 82 / 84, 86 / 88, 90 / 92, 94
  rep(RColorBrewer::brewer.pal(7, "Purples"), each = 2) # 96, 98 / 100, 102 / 104, 106 / 108, 110 / 112, 114 / 116, 118 / 120
)

##Define new color scales as the reviewer suggested
stage.colors <- setNames(c("#0000FF", "#2b8cbe", "#17A589", "#00FF00", "#f4be1d", "#8b5500", "#FF7F00", "#FF1493", "#FF00FF", "#8B008B"), c(unique(obj.smc@meta.data$stage.group)))

##Plot UMAP as in Figure S4A
dpi <- 300
png(file = "mural_cells_UMAP_stage_colored.png", width = 5 * dpi, height = 5 * dpi)
DimPlot(obj.smc, group.by = "stage.nice", cols = stage.colors, pt.size = 2)
dev.off()

DimPlot(obj.smc, label = T)
DimPlot(obj.smc, group.by = "stage.group")
DimPlot(obj.smc, group.by = "RNA_snn_res.3")

##Ok, RNA_snn_res.3 was a better representation of the clustering of these cells
Idents(obj.smc) <- obj.smc@meta.data$RNA_snn_res.3

##Calculate differentially expressed markers between all clusters of the non-skeletal muscle atlas
posMarkers.wilcox <- sapply(levels(Idents(obj)), function(i) FindMarkers(obj,i,assay.type='RNA',test.use="wilcox",only.pos=T, min.cells.group = 1),simplify=F)
posMarkers.roc <- sapply(levels(Idents(obj)), function(i) FindMarkers(obj,i,test.use="roc",only.pos=T),simplify=F)
markers <- list(wilcox = posMarkers.wilcox, roc = posMarkers.roc)
saveRDS(markers, file="~/Box/zfext/annotations_celltype_curated_newMama/mural_cells/obj_seurat/mural-cells_markers.rds")
saveRDS(obj, file="~/Box/zfext/annotations_celltype_curated_newMama/mural_cells/obj_seurat/mural-cells_seurat.rds")

## I. How many pericyte populations do we detect? Can we distinguish them transcriptionally?
##Ok, there are several pieces in the atlas that express the pericyte marker ndufa4l2a. So get those cells out and use them as a separate cluster.

##Plotting the UMAP with clusters - Figure 3A
dpi <- 300
png(file = "mural_cells_UMAP_clusters.png", width = 5 * dpi, height = 5 * dpi)
DimPlot(obj.smc, pt.size = 2) + NoLegend() + NoAxes()
dev.off()


##Plotting some genes on the UMAP projection
FeaturePlot(obj.smc, features = c("ndufa4l2a", "il13ra2", "tesca", "abcc9"))
FeaturePlot(obj.smc, features = c("bgna", "cnn1a", "cln6b", "tdrd9"))
FeaturePlot(obj.smc, features = c("bgna", "ndufa4l2a", "kcnk18", "smpx"))

cluster.names.new <- c("vSMC-aorta", "pericyte_3", "pericyte_2", "mesenchyme_int", "cycling", "neural", "cycling", "int-SMCs_circular", "int-SMCs_longitudinal", "vSMC-transient", "pericyte-transient", "cycling", "vi-SMC", "epicardium", "vSMC-OFT", "HSCs", "cardiac_muscle", "myofibroblasts", "myocardium", "neural", "mesenchyme_int", "HSCs", "pericytes-brain", "vSMC-coro", "neural", "SMC-post", "NCC-precursors", "fibroblasts", "epidermis")
names(cluster.names.new) <- levels(obj.smc)
obj.sm <- RenameIdents(obj.smc, cluster.names.new) 


##II. Compare between the three major classes of non-skeletal muscle cells captured in our atlas - pericytes, vascular, and visceral smooth muscles
##Get cells belonging to  these categories

cells.vascular <- WhichCells(obj.sm, idents = c("2", "15", "11", "21"))
cells.visceral <- WhichCells(obj.smc, idents = c("13", "10", "8", "16", "22", "23"))
cells.pericytes <- WhichCells(obj.smc, idents = c("20", "4", "9"))

obj.sm <- subset(obj.smc, cells = unlist(unique(list(cells.vascular, cells.visceral, cells.pericytes))))

obj.sm@meta.data$smc.class <- NA
obj.sm@meta.data[cells.vascular, "smc.class"] <- "vascular"
obj.sm@meta.data[cells.visceral, "smc.class"] <- "visceral"
obj.sm@meta.data[cells.pericytes, "smc.class"] <- "pericytes"

colors <- setNames(c("#71B000", "#FF6C91", "#00B7E8"), c("vascular", "visceral", "pericytes"))
DimPlot(obj.sm, group.by = "smc.class", cols = colors)

Idents(obj.sm) <- obj.sm@meta.data$smc.class
my_levels <- c("vSMC-transient", "vSMC_aorta", "vSMC-OFT", "vSMC-artery_1", "vSMC-artery_2", "vi-SMCs", "int-SMCs_circular", "int-SMCs_longitudinal", "pericyte_transient", "pericytes-brain", "pericyte_2", "pericyte_3", "NCC-precursors", "mesenchyme_int", "cycling", "unknown", "cardiac_muscle", "myocardium", "epicardium", "HSCs", "fibroblasts", "myofibroblasts?", "SMC-post")
# Relevel object@ident
Idents(obj.sm) <- factor(x = Idents(obj.sm), levels = my_levels)

##Set the smc.class slot as the idents for the SMC object
Idents(obj.sm) <- obj.sm@meta.data$smc.class

genes.smc.peri <- FindAllMarkers(obj.sm, assay = "RNA", only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
library(dplyr)
top15 <- genes.smc.peri %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
gene.list <- unlist(top20$gene)

filter.genes <- function(genes) {
  mt.genes <- grep("^mt-", ignore.case = T, genes, value = T)
  many.genes <- grep("\\(1 of many\\)", ignore.case = T, genes, value = T)
  cox.genes <- grep("^cox", ignore.case = T, genes, value = T)
  rp.genes <- grep("^rp", ignore.case = T, genes, value = T)
  hb.genes <- grep("^hb", ignore.case = T, genes, value = T)
  si.genes <- grep("^si:", ignore.case = T, genes, value = T)
  zgc.genes <- grep("^zgc:", ignore.case = T, genes, value = T)
  loc.genes <- grep("^LOC", ignore.case = T, genes, value = T)
  bx.genes <- grep("^BX", ignore.case = T, genes, value = T)
  cr.genes <- grep("^CT", ignore.case = T, genes, value = T)
  return(setdiff(genes, c(mt.genes, many.genes, cox.genes, rp.genes, hb.genes, si.genes, zgc.genes, loc.genes, cr.genes, bx.genes)))
}

gene.list <- filter.genes(gene.list)
gene.list <- c("wfdc2", "fbln5", "fbp2", "ldha", "slc16a3", "loxa", "desmb", "cald1b", "csrp1b", "nkx2.3", "smtna", "smtnb", "abcc9", "ndufa4l2a", "pdgfrb", "bgnb", "rgs5b", "foxl1")

##Plot the comparative dotolot between pericytes, vascular SMCs and visceral SMCs - Figure 4B
pdf("pericyte_vaSMC_viSMC_dotplot.pdf", width = 8, height = 10)
DotPlot(obj.sm, features = gene.list, dot.scale = 10, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()


## III. Ok, is there a fourth population of pericytes that is at the tip of cluster 21 of the newly reclustered mural-cell atlas?
##In that case, let's take those cells out and compare them to the other pericyte populations
##How distinct are these? I know that these also express acta2 along with ndufa4l2a which is a cell state that has never been observed before. 
## Generally, ndufa4l2a is known to be a pericyte-specific markers while acta2 is known to be a SMC-specific marker


##Create a sub-object with all pericyte populations across the atlas
##Are there three pericyte populations?
##Let's compare the three pericyte populations and the transient state - that's becoz the early pericyte markers are also not known
cells.peri.1 <- WhichCells(obj.smc, idents = c("9"), expression = ndufa4l2a > 1)
cells.peri.2 <- WhichCells(obj.smc, idents = c("20"))
cells.peri.3 <- WhichCells(obj.smc, idents = c("4"))
cells.peri.4 <- WhichCells(obj.smc, idents = c("3"))

cells.all <- unlist(unique(list(cells.peri.1, cells.peri.2, cells.peri.3, cells.peri.4)))

##Create a pericyte sub object in Seurat
obj.peri <- subset(obj.smc, cells = cells.all)

##Assign the cells from the different clusters to different names within a slot in the object - "pericyte.clusters"
obj.peri@meta.data$pericyte.clusters <- NA
obj.peri@meta.data[cells.peri.1, "pericyte.clusters"] <- "pericyte-0"
obj.peri@meta.data[cells.peri.2, "pericyte.clusters"] <- "pericyte_1"
obj.peri@meta.data[cells.peri.3, "pericyte.clusters"] <- "pericytes_2"
obj.peri@meta.data[cells.peri.4, "pericyte.clusters"] <- "mFB"

##Add pericyte.clusters slot as  the idents to the pericyte seurat object
Idents(obj.peri) <- obj.peri@meta.data$pericyte.clusters

library(dplyr)
markers.peri <- FindAllMarkers(obj.peri, assay = "RNA", logfc.threshold = 0.25, only.pos = T, min.pct = 0.25)
top20 <- markers.peri %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
gene.list <- unlist(unique(list(top20$gene, "pdgfrb", "ndufa4l2a", "abcc9", "kcnj8", "plp1b", "bgna", "pld6", "cnn1a", "pthlha", "fbp2", "tdrd9", "cln6b")))

filter.genes <- function(genes) {
  mt.genes <- grep("^mt-", ignore.case = T, genes, value = T)
  many.genes <- grep("\\(1 of many\\)", ignore.case = T, genes, value = T)
  cox.genes <- grep("^cox", ignore.case = T, genes, value = T)
  rp.genes <- grep("^rp", ignore.case = T, genes, value = T)
  hb.genes <- grep("^hb", ignore.case = T, genes, value = T)
  si.genes <- grep("^si:", ignore.case = T, genes, value = T)
  zgc.genes <- grep("^zgc:", ignore.case = T, genes, value = T)
  loc.genes <- grep("^XLOC", ignore.case = T, genes, value = T)
  bx.genes <- grep("^BX", ignore.case = T, genes, value = T)
  cr.genes <- grep("^CT", ignore.case = T, genes, value = T)
  return(setdiff(genes, c(mt.genes, many.genes, cox.genes, rp.genes, hb.genes, si.genes, zgc.genes, loc.genes, cr.genes, bx.genes)))
}

gene.list <- top20$gene
gene.list <- filter.genes(gene.list)

library(viridis)
pdf("pericyte_comparisons_dotplot_5.pdf", width = 12, height = 28)
DotPlot(obj.peri, features = unique(gene.list), dot.scale = 10, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

markers.peri.trans <- markers.peri$gene[which(markers.peri$cluster == "pericyte-0" & markers.peri$pct.2 < 0.08)]
markers.peri.1 <- markers.peri$gene[which(markers.peri$cluster == "pericytes_1" & markers.peri$pct.2 < 0.08)]
markers.peri.2 <- markers.peri$gene[which(markers.peri$cluster == "pericytes_2" & markers.peri$pct.2 < 0.08)]
markers.peri.3 <- markers.peri$gene[which(markers.peri$cluster == "mFB" & markers.peri$pct.2 < 0.08)]

##Plot the full list of differentially expressed pericyte genes across the pericyte populations - Figure S4C
gene.list <- unlist(unique(list(markers.peri.trans, markers.peri.1, markers.peri.2, markers.peri.3)))
gene.list.ordered <- factor(gene.list, levels = unique(c(markers.peri.trans, markers.peri.1, markers.peri.2, markers.peri.3, "acta2", "tagln", "myl9a")))
pdf("pericyte_comparisons_dotplot_5.pdf", width = 12, height = 28)
DotPlot(obj.peri, features = unique(gene.list.ordered), dot.scale = 10, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()


##Make a small dotplot highlighting some of the key differentially expressed markers between the pericyte populations - Figure 4C
gene.list <- c("ndufa4l2a", "abcc9", "rasgef1ba", "pdgfra", "pdgfrb", "ctgfa", "bgna", "dcn", "nt5c1aa", "esama", "epas1a", "adma", "olfml3b", "lmx1bb", "nptna", "sytl4", "pgfb", "acta2", "tagln", "myl9a")
pdf("pericyte_comparisons_dotplot_small_v3.pdf", width = 8, height = 10)
DotPlot(obj.peri, features = gene.list, idents = c("pericyte-0", "pericyte_1", "pericytes_2", "mFB", "vSMCs"), dot.scale = 10, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

##Conclusion: Differential expression of these populations revealed that these cells in cluster 21 differentially expressed ndufa4l2a from the rest of the SMCs in the same branch, but expressed acta2, tagln and all other markers that are expressed in vSMCs. Now as these cells doesn't have as many pericyte markers that are expressed here, maybe these are a SMC-->pericyte transient state (but still very much of a vSMC fate)

## III. Immature pericytes and vSMC markers: So far, immature and mature pericytes cannot be distinguished due to the lack of markers. When researchers have observed pericytes, it has not been clear where they came from and whether they are already mature

FeaturePlot(obj.smc, c("sytl4", "epas1a", "ebf1b", "pgfb")) ## expressed in the pericyte transient state overlapping ndufa4l2a and abcc9 expression and turns off in other pericyte clusters

##Similarly what maybe some markers of immature vSMCs?
FeaturePlot(obj.smc, c("foxc1a", "mef2cb", "slc7a2", "adrb3a")) ## slc7a2 and mef2cb seems to be good markers of this vSMC-transient state. 

##Check what stages these pericyte and SMc transient states belong to
cells.smc.transient <- WhichCells(obj.smc, idents = c("12"))
cells.peri.transient <- WhichCells(obj.smc, idents = c("9"))
table(obj.smc@meta.data[cells.peri.transient, "stage.nice"])
##Both can be found at 5 dpf. 


##IV. How are vascular and visceral SMCs transcriptoinally different from each other? - Figure S5B
##Plot a differential gene expression dotplot between all vascular and visceral SMC populations including  the unknown ones
obj.vi <- subset(obj.smc, idents = c("8", "10", "13", "16", "22", "23", "2", "11", "12", "15", "21"))
genes.vi <- FindAllMarkers(obj.vi, assay = "RNA", only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
my_levels <- c("8", "10", "13", "16", "22", "23", "2", "11", "12", "15", "21")
Idents(obj.vi) <- factor(x = Idents(obj.vi), levels = my_levels)
library(dplyr)
top5 <- genes.vi %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

##Specify gene list
gene.list <- c("kcnk18", "fsta", "actc1a", "acta1b", "acta1a", "cald1b", "CR753876.1", "il13ra2", "tesca", "fhl1a", "csrp1b", "nkx3.3", "prdx1", "rprmb", "hsd11b2", "mdka", "smtnb", "smtna", "zgc:195173", "vwde", "foxl1", "col4a5", "igfbp7", "col6a1", "col6a2", "nrg1", "evx1", "hoxa13a", "hoxa13b", "bambia", "lrrc17", "corin", "rergla", "krt4", "wscd2", "cyt1l", "sfrp2", "wfdc2", "pax1a", "cxcl12a", "vim", "crip1", "ldha", "pgk1", "aldocb", "eno3", "fbp2", "prss35", "rbp4", "acta2", "tagln", "si:dkey-164f24.2", "elnb", "loxa", "si:dkey-57k2.6", "fbln5", "pgam1a", "eno1a", "aldob", "ak1", "pfkpa", "cnn1a")
pdf("allSMC_vascular_visceral_vlnplot_v1.pdf", width = 10, height = 16)
VlnPlot(obj.vi, unique(gene.list), idents = factor(levels(Idents(obj.vi)), levels = c("8", "10", "13", "16", "22", "23", "2", "11", "12", "15", "21")), stack = TRUE, flip=TRUE, combine = FALSE) + NoLegend()
dev.off()


pdf("allSMC_vascular_visceral_dotplot_v1.pdf", width = 18, height = 28)
DotPlot(obj.vi, features = gene.list, dot.scale = 10, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()


##V. Plot gene expression for specific genes on the UMAP projections - Figure 4D, Figure 5B, and Figure S4D,E and S5D
##Ok, now plot the featureplots with the new URD color scheme - Feature Plots shown in Figures 3D and 4B were plotted using this same code
# Define the genes to plot
genes.plot <- c("ndufa4l2a")

# Create a list of plots
fplot_peri <- lapply(genes.plot, function(gene) {
  FeaturePlot(object = obj.smc, features = gene, pt.size = 2.2)
})

# Set the gradient color scale
scale <- scale_color_gradientn(colors = colors_feature_plot, na.value = "#CECECE") 

# Customize the plot theme
theme <- theme(panel.background = element_rect(fill = "white"),
               panel.grid = element_blank(),
               axis.line = element_line(colour = "black"),
               axis.text = element_text(colour = "black"),
               axis.title = element_text(colour = "black"),
               legend.text = element_text(colour = "black"),
               legend.title = element_text(colour = "black"))

# Modify the plots with the custom theme and color scale
fplots <- lapply(fplot_peri, function(plot) {
  plot + scale + theme 
})

dpi <- 300
png(paste0("FeaturePlot_", genes.plot, "_fig.png"), width = 2*dpi, height = 2*dpi)
fplots
dev.off()



##Compare gene expression between pericytes and vsSMCs
##gene.list <- filter.genes(gene.list)
##gene.list <- c("ggt5a", "pdgfab", "pla1a", "myocd", "pdlim3a", "pdlim3b", "rprmb", "nptna", "vcla", "nkx3.3", "nkx2.3", "smtna", "smtnb", "desmb", "pdlim1", "foxc1b", "vim", "fbln5", "wfdc2", "cd248a", "loxl5b", "loxa", "cnn1b", "tagln", "acta2", "ndufa4l2a", "abcc9", "pdgfrb", "rasl12", "kcnj8", "kcne4", "notch2", "notch3", "bgnb", "tm4sf18", "pros1", "c1qtnf5", "rasgef1ba")
##library(viridis)
##pdf("pericyte_vs_SMC_dotplot.pdf", width = 10, height = 16)
##DotPlot(obj.peri, features = gene.list, dot.scale = 10) + coord_flip() +
 ## geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
##  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
##  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
## dev.off()


##VI. Comparing gene expression between the intestinal longitudinal and circular muscles - Figure 5A
obj.vi <- subset(obj.smc, idents = c("8", "10", "13"))
genes.vi <- FindAllMarkers(obj.vi, assay = "RNA", only.pos = T, min.pct = 0.75, logfc.threshold = 0.5)

library(dplyr)
top5 <- genes.vi %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
gene.list <- unlist(top5$gene)

##Create curated gene list
gene.list <- c("foxf2a", "gucy1a1", "kcnk18", "actc1a", "fsta", "smpx", "npnt", "stom", "il13ra2", "tesca", "myl9b", "acta2", "fhl1a", "krt92", "itga4", "postnb", "rgs2", "fhl3b", "pdgfra", "fgfr2", "bambia", "tmem88b", "cxcl12a", "inka1a", "smtna", "smtnb", "desmb")
pdf("visceralSMC_dotplot_v2.pdf", width = 10, height = 16)
DotPlot(obj.vi, features = gene.list, dot.scale = 10, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

FeaturePlot(obj.smc, features = c("smpx", "pln", "slitrk6", "slc1a9"))


#### VII. Analyzing vascular SMC markers 
##Color these SMCs with separate colors
##Find markers of these individual SMC populations
##vSMC-aorta markers
FeaturePlot(obj.smc, features = c("sfrp2", "pax1a", "si:dkeyp-106c3.1", "htra3b"))

genes.smc.1 <- markers.smc$gene[markers.smc$cluster == "2"]
genes.smc.2 <- markers.smc$gene[markers.smc$cluster == "12"]
genes.smc.3 <- markers.smc$gene[markers.smc$cluster == "15"]
genes.smc.4 <- markers.smc$gene[markers.smc$cluster == "11"]
genes.smc.5 <- markers.smc$gene[markers.smc$cluster == "21"]
genes.smc.6 <- markers.smc$gene[markers.smc$cluster == "vSMC-artery"]

##Make a slot in the main object to highlight the vascular SMC clusters
obj.smc@meta.data$smc_types <- NA
obj.smc@meta.data[WhichCells(obj.smc, idents = c("12")), "smc_types"] <- "vSMC-transient"
obj.smc@meta.data[WhichCells(obj.smc, idents = c("2")), "smc_types"] <- "vSMC-aorta"
obj.smc@meta.data[WhichCells(obj.smc, idents = c("15")), "smc_types"] <- "vSMC-OFT"
obj.smc@meta.data[WhichCells(obj.smc, idents = c("11")), "smc_types"] <- "vSMC-artery_1"
obj.smc@meta.data[WhichCells(obj.smc, idents = c("21")), "smc_types"] <- "vSMC-artery_2"

DimPlot(obj.smc, group.by = "smc_types") + NoLegend() + NoAxes()

obj.smc@meta.data$smc_types <- NA
obj.smc@meta.data[WhichCells(obj.smc, idents = c("13")), "smc_types"] <- "viSMC_progenitors"
obj.smc@meta.data[WhichCells(obj.smc, idents = c("10")), "smc_types"] <- "int-SMCs_longitudinal"
obj.smc@meta.data[WhichCells(obj.smc, idents = c("8")), "smc_types"] <- "int-SMCs_circular"

colors <- setNames(c("#FF0000", "#023858", "#800000"), c("viSMC_progenitors", "int-SMCs_longitudinal", "int-SMCs_circular"))
DimPlot(obj.smc, group.by = "smc_types", cols = colors) + NoLegend() + NoAxes()

##Get the non-neural crest trajectory of mFB and vSMCs
DimPlot(obj.smc, label = T)
obj.smc@meta.data$non_ncc <- NA
obj.smc@meta.data[WhichCells(obj.smc, idents = c("11")), "non_ncc"] <- "vSMC_artery_1"
obj.smc@meta.data[WhichCells(obj.smc, idents = c("21")), "non_ncc"] <- "vSMC_artery_2"
obj.smc@meta.data[WhichCells(obj.smc, idents = c("3")), "non_ncc"] <- "mFB"

p <- DimPlot(obj.smc, group.by = "non_ncc")
LabelClusters(p, id = "ident", color = unique(ggplot_build(p)$data[[1]]$colour))

##Check co-expression between cluster 11 and cluster 3
obj <- subset(obj.smc, idents = c("11", "3"))
FeatureScatter(obj, feature1 = "acta2", feature2 = "fgf10a", pt.size = 4.0)

obj.smc <- subset(obj1, idents = c("11", "3"))
FeatureScatter(obj.smc, feature1 = "acta2", feature2 = "foxf2a", pt.size = 4.0)
FeatureScatter(obj.smc, feature1 = "ndufa4l2a", feature2 = "foxf2a", pt.size = 4.0)

##Find some more markers for intestinal SMCs that maybe spatially regionalized
##Look for markers overlapping il13ra2 but not overlapping tesca (tesca localized to mid-posterior intestine)
DimPlot(obj.smc, label = T)
FeaturePlot(obj.smc, features = c("rapgef5b", "cilp", "sptb", "npb", "htr2b", "adamtsl2", "kcnk4a", "pnck", "trpc1", "mpped1", "jph1b", "tesca"))




##After the number of ndufa4l2a and epas1a expressing cells were counted from the microscopy data, the following code was used to plot the relevant graphs shown in Figure 4H and 4L. 

##Load the pericyte proportion table
dt <- readxl::read_xlsx("~/Box/Farrell Lab/pericytes_SMC_HCR/2023-03-13_pericyte_proportions.xlsx")
colnames(dt) <- c("animal", "forebrain", "hindbrain", "pharyngeal_arches", "eye")
# Melt the data frame to convert the data from wide format to long format
library(reshape2)
df_long <- melt(dt, id.vars = "animal", variable.name = "region", value.name = "proportion")
# Reshape data from wide to long format
df_long <- tidyr::gather(dt, key = "region", value = "proportion", -animal)


# Calculate mean and standard error of mean
summary_df <- aggregate(proportion ~ region, data = df_long, FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))

# Create a line graph using ggplot2
library(viridisLite)
# calculate means and standard errors
# Jitter plot with mean and error bars - Figure 4I
ggplot(df_long, aes(x = region, y = proportion, color = animal)) +
  geom_jitter(width = 0.2, height = 0, size = 2) +
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar", width = 0.2, color = "black") +
  scale_color_viridis_d() +
  labs(x = "Region", y = "Proportion of cells expressing gene") +
  theme_classic()

##Plot the number of animals that had either 1, 2, 3 or 4 epas1a+ pericytes - Figure 4M
df <- readxl::read_xlsx("~/Box/Farrell Lab/pericytes_SMC_HCR/2023-03-13_pericyte_2_numbers_table.xlsx")
barplot(df$`number of animals`, names.arg = df$type, horiz = T)




