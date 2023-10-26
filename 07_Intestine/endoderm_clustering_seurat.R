#This code is directed towards using the endoderm tissue-specific atlas to understand differential gene expression across specific cell types. This code is used to generate the following figure panels: Figure 5A-C, E-F" and Figure S10A, B

##Load libraries
library(Seurat)
library(URD)

##Load endoderm seurat object
obj1 <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/endoderm/obj_seurat/endoderm_seurat.rds")

##Load endoderm markers object
markers <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/endoderm/obj_seurat/endoderm_markers.rds")

save.path <- "~/Box/zfext/results/06b-URD/obj_subsets/endoderm_superset_new/"
plot.path <- "~/Box Sync/zfext/results/06b-URD/plot_subsets/"

############################ PLOTTING THE ENDODERM SINGLE-CELL ATLAS ####################################

##Use the same stage colors that have been used in Figure 1A, B - keep that consistent throughout
stage.colors.new <- c(
  colorRampPalette(c("#f4be1d", "#efb918", "#bd8700", "#b37d00", "#a97300", "#9f6900", "#955f00", "#8b5500"))(8), # 3, 4, 5, 6, 7, 8, 9, 10
  rep(RColorBrewer::brewer.pal(6, "Oranges"), each = 2), # 24, 26 / 28, 30 / 32, 34 / 36, 38 / 40, 42 / 44, 46
  colorRampPalette(c("#FFC0CB", "#FFB6C1", "#FF69B4", "#DA70D6", "#BA55D3", "#C71585", "#FF1493"))(7), # 11, 12, 1`4, 16, 18, 21`
  rep(RColorBrewer::brewer.pal(6, "Greens"), each = 2), # 48, 50 / 52, 54 / 56, 58 / 60, 62 / 64, 66 / 68, 70
  rep(RColorBrewer::brewer.pal(6, "Blues"), each = 2),  # 72, 74 / 76, 78 / 80, 82 / 84, 86 / 88, 90 / 92, 94
  rep(RColorBrewer::brewer.pal(7, "Purples"), each = 2) # 96, 98 / 100, 102 / 104, 106 / 108, 110 / 112, 114 / 116, 118 / 120
)

stage.colors.group <- c(
  colorRampPalette(c("#f4be1d", "#efb918", "#bd8700", "#b37d00", "#a97300", "#9f6900", "#955f00", "#8b5500"))(8), # 3, 4, 5, 6, 7, 8, 9, 10
  rep(RColorBrewer::brewer.pal(6, "Oranges"), each = 2), # 24, 26 / 28, 30 / 32, 34 / 36, 38 / 40, 42 / 44, 46
  colorRampPalette(c("#FFC0CB", "#FFB6C1", "#FF69B4", "#DA70D6", "#BA55D3", "#C71585", "#FF1493"))(7), # 11, 12, 1`4, 16, 18, 21`
  rep(RColorBrewer::brewer.pal(6, "Greens"), each = 2), # 48, 50 / 52, 54 / 56, 58 / 60, 62 / 64, 66 / 68, 70
  rep(RColorBrewer::brewer.pal(6, "Blues"), each = 2),  # 72, 74 / 76, 78 / 80, 82 / 84, 86 / 88, 90 / 92, 94
  rep(RColorBrewer::brewer.pal(7, "Purples"), each = 2) # 96, 98 / 100, 102 / 104, 106 / 108, 110 / 112, 114 / 116, 118 / 120
)


##Plot  the endoderm UMAP plot colored by stage (hpf) - Figure S8A
dpi <- 300
png(file = "endoderm_stage.nice.png", width = 5 * dpi, height = 5 * dpi)
DimPlot(obj1, group.by = "stage.nice", cols = stage.colors.new, pt.size = 1.2)
dev.off()

cols.stage <- setNames(c("#B0C4DE", "#708090", "#0000FF", "#2b8cbe", "#17A589", "#00FF00", "#f4be1d",
                                  "#8b5500", "#FF7F00", "#F08080", "#FF1493", "#FF00FF", "#8B008B"), c(unique(obj1@meta.data$stage.group)))
dpi <- 300
png(file = "endoderm_stage.group_v2.png", width = 5 * dpi, height = 5 * dpi)
DimPlot(obj1, group.by = "stage.group", cols = cols.stage, pt.size = 1.2) + NoAxes() + NoLegend()
dev.off()


##Plot the endoderm UMAP plot colored by clusters
dpi <- 300
png(file = "endoderm_UMAP_plot.png", width = 5 * dpi, height = 5 * dpi)
DimPlot(obj1, label = T, pt.size = 1) + theme(legend.position="none") + NoAxes()
dev.off()

##Define a broader level of clustering for downstream purposes.
##broad_clust_ids <- c("hepatocytes", "hepatocytes", "hepatocytes", "hepatocytes", "hepatocytes", "hepatocytes", "hepatocytes", "cholangiocytes", "progenitor_intestine-1", "enterocyte", "enterocyte", "enterocyte", "enterocyte", "progenitor_intestine-2", "secretory", "enterocyte", "secretory", "secretory",  "endocrine pancreas", "endocrine pancreas", "endocrine pancreas", "exocrine pancreas", "exocrine pancreas", "progenitors", "progenitors", "progenitors", "progenitors", "progenitors", "pharynx", "esophagus", "cloaca", "pneumatic-duct/SB", "thyroid", "doublets")
                     
#names(broad_clust_ids) <- levels(obj1)
#obj2 <- RenameIdents(obj1, broad_clust_ids)
#DimPlot(obj2, label = T)

#endo_colors <- c("#c89721", "#A0522D", "#48D1CC", "#9288eb", "#5b4e9f", "#02818a", "#f392be", "#fb75b1", "#a84b99", "#FF7F50", "#32CD32")
#endo.organ <- c("liver", "intestine", "exocrine pancreas", "progenitors", "pharynx", "pneumatic-duct/SB", "esophagus", "cloaca", "endocrine pancreas", "thyroid", "doublets")
#endo.tissue.cols <- setNames(endo_colors, endo.organ)
#DimPlot(obj2, label = T, cols = endo.tissue.cols)

##Plot UMAP plot colored by cluster for Figure 6A.
##Define the colors for each tissue group within the endoderm atlas
cluster_ids <- c("HP-4", "HP-3", "HP-2", "HP-1", "hepatoblasts", "late-HBs", "hepatointestinal", "cholangiocytes", "int_progenitors-1", "enterocyte-1", "enterocyte-2", "enterocyte-3", "posterior_LREs", "int_progenitors-2", "goblet-cells", "best4+_enterocytes", "EECs", "tuft-like", 
"ε-cells", "α/β/δ-cells", "PP-cells", "exo_panc_progenitors", "centroacinar", "progenitor_1", "progenitor_2", "progenitor_3", "progenitor_4",  "hepatopancreatic progenitors", "pharynx", "esophagus", "cloaca", "pneumatic-duct", "thyroid", "doublets")

liver_colors <- c("#D2B48C", "#DEB887", "#F4A460", "#CD853F", "#D2691E", "#A0522D")
liver_names <- c("HP-4", "HP-3", "HP-2", "HP-1", "hepatoblasts", "ad-HBs")
liver.cols <- setNames(liver_colors, liver_names)

intestine_colors <- c("#ADD8E6", "#B0E0E6", "#87CEEB", "#A6BDDB", "#74A9CF", "#3690C0", "#3182bd", "#023858", "#08519c", "#000080")
intestine_names <- c("ant-int_progenitors", "enterocyte-1", "enterocyte-2", "enterocyte-3", "posterior_LREs", "post_int_progenitors", "goblet-cells", "best4+_otop2+", "EECs", "tuft-like")
intestine.cols <- setNames(intestine_colors, intestine_names)

endo_panc_colors <- c("#FA8072", "#FF6347", "#FF0000")
endo_panc_names <- c("δ/ε-cells","α/β-cells", "PP-cells")
endo.panc.cols <- setNames(endo_panc_colors, endo_panc_names)

exo_panc_colors <- c("#DDA0DD", "#BA55D3")
exo_names <- c("exo_panc_progenitors", "centroacinar")
exo.panc.cols <- setNames(exo_panc_colors, exo_names)

colors.rest <- c("#00CED1", "#FF1493", "#008000", "#8FBC8F", "#800000", "#00FA9A", "#FF00FF", "#008080", "#DB7093")
names.rest <- c("progenitor_1", "progenitor_2", "progenitor_3", "progenitor_4", "pharynx", "esophagus", "cloaca", "hepatopancreatic progenitors", "sftpba+_cells")
rest.cols <- setNames(colors.rest, names.rest)

##make a list of all colors
all.colors <- c(liver.cols, intestine.cols, endo.panc.cols, exo.panc.cols, rest.cols)

##Adjust levels such that cluster numbers are associated with cluster names
names(cluster_ids) <- levels(obj1)
obj.endo <- RenameIdents(obj1, cluster_ids)
DimPlot(obj.endo, label = T)

##Plot the UMAP
DimPlot(obj.endo, label = F, cols = all.colors)


##Plot the differentially expressed transcription factors expressed between broader cell types within the endodermal atlas - Figure S8B
##Plotting the broader clusters in UMAP
endo.organ <- c("liver", "intestine", "exocrine pancreas", "endocrine pancreas", "progenitor_1", "progenitor_2", "progenitor_3", "progenitor_4", "pharynx", "esophagus", "cloaca", "hepatopancreatic progenitors", "sftpba+_cells")
endo_colors <- c("#DAA520", "#00BFFF", "#9932CC", "#FF6347", "#00CED1", "#FF1493", "#008000", "#8FBC8F", "#800000", "#00FA9A", "#FF00FF", "#008080", "#DB7093")
endo.tissue.cols <- setNames(endo_colors, endo.organ)
DimPlot(obj2, label = F, cols = endo.tissue.cols)

dpi <- 300
png("endoderm_UMAP_broad_clusters.png", width = 5 * dpi, height = 5 * dpi)
DimPlot(obj2, label = F, cols = endo.tissue.cols, pt.size = 1.2) + theme(legend.position="none") + NoAxes()
dev.off()

stage.colors <-
dpi <- 300
png("endoderm_UMAP_stage.group.png", width = 5 * dpi, height = 5 * dpi)
DimPlot(obj1, group.by = "stage.group", label = F, pt.size = 1.2) + theme(legend.position="none") + NoAxes()
dev.off()

#Labeling the clusters using the same color as the clusters
p <- DimPlot(obj2, label = F, cols = endo.tissue.cols) + theme(legend.position="none") + NoAxes()
LabelClusters(p, id = "ident", color = unique(ggplot_build(p)$data[[1]]$colour))

##Plotting a dotplot between the major cluster groups in the endoderm
##For that first I need to do FindMarkers
markers.endo <- FindAllMarkers(obj.merge, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

##Only focus on TFs
tf.list <- read.delim(file="~/Box/zfext/02-Clustering/2021-03 Iterative Clustering/gene_info/tfs/2021-06-25_zebrafish_LTA_TFs.txt", header = T, sep = "")
markers.tfs.only <- markers.endo[which(rownames(markers.endo) %in% tf.list$Symbol), ]
#Take only top 5 of the TFs expressed in the different endoderm cell types
top5 <- markers.tfs.only %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
gene.list <- unlist(top5$gene)

#Plot DotPlot
p <- DotPlot(obj1, features = gene.list) + coord_flip()

##Let's replicate this dotplot using ComplexHeatMap
library(ComplexHeatmap)
library(tidyverse)
library(circlize)

df<- p$data
head(df)

### the matrix for the scaled expression 
exp_mat <-df %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

row.names(exp_mat) <- exp_mat$features.plot  
exp_mat <- exp_mat[,-1] %>% as.matrix()

head(exp_mat)

## the matrix for the percentage of cells express a gene

percent_mat <- df %>% 
  select(-avg.exp, -avg.exp.scaled) %>%  
  pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame() 

row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()
percent_mat <- percent_mat/100

head(percent_mat)

## the range is from 0 - 100
range(percent_mat)

## these two matrix have the same dimension
dim(exp_mat)
dim(percent_mat)

library(viridis)
library(Polychrome)

Polychrome::swatch(viridis(20))

## get an idea of the ranges of the matrix
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99))

## any value that is greater than 2 will be mapped to yellow
col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(20)[c(1,10, 20)])


cell_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.circle(x=x,y=y,r= percent_mat[i, j]/100 * min(unit.c(w, h)),
              gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}


##Write a function to plot the dotplot
getDotPlot <- function(cl.expr.mat, cl.prop.mat, gene.reorder=TRUE, specified.order=NULL, cl.reorder.by.hc=TRUE, genes=NULL, colors=c("#d9d9d9", "#252525"), point.size.factor=3.5, plot.name=NULL, plot.height=20, plot.width=30, plot.margin=c(8,12,5,8), plot.cex=1.5, max.diag=TRUE, col.space.factor=1, row.space.factor=1){
  if(is.null(genes)){
    genes <- rownames(cl.expr.mat)
  }
  expr <- cl.expr.mat[genes,]
  sd.vec <- apply(expr, 1, sd)
  expr <- expr[which(sd.vec>0),]
  if(gene.reorder){
    print("Reorder genes based on hierarchical clustering")
    hc.gene <- hclust(as.dist(1-cor(t(expr))))
    gene.order = hc.gene$order
  }else{
    print("No reordering for genes")
    gene.order <- rownames(expr)
  }
  expr <- expr[gene.order,] 
  
  if(!is.null(specified.order)){
    print("Reorder clusters according to user provided order")
    expr <- expr[,specified.order]
  }else{
    if(cl.reorder.by.hc){
      print("Reorder clusters based on hierarchical clustering of provided expression matrix")
      hc <- hclust(as.dist(1-cor(expr)))
      expr <- expr[,hc$order]
    }else{
      print("Do not reorder clusters")
      expr <- expr
    }
  }
  
  if(max.diag){
    max.cl.idx <- colnames(expr)[apply(expr, 1, which.max)]
    idx <- unlist(lapply(colnames(expr), function(i){
      which(max.cl.idx==i)
    }))
    expr <- expr[idx,]
  }
  scale.cl.expr <- t(scale(t(expr)))
  scale.cl.expr <- rbind(scale.cl.expr, rep(c(-1,1), c(ncol(scale.cl.expr)-4,4)))
  colorPal <- grDevices::colorRampPalette(colors)
  color.break.num <- 30
  cellColor <- as.vector(t(apply(scale.cl.expr, 1, function(values){
    adjustcolor(colorPal(color.break.num), alpha=.8)[as.numeric(cut(values, breaks=color.break.num, right=F, include.lowest=T))]
  })))
  cl.num <- ncol(expr)
  prop <- cl.prop.mat[rownames(expr),colnames(expr)]
  new.vec <- rep(0, ncol(prop))
  new.vec[ncol(prop)-0:3] <- seq(from=1,to=0.25,length=4)
  prop <- rbind(prop, new.vec)
  point.size <- as.vector(prop)*point.size.factor
  
  g1 <- c(rownames(expr), "")
  if(is.null(plot.name)){
    plot.name <- "DotPlot_selected_genes.pdf"
  }
  pdf(plot.name, height=plot.height, width=plot.width)
  par(mar=plot.margin, xpd=TRUE)
  plot(c(1, length(g1)), c(1, cl.num), type="n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
  points(rep(seq(length(g1))*col.space.factor, cl.num), 1+rep(seq(cl.num)*row.space.factor, each=length(g1)), pch=16, col=cellColor, cex=point.size)
  mtext(g1, side=1, at=seq(length(g1))*col.space.factor, las=2, cex=plot.cex)
  mtext(sub("Cluster_", "", colnames(expr)), side=2, at=1+seq(cl.num)*row.space.factor, las=1, cex=plot.cex)
  #points(rep((length(g1)+2)*col.space.factor, 4), c(0.5*cl.num-2.9, 0.5*cl.num-1.75, 0.5*cl.num-0.5, 0.5*cl.num+0.75), pch=16, cex=c(0.25,0.5,0.75,1)*point.size.factor,xpd=TRUE)
  text(rep((length(g1)+2)*col.space.factor, 4)+0.5, cl.num-3:0, labels=c(0.25,0.5,0.75,1), cex=plot.cex, pos=4)
  dev.off()
}

##Define the order of the cell types in which you want them to be plotted
cell.type.order <- c("progenitors", "hepatocytes", "cholangiocytes", "enterocyte", "secretory", "progenitor_intestine-1", "progenitor_intestine-2", "pharynx", "esophagus", "cloaca", "pneumatic-duct/SB", "thyroid", "doublets")

endo_colors <- c("#c89721", "#A0522D", "#48D1CC", "#9288eb", "#5b4e9f", "#02818a", "#f392be", "#fb75b1", "#a84b99", "#FF7F50", "#32CD32", "#48D1CC")
endo.organ <- c("liver", "intestine", "exocrine pancreas", "progenitor_1", "progenitor_4", "pharynx", "sftpba+_cells", "esophagus", "cloaca", "progenitor_2", "endocrine pancreas", "progenitor_3")
endo.tissue.cols <- setNames(endo_colors, endo.organ)

##Plot the dotplot showing differentially expressed TFs - Figure S10B
getDotPlot(cl.expr.mat= exp_mat, cl.prop.mat=percent_mat, gene.reorder=FALSE, specified.order=NULL, cl.reorder.by.hc=FALSE, genes=top5$gene, colors=c("#d9d9d9", "#252525"), point.size.factor=5.0, plot.name="Plot_dot_plot_cell_endoderm_RNA_LOG_expr_all_2.pdf", plot.height=100, plot.width=80, max.diag=FALSE, col.space.factor=0.7, row.space.factor=0.12) + coord_flip()


################# CHARACTERIZING THE PNEUMATIC DUCT ##############################################
#Let's create an object with intestine, exo_panc, endo_panc, liver, esophagus, pharynx, cloaca, and sftpba-cells and compare gene expression in a dotplot and featurePlots

##Plot gene expression on the UMAP projection for pneumatic duct markers
##Plot feature plots for Figure 5 using the new URD color scheme

# Define the genes to plot
genes.plot <- c("mnx1")

# Create a list of plots - Figure 5B
fplot_peri <- lapply(genes.plot, function(gene) {
  FeaturePlot(object = obj1, features = gene, pt.size = 2.2, order = F)
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

##Now plot the dotplot
##Define clustering
broad_clust_ids <- c("liver", "liver", "liver", "liver", "liver", "liver", "liver", "liver", "intestine", "intestine", "intestine", "intestine",
                     "intestine", "intestine", "intestine", "intestine", "intestine", "intestine", "endocrine pancreas", "endocrine pancreas",
                     "endocrine pancreas", "exocrine pancreas", "exocrine pancreas", "progenitors", "progenitors", "progenitors", "progenitors",
                     "progenitors", "pharynx", "esophagus", "cloaca", "pneumatic-duct", "thyroid", "doublets")

names(broad_clust_ids) <- levels(obj1)
obj2 <- RenameIdents(obj1, broad_clust_ids)
DimPlot(obj2, label = T)

endo_colors <- c("#c89721", "#A0522D", "#48D1CC", "#9288eb", "#5b4e9f", "#02818a", "#f392be", "#fb75b1", "#a84b99", "#FF7F50", "#32CD32", "#48D1CC")
endo.organ <- c("liver", "intestine", "exocrine pancreas", "progenitor_1", "progenitor_4", "pharynx", "sftpba+_cells", "esophagus", "cloaca", "progenitor_2", "endocrine pancreas", "progenitor_3")
endo.tissue.cols <- setNames(endo_colors, endo.organ)
DimPlot(obj2, label = T, cols = endo_colors)

##Use markers for pneumatic-duct
genes.to.use <- rownames(markers$wilcox$`32`)
genes.to.use <- filter.genes(genes.to.use)
genes.top <- c(head(genes.to.use, 20), "abca3b", "mnx1", "sox2", "anxa5b", "arnt2", "lrata", "slit1b", "sim1b", "sim2", "sim2.1", "sftpba", "ihha", "shha", "aplp2", "cxcl20", "tfa", "ccl20b", "sgms1", "lamb1b")

#Plot DotPlot - Figure 5C
pdf("sftpba_cells_dotplot.pdf", width = 14, height = 20)
DotPlot(obj2, features = genes.top, dot.scale = 16, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

pdf("endoderm_broad_UMAP_clusters.pdf", width = 6, height = 6)
DimPlot(obj2, label = F, cols = endo_colors) + NoLegend() + NoAxes()
dev.off()


########################### CHARACTERIZING BEST4+ ENTEROCYTES ###########################################

##Create a dotplot of all genes differentially expressed in best4+ cells (cluster 16) against other intestinal cells
obj.intestine <- subset(obj.endo, idents = c("int_progenitors-1", "enterocyte-1", "enterocyte-2", "enterocyte-3", "posterior_LREs", "int_progenitors-2", "goblet-cells", "best4+_enterocytes", "EECs", "tuft-like"))

##Calculate markers between between intestinal cells
markers.int <- FindAllMarkers(obj.intestine, assay = "RNA", logfc.threshold = 0.25, only.pos = T, min.pct = 0.25)
##Get markers enriched in best4+ cells
genes.to.plot <- rownames(markers$wilcox$`16`)

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

genes.to.plot <- filter.genes(genes.to.plot)
genes.bo4 <- genes.to.plot[1:50]

#Plot DotPlot - Figure 5E
pdf("best4_otop2_genes_dotplot.pdf", width = 12, height = 16)
DotPlot(obj.intestine, features = genes.bo4, dot.scale = 10, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

## Now plot select genes for best4+ enterocytes on the UMAP projections
##Use the new URD color scheme
# Define the genes to plot
genes.plot <- c("best4", "otop2", "cdx1b")

# Create a list of plots
fplot_peri <- lapply(genes.plot, function(gene) {
  FeaturePlot(object = obj1, features = gene, pt.size = 2.2, order = F)
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


##Plotting dotplots and gene expression on UMAP for LREs
##Do the same dotplot for the posterior lysozome-rich enterocytes - not shown in the paper. 
genes.to.plot <- rownames(markers$wilcox$`15`)

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

genes.to.plot <- filter.genes(genes.to.plot)
genes.lre <- genes.to.plot[1:50]

#Plot DotPlot
pdf("posterior_LREs_genes_dotplot.pdf", width = 20, height = 24)
DotPlot(obj.intestine, features = genes.lre) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

saveRDS(obj1, "~/Box/zfext/annotations_celltype_curated/endoderm/obj_seurat/endoderm_seurat.rds")

##Plot feature plots
gridExtra::grid.arrange(grobs = lapply(c("cdx1b", "pbx3a", "mafa", "atoh1b"), plotTree,
                                       object = obj.build, label.x = F, plot.cells = T), ncol = 2)
gridExtra::grid.arrange(grobs = lapply(c("foxa3", "ptf1a", "cdx1b", "arxa"), plotTree,
                                       object = obj.build, label.x = F, plot.cells = T), ncol = 2)




############## COMPARE BEST4+ CELLS TO ABSORPTIVE AND SECRETORY INTESTINAL CELLS #####################
## As reviewer 1 asked, the identity of best4+ cells remain debatable - Wen et al. 2020 calls best4+ cells as ionocytes, Willms et al., 2022 calls best4+ cells absorptive, the human papers also call them enterocytes. But how do they compare to absorptive and secretory cell types? 

##To do this, I need to get a list of genes that define ionocytes, absorptive cells, and secretory cells. So first get a unique signature of ionocytes

##Load mama
mama <- readRDS("~/Box/zfext/02-Clustering/mama+ds_integrated/merged_mama_ds_slim_seurat_25Jan2022.rds")

# Load cell annotations to introduce for downsample purposes
cell.annot.path <- "~/Box/zfext/annotations_celltype_curated_newMama/celltype_annotations_df_full.tsv"
message(Sys.time(), ": Loading ", cell.annot.path)
cell.annot <- read.table(cell.annot.path, sep = " ", header = T, stringsAsFactors = F)

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


##Get all absorptive cell types
cells.absorptive <- WhichCells(obj1, idents = c("10", "11", "12"))
cells.secretory <- WhichCells(obj1, idents = c("15", "17", "18"))
cells.best4 <- WhichCells(obj1, idents = c("16"))

##Get genes for absorptive intestinal cells that are unique to them - what genes are characteristic of absorptive cells?

##Assign absorptive cells as an identity and non-absorptive cells as a different identity for comparisons
obj.intestine@meta.data$cell.groups <- NA
obj.intestine@meta.data[cells.absorptive, "cell.groups"] <- "absorptive"
obj.intestine@meta.data[setdiff(rownames(obj.intestine@meta.data), cells.absorptive), "cell.groups"] <- "non-absorptive"

##Set the cell.groups clustering as idents
Idents(obj.intestine) <- obj.intestine@meta.data$cell.groups

##Get differential markers for intestine absorptive cells
##Calculate markers for the 2 groups
m.abs <- FindAllMarkers(obj.intestine, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)


##Get the markers that are specific to the absorptive cluster
markers.absorptive <- m.abs[which(m.abs$cluster == "absorptive"), "gene"]
##Save the list of absorptive markers
write(markers.absorptive, "~/Box/zfext/annotations_celltype_curated_newMama/endoderm/best4_markers/markers_absorptive_intestine_vs_intestine.txt")
markers.absorptive <- scan("~/Box/zfext/annotations_celltype_curated_newMama/endoderm/best4_markers/markers_absorptive_intestine_vs_intestine.txt", what = "character")

##Kidney cells (proximal tubules) are also absorptive. Maybe look at the profile of kidney cells that are also absorptive
#cells.kidney <- WhichCells(mama, idents = c("pron.1", "pron.14", "pron.5", "pron.7"))
#mat.kidney <- mama@assays$RNA@data[, cells.kidney]
#genes.absorptive.kidney <- unique(rownames(mat.kidney)[which(mat.kidney > 0.25)])

##Take the intersection of intestine absorptive genes and kidney absorptive genes
#genes.absorptive.specific <- intersect(markers.absorptive, genes.absorptive.kidney)

##Do the same for secretory cells
obj.intestine@meta.data$cell.groups <- NA
obj.intestine@meta.data[cells.secretory, "cell.groups"] <- "secretory"
obj.intestine@meta.data[setdiff(rownames(obj.intestine@meta.data), cells.secretory), "cell.groups"] <- "non-secretory"

##Set the cell.groups clustering as idents
Idents(obj.intestine) <- obj.intestine@meta.data$cell.groups

##Get differential markers for intestine absorptive cells
##Calculate markers for the 2 groups
m.sec <- FindAllMarkers(obj.intestine, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

##Get the markers that are specific to the absorptive cluster
markers.secretory <- m.sec[which(m.sec$cluster == "secretory"), "gene"]
##Save the list of absorptive markers
write(markers.secretory, "~/Box/zfext/annotations_celltype_curated_newMama/endoderm/best4_markers/markers_secretory_intestine_vs_intestine.txt")
markers.secretory <- scan("~/Box/zfext/annotations_celltype_curated_newMama/endoderm/best4_markers/markers_secretory_intestine_vs_intestine.txt", what = "character")

#Get other secretory cells such as pancreas, mucous cells, otic/lateral line etc. 
##Maybe use agr2 to culminate all secretory cell types - I used Daniocell to pick clusters that express agr2 higher than 2 and use those to compare secretory cell types
Idents(mama) <- mama@meta.data$clust
cells.secretory.other <- WhichCells(mama, expression = agr2 >= 1)
cells.secretory.non.intestine <- setdiff(cells.secretory.other, cells.secretory.int)

##Find genes that are expressed in these cells over a fold change of 0.25
mat.secretory <- mama@assays$RNA@data[, cells.secretory.non.intestine]
genes.secretory.other <- unique(rownames(mat.secretory)[which(mat.secretory > 1)])

##Now take an intersection of genes that are expressed in intestinal secretory cells and those in other secretory cells
genes.secretory.specific <- intersect(genes.secretory.other, markers.secretory)
write(genes.secretory.specific, "~/Box/zfext/annotations_celltype_curated_newMama/endoderm/best4_markers/markers_secretory_specific.txt")



##Now get the ionocytes from many different sources - ionocytes from kidney and skin
Idents(mama) <- mama@meta.data$clust
cells.ionocytes <- WhichCells(mama, idents = c("iono.5", "iono.10", "iono.17", "iono.12", "iono.13", "iono.20"))

##To get the ionocyte specific genes, we need to compare the gene expression of the above clusters against the periderm where they are located

##Change identity assignment
Idents(mama) <- mama@meta.data$tissue.annot
cells.periderm <- WhichCells(mama, idents = c("periderm"))
cells.epidermis <- WhichCells(mama, idents = c("epidermis"))

##Load sample object
sample <- "ionocytes"
obj.iono <- readRDS(paste0(file = "~/Box/zfext/annotations_celltype_curated_newMama/", sample, "/obj_seurat/", sample, "_seurat.rds"))

##Get all ionocyte genes
obj.iono@meta.data$cell.type <- NA
obj.iono@meta.data[cells.ionocytes, "cell.type"] <- "ionocytes"
obj.iono@meta.data[setdiff(WhichCells(obj.iono), cells.ionocytes), "cell.type"] <- "mucous_and_progenitors"

##Set Ident
Idents(obj.iono) <- obj.iono@meta.data$cell.type
                                 
##Get markers specific to ionocytes from the ionocyte/mucous-secreting dataset
iono.markers <- FindAllMarkers(obj.iono, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

##Genes specific to ionocytes when compared to periderm (as most ionocytes are present in the skin)
markers.iono.sub <- iono.markers[which(iono.markers$cluster == "ionocytes" & iono.markers$avg_log2FC >= 0.25), "gene"]

##Now as most ionocytes are in the periderm or epidermis, we need to substract the periderm/epidermis signature
##Find genes expressed in periderm
mat.periderm <- mama@assays$RNA@data[, cells.periderm]
genes.periderm <- unique(rownames(mat.periderm)[which(mat.periderm > 0.25)])

##Find genes expressed in epidermis
mat.epidermis <- mama@assays$RNA@data[, cells.epidermis]
genes.epidermis <- unique(rownames(mat.epidermis)[which(mat.epidermis > 0.25)])

##Take union of genes expressed in periderm and epidermis and substract that from the ionocyte genes
genes.skin.total <- unlist(unique(list(genes.periderm, genes.epidermis)))
##Get genes that are ionocyte specific and not skin related
genes.iono.specific <- setdiff(markers.iono.sub, genes.skin.total)
write(genes.iono.specific, "~/Box/zfext/annotations_celltype_curated_newMama/endoderm/best4_markers/markers_ionocyte_specific.txt")

##Genes expressed in absorptive and secretory cells
mat.abs <- obj1@assays$RNA@data[, cells.absorptive]
mat.sec <- obj1@assays$RNA@data[, cells.secretory]

#Genes expressed in absorptive and secretory cells in the intestine based on expression threshold > 0.25
genes.absorptive <- unique(rownames(mat.abs)[which(mat.abs >= 0.25)])
genes.secretory <- unique(rownames(mat.sec)[which(mat.sec >= 0.25)])
genes.common.abs.sec <- intersect(genes.absorptive, genes.secretory)

genes.absorptive.only <- setdiff(genes.absorptive, genes.common.abs.sec)
genes.secretory.only <- setdiff(genes.secretory, genes.common.abs.sec)

genes.absorptive.only <- unlist(unique(list(genes.absorptive.only, markers.absorptive.high)))
genes.secretory.only <- unlist(unique(list(genes.secretory.only, markers.secretory.high)))

##Ok, so at this point we have extracted a list of markers specific to absorptive enterocytes (within the intestine), secretory  cells (within intestine as well as cumulative), and genes specific to ionocytes after substraction of periderm and epidermis genes

##Now I need to compare the 3 sets and make sure that there are no genes that overlap. 
## 1. For ionocytes: genes are specific to skin ionocytes and excludes periderm and epidermis genes
## 2. For absorptive: genes specific to absorptive cells in the intestine as well as shared with other absorptive cell types (e.g., kidney)
## 3. For secretory: genes expressed exclusively in secretory cells in the intestine and also shared across other secretory cell types. 

##Ok, only keep markers in absorptive, secretory and ionocyte cell types that are expressed at a certain level
markers.absorptive.high <- m.abs[which(m.abs$cluster == "absorptive" & m.abs$avg_log2FC >= 0.25), "gene"]
markers.secretory.high <- m.sec[which(m.sec$cluster == "secretory" & m.sec$avg_log2FC >= 0.25), "gene"]
markers.iono.high <- genes.iono.specific

##Save the lists of genes for each category
write(markers.absorptive.high, "~/Box/zfext/annotations_celltype_curated_newMama/endoderm/best4_markers/markers_absorptive_highly_expressed.txt")
write(markers.secretory.high, "~/Box/zfext/annotations_celltype_curated_newMama/endoderm/best4_markers/markers_secretory_highly_expressed.txt")
write(markers.iono.high, "~/Box/zfext/annotations_celltype_curated_newMama/endoderm/best4_markers/markers_ionocyte_highly_expressed.txt")


##Absorptive only, secretory only, and ionocyte only genes
genes.abs.only <- setdiff(markers.absorptive.high, unlist(unique(list(markers.secretory.high, genes.iono.specific))))
genes.sec.only <- setdiff(markers.secretory.high, unlist(unique(list(markers.absorptive.high, genes.iono.specific))))
genes.iono.only <- setdiff(genes.iono.specific, unlist(unique(list(markers.absorptive.high, markers.secretory.high))))

##Save the lists of genes for each category
write(genes.abs.only, "~/Box/zfext/annotations_celltype_curated_newMama/endoderm/best4_markers/markers_intestine_absorptive_only.txt")
write(genes.sec.only, "~/Box/zfext/annotations_celltype_curated_newMama/endoderm/best4_markers/markers_intestine_secretory_only.txt")
write(genes.iono.only, "~/Box/zfext/annotations_celltype_curated_newMama/endoderm/best4_markers/markers_ionocyte_only.txt")


##Now get the genes expressed in best4 cells in the intestine
cells.best4 <- WhichCells(obj1, idents = c("16"))
mat.best4 <- obj1@assays$RNA@data[, cells.best4]
genes.best4 <- unique(rownames(mat.best4)[which(mat.best4 > 0.25)]) ##2855

##Get genes expressed differentially in best4 cells among other cells in the intestine
markers.best4 <- markers$wilcox$`16`

##Now we have generated all different marker combinations including genes 

##Out of these 2855 genes, how many overlap with absorptive and secretory cells? 
genes.best4.abs.common <- intersect(genes.best4, genes.abs.only)
genes.best4.sec.common <- intersect(genes.best4, genes.sec.only)
genes.best4.iono.common <- intersect(genes.best4, genes.iono.only)

##Now define a threshold at which best4 cells express these genes
##Subset the gene expression matrix for these genes 
mat.best4.abs.common <- mat.best4[genes.best4.abs.common, ]
mat.best4.sec.common <- mat.best4[genes.best4.sec.common, ]
mat.best4.iono.common <- mat.best4[genes.best4.iono.common, ]

##With the different ranges of gene expression, I set a gene expression threshold of 3.5 - i.e., only genes considered to be expressed more than a FC of 3.5 will be considered expressed in best4 cells

genes.best4.abs.final <- unique(rownames(mat.best4.abs.common)[which(mat.best4.abs.common >= 0.75)])
genes.best4.sec.final <- unique(rownames(mat.best4.sec.common)[which(mat.best4.sec.common >= 1)])
genes.best4.iono.final <- unique(rownames(mat.best4.iono.common)[which(mat.best4.iono.common >= 0.5)])

##Save the shared gene lists
write(genes.best4.abs.common, "~/Box/zfext/annotations_celltype_curated_newMama/endoderm/best4_markers/markers_absorptive_shared_best4.txt")
write(genes.best4.sec.common, "~/Box/zfext/annotations_celltype_curated_newMama/endoderm/best4_markers/markers_secretory_shared_best4.txt")
write(genes.best4.iono.common, "~/Box/zfext/annotations_celltype_curated_newMama/endoderm/best4_markers/markers_ionocyte_shared_best4.txt")

##PLOT THE GENES

colors_feature_plot <- defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)

# Create a list of feature plots for absorptive specific genes
fplot_abs <- lapply(genes.secretory.only, function(gene) {
  FeaturePlot(object = obj1, features = gene, pt.size = 2, order = F)
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
fplots <- lapply(fplot_abs, function(plot) {
  plot + scale + theme 
})

##Now combine the feature plots into one document
##combined.feature.plots <- do.call(gridExtra::grid.arrange, c(fplots, ncol = 5))

# Step 6: Combine the histograms into multiple pages with 50 plots per page
num_plots <- length(fplots)
num_pages <- ceiling(num_plots / 25)

##Save the above histograms in a single pdf

pdf("~/Box/zfext/annotations_celltype_curated_newMama/endoderm/best4_markers/secretory_only_genes_endoderm.pdf", height = 48, width = 36)

for (page in 1:num_pages) {
  start_plot <- (page - 1) * 25 + 1
  end_plot <- min(page * 25, num_plots)
  plot_subset <- fplots[start_plot:end_plot]
  combined_plots <- do.call(grid.arrange, c(plot_subset, ncol = 5))  # Adjust ncol as needed
  print(combined_plots)
}

dev.off()

##I guess we need to get rid of genes that are expressed universally across the whole dataset that are also shared between best4+ cells and absorptive and secretory cell types
##Subset mama data slot to only the best4+ vs absorptive shared genes 
mat.mama.abs <- mama@assays$RNA@data[is.na(genes.best4.abs.final), ]
mat.mama.sec <- mama@assays$RNA@data[is.na(genes.best4.sec.final), ]
mat.mama.iono <- mama@assays$RNA@data[is.na(genes.best4.iono.final), ]

obj.merge <- subset(mama, cells = unlist(unique(list(cells.absorptive, cells.secretory, cells.ionocytes, cells.best4))))

##Plot a dotplot for the genes differentially expressed in absorptive, secretory, and ionocyte cell types
obj.merge@meta.data$cell.groups <- NA
obj.merge@meta.data[cells.absorptive, "cell.groups"] <- "absorptive"
obj.merge@meta.data[cells.secretory, "cell.groups"] <- "secretory"
obj.merge@meta.data[cells.ionocytes, "cell.groups"] <- "ionocytes"
obj.merge@meta.data[cells.best4, "cell.groups"] <- "best4+_cells"
mama@meta.data[setdiff(WhichCells(mama), unlist(unique(list(cells.absorptive, cells.secretory, cells.ionocytes, cells.best4)))), "cell.groups"] <- "rest"

obj.intestine@meta.data$cell.groups <- NA
obj.intestine@meta.data[cells.absorptive, "cell.groups"] <- "absorptive"
obj.intestine@meta.data[cells.secretory, "cell.groups"] <- "secretory"
obj.intestine@meta.data[cells.best4, "cell.groups"] <- "best4+_cells"
obj.intestine@meta.data[setdiff(WhichCells(obj.intestine), unlist(unique(list(cells.absorptive, cells.secretory, cells.best4)))), "cell.groups"] <- "rest"

Idents(obj.intestine) <- obj.intestine@meta.data$cell.groups


Idents(mama) <- mama@meta.data$cell.groups
genes.to.plot <- unlist(unique(list(genes.absorptive.only[1:10], genes.secretory.only[1:10], genes.iono.only[1:10])))
genes.to.plot <- unlist(unique(list(genes.abs.only[1:10], "foxi3b", "rap1gap", "urahb", "zgc:92275", "calm1b", "si:ch73-359m17.9", "adcyap1a", "penka", "agr2", "sec11a", "xbp1", "sec62", "fev", "mir375-2", genes.iono.only[1:10])))

pdf("absorptive_vs_secretory_vs_ionocytes_specific_genes.pdf", width = 12, height = 24)
DotPlot(mama, features = genes.to.plot, idents = c("absorptive", "secretory", "ionocytes"), dot.scale = 16, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

##Plot genes overlapping with best4+ cells and absorptive cells
##Subset absorptive marker dataframe based on the shared genes between best4 cells and absorptive
markers.absorptive.df <- m.abs[genes.best4.abs.common, ]
##Sort the gene lists based on the avg_log2FC
markers.absorptive.df <- markers.absorptive.df[order(-markers.absorptive.df$avg_log2FC), ]
genes.abs.best4.plot <- markers.absorptive.df$gene[1:30]

pdf("best4_vs_absorptive_shared_genes.pdf", width = 16, height = 32)
DotPlot(obj.merge, features = c("fabp1b.1", "fabp2", "ada", "adipor2", "plac8.1", "acox1", "acox3", "acaa2", "acsl4a", "akap7", "ca4b", "anpepb", "vil1", "crip1", "zgc:77748", "si:ch211-71m22.1", "sec14l8", "cers3a", "ap1s3a", "aifm4", "si:ch211-133l5.7", "zgc:198329", "zgc:112146", "zgc:136472", "slc7a8a", "atp1a1a.4", "stoml3b", "muc13b", "eps8l3a", "zgc:172079", "gpd1c", "pygb", "apoa4b.1", "apoa4b.2", "aqp8a.1", "aqp8a.2", "ace2"), dot.scale = 16, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

pdf("best4_vs_secretory_shared_genes.pdf", width = 16, height = 32)
DotPlot(obj.merge, features = c("adgrg4b", "tmem176l.2", "dll4", "calml4a", "krt18a.1", "si:busm1-57f23.1", "CR318588.3", "fhl1b", "aqp3a", "sh3bgrl2", "litaf", "BX908782.3", "si:ch211-202h22.9", "tmem98", "socs3a", "tmem59"), dot.scale = 16, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

pdf("best4_vs_ionocytes_shared_genes.pdf", width = 16, height = 32)
DotPlot(obj.merge, features = c("zgc:193726", "si:ch211-39i22.1", "ca2", "sstr5", "slc6a6b", "skap2", "si:ch211-147d7.5", "tpm1", "si:dkey-112a7.4", "lrba", "CABZ01113812.1"), dot.scale = 20, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

##Ok what do best4 cells secrete? 
genes.best4.mucin <- c("muc13b", "galnt12", "galnt2", "galnt6", "b3gnt7l", "b4galt6", "st6galnac1.1", "st3gal7")

DotPlot(obj1, features = genes.best4.mucin, idents = c("10", "11", "12", "13", "15", "16", "17", "18"), dot.scale = 20, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

# Create the Venn diagram
library(VennDiagram)
venn.diagram(
  x = list("best4+ cells" = genes.best4, "absorptive" = genes.abs.only, "secretory" = genes.sec.only, "ionocytes" = genes.iono.only),
  filename = "best4_gene_overlapping_abs_sec_iono",
  col = "transparent",
  fill = c("orange", "skyblue1", "pink1", "mediumorchid1"),
  alpha = 0.5,
  cex = 1.5,
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  cat.col = "black",
  cat.pos = 0,
  margin = 0.05,
  main = "best4 genes overlap with abs/sec/iono"
)

gene.expression <- as.matrix(obj.intestine@assays$RNA@data)

# Subset the matrix to include only the columns (clusters) of interest
clusters <- levels(Idents(obj.intestine))


correlation_matrix <- cor(gene.expression, method = "pearson")
saveRDS(correlation_matrix, "~/Box/zfext/annotations_celltype_curated_newMama/endoderm/best4_markers/intestine_correlation_matrix.rds")
saveRDS(obj.intestine, "~/Box/zfext/annotations_celltype_curated_newMama/endoderm/best4_markers/intestine_seurat.rds")

library(corrplot)
corrplot(correlation_matrix, 
         method = "color",
         type = "upper",
         tl.cex = 0.7,
         addCoef.col = "black",
         number.cex = 0.7,
         title = "Correlation Plot")

