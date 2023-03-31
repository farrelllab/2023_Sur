#This code is directed towards using the endoderm tissue-specific atlas to understand differential gene expression across specific cell types. This code is used to generate the following figure panels: Figure 5A-C, E-F" and Figure S10A, B

##Load libraries
library(Seurat)
library(URD)

obj1 <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/endoderm/obj_seurat/endoderm_seurat.rds")
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

##Plot  the endoderm UMAP plot colored by stage (hpf) - Figure S10A
dpi <- 300
png(file = "endoderm_stage.nice.png", width = 5 * dpi, height = 5 * dpi)
DimPlot(obj1, group.by = "stage.nice", cols = stage.colors.new, pt.size = 1.2)
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

all.colors <- c(liver.cols, intestine.cols, endo.panc.cols, exo.panc.cols, rest.cols)

names(cluster_ids) <- levels(obj1)
obj.endo <- RenameIdents(obj1, cluster_ids)
DimPlot(obj.endo, label = T)

DimPlot(obj.endo, label = F, cols = all.colors)


##Plot the differentially expressed transcription factors expressed between broader cell types within the endodermal atlas
##Plotting the broader clusters in UMAP
endo.organ <- c("liver", "intestine", "exocrine pancreas", "endocrine pancreas", "progenitor_1", "progenitor_2", "progenitor_3", "progenitor_4", "pharynx", "esophagus", "cloaca", "hepatopancreatic progenitors", "sftpba+_cells")
endo_colors <- c("#DAA520", "#00BFFF", "#9932CC", "#FF6347", "#00CED1", "#FF1493", "#008000", "#8FBC8F", "#800000", "#00FA9A", "#FF00FF", "#008080", "#DB7093")
endo.tissue.cols <- setNames(endo_colors, endo.organ)
DimPlot(obj2, label = F, cols = endo.tissue.cols)

dpi <- 300
png("endoderm_UMAP_broad_clusters.png", width = 5 * dpi, height = 5 * dpi)
DimPlot(obj2, label = F, cols = endo.tissue.cols, pt.size = 1.2) + theme(legend.position="none") + NoAxes()
dev.off()

dpi <- 300
png("endoderm_UMAP_stage.group.png", width = 5 * dpi, height = 5 * dpi)
DimPlot(obj1, group.by = "stage.group", label = F, pt.size = 1.2) + theme(legend.position="none") + NoAxes()
dev.off()

#Labeling the clusters using the same color as the clusters
p <- DimPlot(obj2, label = F, cols = endo.tissue.cols) + theme(legend.position="none") + NoAxes()
LabelClusters(p, id = "ident", color = unique(ggplot_build(p)$data[[1]]$colour))

##Plotting a dotplot between the major cluster groups in the endoderm
##For that first I need to do FindMarkers
markers.endo <- FindAllMarkers(obj2, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

##Only focus on TFs
tf.list <- read.delim(file="~/Box/zfext/02-Clustering/2021-03 Iterative Clustering/gene_info/tfs/2021-06-25_zebrafish_LTA_TFs.txt", header = T, sep = "")
markers.tfs.only <- markers.endo[which(rownames(markers.endo) %in% tf.list$Symbol), ]
#Take only top 5 of the TFs expressed in the different endoderm cell types
top5 <- markers.tfs.only %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
gene.list <- unlist(top5$gene)

#Plot DotPlot
p <- DotPlot(obj2, features = gene.list) + coord_flip()

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
genes.plot <- c("sftpba", "sim1b", "ihha")

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
genes.to.use <- rownames(markers.endo$wilcox$`32`)
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
genes.to.plot <- rownames(markers.endo$wilcox$`16`)

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




