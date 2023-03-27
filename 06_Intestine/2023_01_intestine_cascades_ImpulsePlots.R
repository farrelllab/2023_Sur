library(Seurat)
library(URD)

#Specify paths
save.path <- "~/Box/zfext/results/06b-URD/cascades/intestine"
plot.path <- "~/Box/zfext/results/06b-URD/cascades/intestine/spline_heatmaps/"

sample = "intestine"

#Load URD tree object
obj.tree <- readRDS(file="~/Box/zfext/results/06b-URD/obj_subsets/2022-11_intestine_trajectory/obj/intestine_TREE.rds")
obj.build <- readRDS(file=paste0("~/Box/zfext/results/06b-URD/FDL/", sample, "/obj/", sample, "_urd_treeView2_modified_final.rds"))

#Differential expression with precision-recall along URD dendrogram
#Genes are considered differentially expressed if they are expressed in atleast 10% of cells in the trajectory segment under consideration and their mean expression was upregulated 1.5X compared to the sibling and the gene was 1.25x better than a random classifier for the population as defined by the area under a precision-recall curve

#Determine tips to run DE for
tips.to.run <- as.character(obj.tree@tree$segment.names)
genes.use <- NULL #Calculate for all genes

#Calculate the markers for each of the populations
gene.markers <- list()
for(tipn in 1:length(tips.to.run)) {
  tip <- tips.to.run[[tipn]]
  print(paste0(Sys.time(), ":", tip))
  markers <- aucprTestAlongTree(obj.tree, pseudotime = "pseudotime", tips = tip, log.effect.size = 0.4, auc.factor = 0.6, max.auc.threshold = 0.85, frac.must.express = 0.1, frac.min.diff = 0, genes.use = genes.use, root = "16", only.return.global = F, must.beat.sibs = 0.6, report.debug = T)
  saveRDS(markers, paste0(save.path, "/aucpr_markers/", tip, ".rds"))
  gene.markers[[tip]] <- markers
}

saveRDS(gene.markers, paste0(save.path, "/aucpr_markers/", "intestine_gene_markers_allTips.rds"))
gene.markers <- readRDS(paste0(save.path, "/aucpr_markers/", "intestine_gene_markers_allTips.rds"))

# Separate actual marker lists from the stats lists
gene.markers.de <- lapply(gene.markers, function(x) x[[1]])
gene.markers.stats <- lapply(gene.markers[1:8], function(x) x[[2]])
names(gene.markers.de) <- names(gene.markers)
names(gene.markers.stats) <- names(gene.markers)

## Examine the relationsip between DE genes and library complexity
# Compile all comparison stats into a single table
all.de.stats <- do.call("rbind", gene.markers.stats)
all.de.stats$tip <- substr(rownames(all.de.stats),1,nchar(rownames(all.de.stats))-2)
# Do a few plots
p1 <- ggplot(all.de.stats, aes(x=pt.1.mean, y=pt.2.mean)) + geom_point() + theme_bw() + geom_abline(slope = 1, intercept=0, col='red', lty=2) + labs(x="Mean Pseudotime (Group 1)", y="Mean Pseudotime (Group 2)")
p2 <- ggplot(all.de.stats, aes(x=genes.1.mean, y=genes.2.mean)) + geom_point() + theme_bw() + geom_abline(slope = 1, intercept=0, col='red', lty=2) + labs(x="Mean Detected Genes (Group 1)", y="Mean Detected Genes (Group 2)")
p3 <- ggplot(all.de.stats, aes(x=trans.1.mean, y=trans.2.mean)) + geom_point() + theme_bw() + geom_abline(slope = 1, intercept=0, col='red', lty=2) + labs(x="Mean Transcripts (Group 1)", y="Mean Transcripts (Group 2)")
cowplot::plot_grid(p1,p2,p3, ncol = 3)


# Determine temporal gene expression with impulse fitting
#Impulse fits
gene.cascades <- lapply(tips.to.run, function(tip) {
  print(paste0(Sys.time(), ": Impulse Fit ", tip))
  seg.cells <- cellsAlongLineage(obj.tree, tip, remove.root=F)
  casc <- geneCascadeProcess(object = obj.tree, pseudotime='pseudotime', cells = seg.cells, 
                             genes= rownames(gene.markers.de[[tip]]), 
                             moving.window=5, cells.per.window=18, 
                             pseudotime.per.window = 0.01, verbose = T, verbose.genes = T)
  return(casc)
})
names(gene.cascades) <- tips.to.run
saveRDS(gene.cascades, file = "./data/allGene.cascades.rds")


##As the tuft-cell branch comprises very few cells, try doing this manually; not sure this will work
cascade.tuft <- geneCascadeProcess(object = obj.tree, pseudotime = 'pseudotime', cells = cellsAlongLineage(obj.tree, "9", remove.root = F), 
                                   genes = rownames(gene.markers.de[[4]]),
                                   moving.window = 2, cells.per.window = 1,
                                   pseudotime.per.window = 0.01, verbose = T, verbose.genes = T)

saveRDS("~/Box/zfext/results/06b-URD/cascades/intestine/obj_cascades/intestine_tuft-like_cascade.rds")


cascade.B4O2<- geneCascadeProcess(object = obj.tree, pseudotime = 'pseudotime', cells = cellsAlongLineage(obj.tree, "5", remove.root = F), 
                                   genes = c("atoh1b", "ascl1a", "tnfrsf11a", "dacha", "pbx3a", "best4", "otop2", "cftr", "meis1b", "tox2", "sctr", "tacr2"),
                                   moving.window = 5, cells.per.window = 18,
                                   pseudotime.per.window = 0.01, verbose = T, verbose.genes = T)

geneCascadeImpulsePlots(cascade.B4O2, genes = c("tnfrsf11a", "dacha", "pbx3a", "best4", "otop2", "cftr", "ascl1a", "atoh1b"))

#Plotting gene expression
genes.plot <- c("tnfrsf11a", "meis1b", "dacha", "pbx3a", "best4", "otop2", "cftr", "atoh1b", "ascl1a")
plotSmoothFit(smoothed.fit = cascade.B4O2, genes = genes.plot, scaled = T, multiplot = F, alpha.data = 0.2, alpha.smooth = 1.6, lwd.smooth = 1.4)

#Load endoderm best4/otop2 cascade object
cascade.bo4 <- readRDS("~/Box/zfext/results/06b-URD/cascades/intestine/obj_cascades/intestine_B4O2_cascade.rds")
geneCascadeImpulsePlots(cascade.bo4, genes = c("tnfrsf11a", "dacha", "pbx3a", "best4", "otop2", "cftr", "ascl1a", "atoh1b"))
png("ImpulsePlot_cftr.png", width = 2*dpi, height = 2*dpi)
geneCascadeImpulsePlots(cascade.bo4, genes = c("atoh1b"))
dev.off()

cascade.lre <- readRDS("~/Box/zfext/results/06b-URD/cascades/intestine/obj_cascades/intestine_LREs_cascade.rds")
geneCascadeImpulsePlots(cascade.lre, genes = c("atoh1b", "tfeb", "amn", "dab2", "lgmn", "cubn"))
geneCascadeImpulsePlots(cascade.lre, genes = c("atoh1b", "foxd2", "cdx1a", "mafa", "mafbb", "tfeb","satb2","cubn", "amn", "lgmn", "dab2", "skilb"))

cascade.ec1 <- readRDS("~/Box/zfext/results/06b-URD/cascades/intestine/obj_cascades/intestine_EC-1_cascade.rds")
fire.with.grey <- c("#CECECE", "#DDC998", RColorBrewer::brewer.pal(9, "YlOrRd")[3:9])
pond.with.grey <- c("#CECECE", "#CBDAC2", RColorBrewer::brewer.pal(9, "YlGnBu")[3:9])
gene.num <- nrow(cascade.bo4$scaled.expression)
cols <- (scales::gradient_n_pal(RColorBrewer::brewer.pal(9, "YlOrRd")))(seq(0, 1,
                                                                            length.out = 50))

tf.list <- read.delim(file="~/Box/zfext/02-Clustering/2021-03 Iterative Clustering/gene_info/tfs/2021-06-25_zebrafish_LTA_TFs.txt", header = T, sep = "")
gene.num <- nrow(cascade.bo4$scaled.expression)
genes <- rownames(cascade.bo4$scaled.expression)
bo4.tfs <- genes[which(genes %in% tf.list$Symbol)]
bo4.ph <- c("tnfrsf11a", "best4", "otop2", "sctr", "pbx3a", "tacr2", "cftr", "gucy2c")

anno <- list(green=intersect(bo4.tfs, genes),
             blue=intersect(bo4.ph, genes))

dpi <- 300
pdf("best4_otop2_gene_cascade_2_red_yellow.pdf", w = 28, h = 22)
geneCascadeHeatmap(cascade.bo4, annotation.list = anno, row.font.size = 0.002*gene.num, color.scale = fire.with.grey)
dev.off()

gene.num <- nrow(cascade.bo4$scaled.expression)
genes <- rownames(cascade.bo4$scaled.expression)
bo4.tfs <- genes[which(genes %in% tf.list$Symbol)]
bo4.ph <- c("best4", "otop2", "sctr", "tacr2", "cftr", "gucy2c")

anno <- list(green=intersect(bo4.tfs, genes),
             blue=intersect(bo4.ph, genes))

##Only focus on TFs
tf.list <- read.delim(file="~/Box/zfext/02-Clustering/2021-03 Iterative Clustering/gene_info/tfs/2021-06-25_zebrafish_LTA_TFs.txt", header = T, sep = "")

gene.num <- nrow(cascade.lre$scaled.expression)
cols <- (scales::gradient_n_pal(RColorBrewer::brewer.pal(9, "YlOrRd")))(seq(0, 1,
                                                                            length.out = 50))
gene.num <- nrow(cascade.lre$scaled.expression)
genes <- rownames(cascade.lre$scaled.expression)
lre.tfs <- genes[which(genes %in% tf.list$Symbol)]
lre.genes <- c("cubn", "dab2", "amn", "lgmn", "hexa", "dpp", "lrp2b")
anno <- list(green=intersect(lre.tfs, genes),
             blue=intersect(lre.genes, genes))

##Compare the TFs with those captured from Seurat analysis
markers <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/endoderm/obj_seurat/endoderm_markers.rds")
genes.lre <- rownames(markers$wilcox$`13`)[which(markers$wilcox$`13`$pct.2 <= 0.2 & markers$wilcox$`13`$avg_log2FC >= 0.25)]
genes.lre.tfs <- genes.lre[which(genes.lre %in% tf.list$Symbol)]

genes.prog <- rownames(markers$wilcox$`14`)[which(markers$wilcox$`14`$pct.2 <= 0.3 & markers$wilcox$`14`$avg_log2FC >= 0.25)]
genes.prog.tfs <- genes.prog[which(genes.prog %in% tf.list$Symbol)]

intersect(genes.prog.tfs, genes.lre.tfs)

dpi <- 300
pdf("post_LRE_gene_cascade.pdf", w = 28, h = 22)
geneCascadeHeatmap(cascade.lre, annotation.list = anno, row.font.size = 0.002*gene.num)
dev.off()

##Plot spline curves for LREs
plotSmoothFit(smoothed.fit = cascade.lre, genes = c("atoh1b", "foxd2", "cdx1a", "mafa", "mafbb", "tfeb","satb2","cubn", "amn", "lgmn", "dab2", "skilb"), scaled = T, multiplot = F,alpha.data = 0.2, alpha.smooth = 1.6, lwd.smooth = 1.4)
plotSmoothFit(smoothed.fit = cascade.lre, genes = c("cubn", "amn", "lgmn", "dab2", "cdx1a", "cdx1b", "skilb"), scaled = T, multiplot = F, alpha.data = 0.2, alpha.smooth = 1)
splines <- list(cascade, cascade.lre)
names(splines) <- c("best4+/otop2+", "post-LREs")

plotSmoothFit(smoothed.fit = cascade.lre, genes = lre.tfs[220:240], scaled = T, multiplot = F, alpha.data = 0.2, alpha.smooth = 1)

#So the question is whether best4/otop2+ cells originate from the atoh1+ progenitors
##First plot a UMAP plot with the cells in question
plotTree(obj.tree, label.segments = T)
gridExtra::grid.arrange(grobs = lapply(c("atoh1b", "tnfrsf11a", "ascl1a", "pou2f3"), plotTree,
                                       object = obj.tree, label.x = F, plot.cells = T), ncol = 2)

cells.in.post.prog <- cellsInCluster(obj.tree, "segment", c("39", "34", "29", "16", "13"))
obj.tree <- groupFromCells(obj.tree, group.id = "cells.in.post.prog", cells = cells.in.post.prog)
plotTreeHighlight(obj.tree, "cells.in.post.prog", "TRUE", highlight.size = 1, highlight.alpha = 0.5)

cells.prog <- WhichCells(obj1, idents = c("24"))
cells.bo4 <- WhichCells(obj1, idents = c("29"))
cells.gob <- WhichCells(obj1, idents = c("31"))
cells.lre <- WhichCells(obj1, idents = c("15"))
cells.tuft <- WhichCells(obj1, idents = c("35"))
cells.eec <- WhichCells(obj1, idents = c("32"))

obj1@meta.data$cells.post <- NA
obj1@meta.data[cells.prog, "cells.post"] <- "progenitor_2"
obj1@meta.data[cells.bo4, "cells.post"] <- "best4+/otop2+"
obj1@meta.data[cells.gob, "cells.post"] <- "goblet-cells"
obj1@meta.data[cells.lre, "cells.post"] <- "posterior_LREs"
obj1@meta.data[cells.tuft, "cells.post"] <- "tuft-like"
obj1@meta.data[cells.eec, "cells.post"] <- "EECs"

post.cols <- c("#a84b99", "#FF7F50", "#32CD32", "#48D1CC", "#85AC00", "#A0522D")
post.cells <- c("goblet-cells", "progenitor_2", "best4+/otop2+", "posterior_LREs", "EECs", "tuft-like")
colors <- setNames(post.cols, post.cells)
DimPlot(obj1, group.by = "cells.post", cols = colors) 

png("endoderm_posterior_cells.png", width = 5*dpi, height = 5*dpi)
#Labeling the clusters using the same color as the clusters
DimPlot(obj1, group.by = "cells.post", cols = colors, pt.size = 1) + theme(legend.position="none") + NoAxes()
dev.off()

##Find differentially expressed markers for tuft-cells
markers.tuft <- FindMarkers(obj1, ident.1 = "35", ident.2 = setdiff(levels(Idents(obj1)), "35"), logfc.threshold = 0.25, min.diff.pct = 0.25)

png("endoderm_gene_exp.png", width = 5*dpi, height = 5*dpi)
FeaturePlot(obj1, features = c("atoh1b", "tnfrsf11a", "ascl1a", "sox4b"), cols = c("lightgrey", "red"), pt.size = 1.2) + theme(legend.position="none") + NoAxes()
FeaturePlot(obj1, features = c("best4", "otop2", "cftr", "tacr2"), cols = c("lightgrey", "red"), pt.size = 1.2)  + theme(legend.position="none") + NoAxes()
FeaturePlot(obj1, features = c("agr2", "muc2.1", "fmn1", "fer1l6"), cols = c("lightgrey", "red"), pt.size = 1.2)  + theme(legend.position="none") + NoAxes()
FeaturePlot(obj1, features = c("amn", "dab2", "lgmn", "cubn"), cols = c("lightgrey", "red"), pt.size = 1.2)  + theme(legend.position="none") + NoAxes()
FeaturePlot(obj1, features = c("pou2f3", "sox8b", "tnfrsf11b", "rgs13"), cols = c("lightgrey", "red"), pt.size = 0.6)  + theme(legend.position="none")
dev.off()

plotSmoothFit(smoothed.fit = cascade.bo4, genes = c("tox2", "tnfrsf11a", "meis1b", "ascl1a", "dacha", "pbx3a", "best4", "otop2", "cftr"), scaled = T, multiplot = F, alpha.data = 0.2, alpha.smooth = 1.6, lwd.smooth = 1.4)
genes.plot <- c("tox2", "tnfrsf11a", "pbx3a", "dacha", "best4", "otop2", "ascl1a", "mafa", "tfeb", "cubn", "amn")
plotSmoothFitMultiCascade(smoothed.fit = splines, genes = genes.plot, scaled = T, alpha.data = 0.2, alpha.smooth = 1)

geneCascadeImpulsePlots(cascade, genes = genes.plot, ymin0 = TRUE)

cascade.lre <- readRDS("~/Box/zfext/results/06b-URD/cascades/endoderm/obj_cascades/endoderm_post-LREs_cascade.rds")
plotSmoothFit(smoothed.fit = cascade.lre, genes = c("atoh1b", "ascl1a", "mafa", "tfeb", "cubn", "amn"), scaled = T, multiplot = F, alpha.data = 0.2, alpha.smooth = 1)
splines <- list(cascade, cascade.lre)
names(splines) <- c("best4+/otop2+", "post-LREs")

plotSmoothFit(smoothed.fit = cascade.lre, genes = c("atoh1b", "ascl1a", "mafa", "tfeb", "cubn", "amn", "ctsz", "hexa", "dpp", "lgmn"), scaled = T, multiplot = F, alpha.data = 0.2, alpha.smooth = 1)

##Get cascades for endocrine pancreas and EECs
cascade.panc <- readRDS("~/Box/zfext/results/06b-URD/cascades/endoderm/obj_cascades/endoderm_endo-panc_cascade.rds")
cascade.eec <- readRDS("~/Box/zfext/results/06b-URD/cascades/endoderm/obj_cascades/endoderm_EECs_cascade.rds")
cascades <- list(cascade.panc, cascade.eec)
names(cascades) <- c("endo_panc", "EECs")
genes.plot <- c("lmx1ba", "atoh1b", "penka", "pax6b", "pax4", "neurod1")
plotSmoothFit(smoothed.fit = cascade.panc, genes = genes.plot, scaled = T, multiplot = F, alpha.data = 0.2, alpha.smooth = 1)
plotSmoothFit(smoothed.fit = cascade.eec, genes = genes.plot, scaled = T, multiplot = F, alpha.data = 0.2, alpha.smooth = 1)

plotSmoothFitMultiCascade(smoothed.fit = cascades, genes = genes.plot, scaled = T, alpha.data = 0.2, alpha.smooth = 1)
gene.endo.order <- rownames(cascade.panc$timing)[order(cascade.panc$timing$time.on)]

gene.num <- nrow(cascade.panc$scaled.expression)
cols <- (scales::gradient_n_pal(RColorBrewer::brewer.pal(9, "YlOrRd")))(seq(0, 1,
                                                                            length.out = 50))
dpi <- 300
pdf("endo_panc_cascade.pdf", w = 12, h = 0.08*gene.num)
gplots::heatmap.2(x = as.matrix(cascade.panc$scaled.smooth), Rowv = F,
                  Colv = F, dendrogram = "none", trace = "none", density.info = "none",
                  key = F, cexCol = 0.8, cexRow = 0.15, margins = c(8, 8), lwid = c(0.3, 4), lhei = c(0.3,
                                                                                                      4), labCol = NA)
dev.off()

anno <- list(c(blue = endo.tfs))
dpi <- 300
pdf("endo_panc_cascade.pdf", w = 20, h = 22)
geneCascadeHeatmap(cascade.panc, color.scale = cols, annotation.list = anno, row.font.size = 0.008*gene.num)
dev.off()

geneCascadeHeatmap <- function(cascade, color.scale=RColorBrewer::brewer.pal(9, "YlOrRd"), add.time=NULL, times.annotate=seq(0,1,0.1), title="", annotation.list=NULL, row.font.size=1, max.font.size=0.9) {
  # Correct for NA timings
  timing <- cascade$timing
  timing[intersect(which(is.na(timing$time.on)), which(is.infinite(timing$time.off))), "time.on"] <- Inf
  gene.order <- order(timing$time.on, timing$time.off, na.last=F)
  cols <- scales::gradient_n_pal(color.scale)(seq(0,1,length.out = 50))
  if (!is.null(add.time)) {
    time <- unlist(lapply(cascade$pt.windows, function(cells) mean(object@meta[cells, add.time])))
  } else {
    time <- as.numeric(names(cascade$pt.windows))
  }
  time.lab <- rep("", length(time))
  for (annotate in times.annotate) {
    gt <- which(time >= annotate)
    if (length(gt)>0) time.lab[min(gt)] <- as.character(annotate)
  }
  if (!is.null(annotation.list)) {
    annot <- data.frame(
      gene=rownames(cascade$timing)[gene.order],
      type=factor(rep(NA, length(gene.order)), levels=unique(names(annotation.list))),
      row.names=rownames(cascade$timing)[gene.order],
      stringsAsFactors=F
    )
    for (l in names(annotation.list)) {
      annot[annotation.list[[l]], "type"] <- l
    }
    gplots::heatmap.2(as.matrix(cascade$scaled.expression[gene.order,]), Rowv=F, Colv=F, dendrogram="none", col=cols, trace="none", density.info="none", key=F, labCol=time.lab, RowSideColors=as.character(annot$type), cexCol=0.8, cexRow=0.15, margins = c(8,8), lwid=c(0.3,4), lhei=c(0.3, 4))
  } else {
    gplots::heatmap.2(as.matrix(cascade$scaled.expression[gene.order,]), Rowv=F, Colv=F, dendrogram="none", col=cols, trace="none", density.info="none", key=F, labCol=time.lab, cexCol=0.8, cexRow=0.7, margins = c(8,8), lwid=c(0.3,4), lhei=c(0.3, 4))
  }
  title(title, line=1, adj=0.4)
}

geneCascadeImpulsePlots(cascade.lre, genes = c("ascl1a", "satb2", "hoxc13a", "amn"))
png("ImpulsePlot_amn.png", width = 2*dpi, height = 2*dpi)
geneCascadeImpulsePlots(cascade.lre, genes = c("amn"))
dev.off()

cascade.lre <- readRDS("~/Box/zfext/results/06b-URD/cascades/endoderm/obj_cascades/endoderm_goblet-cells_cascade.rds")
geneCascadeImpulsePlots(cascade.lre, genes = c("sox4a", "atoh1b", "sall4", "agr2"))
png("ImpulsePlot_agr2.png", width = 2*dpi, height = 2*dpi)
geneCascadeImpulsePlots(cascade.lre, genes = c("agr2"))
dev.off()

png("best4_otop2_spline_curve.png", width = 5 * dpi, height = 5*dpi)
plotSmoothFit(smoothed.fit = cascade.bo4, genes = c("tnfrsf11a", "dacha", "pbx3a", "best4", "otop2", "cftr"), scaled = T, multiplot = F, alpha.data = 0.2, alpha.smooth = 1.2, lwd.smooth = 1.6)
dev.off()



########################### SPLINE CURVES ALONG URD TRJAECTORY BRANCHES ########################################
## Let's look at some of the spline curves of the tuft-like genes and other cell types too. 
#Calculate smoothed spline fits of their expression
tip.to.use <- "B4O2"
spline <- geneSmoothFit(obj.build, method = "spline", pseudotime = "pseudotime", cells = cellsAlongLineage(obj.build, tip.to.use), genes = c(rownames(gene.markers$B4O2$diff.exp), "atoh1b", "tnfrsf11a"), moving.window = 2, cells.per.window = 5, spar = 0.9)

#Finding dynamic genes who change their expression values considerably
#Which genes change their actual mean expression value by at least 0.75
change.real <- apply(spline$mean.smooth, 1, function(x) diff(range(x)))

genes.mean <- names(which(change.real >= 0.75))

#Which genes are well fit by the spline curves? (Noise is poorly fit)
spline.fit <- apply(spline$scaled.smooth - spline$scaled.expression.red, 1, function(i) sum(i^2))

spline.fit.norm <- spline.fit/ncol(spline$scaled.smooth)
wellfit <- names(which(spline.fit.norm <= 0.01))

#Which genes change their log2 mean expression value sufficiently? At least 33% and requiring more change the more poorly fit the data is by its spline curve.
change.scale <- apply(spline$scaled.smooth, 1, function(x) diff(range(x)))

genes.scale <- names(which(change.scale >= spline.fit.norm * 0.15/0.02 + 0.33))

#Ensure that genes are fit by the spline curve significantly better than a flat line of slope 0. (Weighted by distance to next point in pseudotime to compensate for point density)
w <- ncol(spline$scaled.smooth) - 1
weight <- diff(as.numeric(colnames(spline$scaled.smooth))) * 1000
spline.fit.weighted <- apply(spline$scaled.smooth[, 1:w] - spline$scaled.expression.red[, 1:w], 1, function(i) sum(weight * i^2))

flat.fit.weighted <- apply(spline$scaled.expression.red[, 1:w], 1, function(x) sum(weight * (x - mean(x)) ^ 2))
spline.fit.ratio <- log2(flat.fit.weighted/spline.fit.weighted)
spline.fit.betterthanflat <- names(which(spline.fit.ratio > 0.25))

#Take intersection of these genes and use them as the "varying genes" in the heatmap and analysis
genes <- intersect(intersect(intersect(genes.scale, genes.mean), wellfit), spline.fit.betterthanflat)

##Combine pieces of spline fits into a single list of multi-plotting identify the pesudotime of the branchpoint
pt.crop <- as.numeric(unlist(obj.tree@tree$segment.pseudotime.limits)[1])
#Crop according to the pseudotime of the branchpoint
only.spline <- cropSmoothFit(spline, pt.min = 0.1)

splines <- list(only.spline)
names(splines) <- "B4O2"

#Plotting gene expression
genes.plot <- c("tnfrsf11a", "tox2", "dacha", "pbx3a", "best4", "otop2", "cftr", "atoh1b", "ascl1a")

plotSmoothFit(smoothed.fit = spline, genes = genes.plot, scaled = T, multiplot = F, alpha.data = 0.2, alpha.smooth = 1)
plotSmoothFitMultiCascade(smoothed.fit = splines, genes = genes.plot, scaled = T, alpha.data = 0.2, alpha.smooth = 1)

# Reduce the data in the three splines, such that each block is minimum 0.01
# pseudotime; in this case, each block will be >= 5 cells, >= pseudotime 0.01
splines.reduced <- lapply(names(splines), function(n) {
  s <- spline
  l <- tolower(substr(n, 1, 1))
  colnames(s$mean.expression) <- paste0(l, as.character(round(as.numeric(colnames(s$mean.expression)),
                                                              digits = 2)))
  s$mean.expression <- matrixReduce(s$mean.expression)
  colnames(s$mean.smooth) <- paste0(l, as.character(round(as.numeric(colnames(s$mean.smooth)), digits = 2)))
  s$mean.smooth <- matrixReduce(s$mean.smooth)
  return(s)
})

# Combine the spline into a table & rescale maximum across all regions
me <- do.call("cbind", lapply(splines.reduced, function(s) s$mean.expression))
ms <- do.call("cbind", lapply(splines.reduced, function(s) s$mean.smooth))
se <- sweep(me, 1, apply(me, 1, min), "-")
se <- sweep(se, 1, apply(se, 1, max), "/")
ss <- sweep(ms, 1, apply(ms, 1, min), "-")
ss <- sweep(ss, 1, apply(ss, 1, max), "/")

# Get all genes
genes.plot <- rownames(cascade.bo4$scaled.expression)

#Filter heatmap genes
#Removes undesired (mitochondrial, ribosomal, tandem duplicated genes) from heatmaps for presentation purposes
#Genes (Character vector) genes to check
#Returns genes with undesired genes removed

filter.genes <- function(genes) {
  mt.genes <- grep("^mt-", ignore.case = T, genes, value = T)
  many.genes <- grep("\\(1 of many\\)", ignore.case = T, genes, value = T)
  cox.genes <- grep("^cox", ignore.case = T, genes, value = T)
  rp.genes <- grep("^rp", ignore.case = T, genes, value = T)
  hb.genes <- grep("^hb", ignore.case = T, genes, value = T)
  si.genes <- grep("^si:", ignore.case = T, genes, value = T)
  zgc.genes <- grep("^zgc:", ignore.case = T, genes, value = T)
  return(setdiff(genes, c(mt.genes, many.genes, cox.genes, rp.genes, hb.genes, si.genes, zgc.genes)))
}

genes.plot <- filter.genes(genes.plot)

# Hierarchical cluster based on smoothed expression
h.ss <- hclust(dist(as.matrix(ss[genes.plot, ])), method = "ward.D2")
h.ss.d <- as.dendrogram(h.ss)
k = 6 # Cluster number chosen by eye.
h.ss.clust <- cutree(h.ss, k = k)
# Get cluster order as it will be in the heatmap
clust.order <- unique(h.ss.clust[h.ss$order])
h.ss.clust.ord <- plyr::mapvalues(from = clust.order, to = 1:k, x = h.ss.clust)
# Generate cluster color vector
cluster.h <- seq(0, 1, length.out = k + 1)
cluster.s <- rep(c(1, 0.75), length = k)
cluster.v <- rep(c(1, 0.75), length = k)
cluster.colors <- hsv(h = cluster.h[1:k], s = cluster.s, v = cluster.v)
h.ss.clust.col <- plyr::mapvalues(from = as.character(1:k), to = cluster.colors,
                                  x = h.ss.clust.ord)
cols <- (scales::gradient_n_pal(RColorBrewer::brewer.pal(9, "YlOrRd")))(seq(0, 1,
                                                                            length.out = 50))
# Generate the actual heatmap and save as a PDF.
pdf(paste0(plot.path, "best4_otop2_heatmap_spline.pdf"), width = 17, height = 22)
# Plot the heatmap
gplots::heatmap.2(as.matrix(ss[genes.plot[h.ss$order], ]), Rowv = F, RowSideColors = h.ss.clust.col[h.ss$order],
                  Colv = F, dendrogram = "none", col = cols, trace = "none", density.info = "none",
                  key = F, cexCol = 0.6, cexRow = 0.3, margins = c(8, 8), lwid = c(0.3, 4), lhei = c(0.3,
                                                                                                   4), labCol = NA)
dataset <- "best4+/otop2+"  
# Put a title on it
title(dataset, line = -1, adj = 0.48, cex.main = 4)

dev.off()

