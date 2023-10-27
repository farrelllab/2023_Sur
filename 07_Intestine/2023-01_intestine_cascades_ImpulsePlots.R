##This code was written to calculate the sequence of gene expression changes along the different branches of the URD trajectory especially best4+ enterocytes and posterior lysosome-rich enterocytes. This code was used to generate figure panels: Figure 7F, I and Figure S10. 

##Load libraries
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


##Define geneCascadeHeatmap function
##geneCascadeHeatmap function in URD  was slightly modified. 
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


##As the tuft-cell branch comprises very few cells, try doing this manually; not sure this will work
cascade.tuft <- geneCascadeProcess(object = obj.tree, pseudotime = 'pseudotime', cells = cellsAlongLineage(obj.tree, "9", remove.root = F), 
                                   genes = rownames(gene.markers.de[[4]]),
                                   moving.window = 2, cells.per.window = 1,
                                   pseudotime.per.window = 0.01, verbose = T, verbose.genes = T)

saveRDS("~/Box/zfext/results/06b-URD/cascades/intestine/obj_cascades/intestine_tuft-like_cascade.rds")

## Fit  the gene expression cascades to an impulse model for the best4 enterocytes
cascade.B4O2<- geneCascadeProcess(object = obj.tree, pseudotime = 'pseudotime', cells = cellsAlongLineage(obj.tree, "5", remove.root = F), 
                                   genes = rownames(gene.markers.de[[3]]),
                                   moving.window = 5, cells.per.window = 18,
                                   pseudotime.per.window = 0.01, verbose = T, verbose.genes = T)

##Plot the Impulse plots - Figure 7F
geneCascadeImpulsePlots(cascade.B4O2, genes = c("tnfrsf11a", "dacha", "pbx3a", "best4", "otop2", "cftr", "ascl1a", "atoh1b", "sox4a", "sox4b"))

#Plotting gene expression
genes.plot <- c("tnfrsf11a", "meis1b", "dacha", "pbx3a", "best4", "otop2", "cftr", "atoh1b", "ascl1a", "sox4a", "sox4b", "cbx3b", "onecut1")
plotSmoothFit(smoothed.fit = cascade.B4O2, genes = genes.plot, scaled = T, multiplot = F, alpha.data = 0.2, alpha.smooth = 1.6, lwd.smooth = 1.4)

#Load endoderm best4/otop2 cascade object
cascade.bo4 <- readRDS("~/Box/zfext/results/06b-URD/cascades/intestine/obj_cascades/intestine_B4O2_cascade.rds")
geneCascadeImpulsePlots(cascade.bo4, genes = c("tnfrsf11a", "dacha", "pbx3a", "best4", "otop2", "cftr", "ascl1a", "atoh1b", "sox4a", "sox4b", "cbx3b", "onecut1"))
png("ImpulsePlot_cftr.png", width = 2*dpi, height = 2*dpi)
geneCascadeImpulsePlots(cascade.bo4, genes = c("atoh1b"))
dev.off()


##load the transcription factor list 
tf.list <- read.delim(file="~/Box/zfext/02-Clustering/2021-03 Iterative Clustering/gene_info/tfs/2021-06-25_zebrafish_LTA_TFs.txt", header = T, sep = "")
gene.num <- nrow(cascade.bo4$scaled.expression)
genes <- rownames(cascade.bo4$scaled.expression)
bo4.tfs <- genes[which(genes %in% tf.list$Symbol)]
bo4.ph <- c("tnfrsf11a", "best4", "otop2", "sctr", "pbx3a", "tacr2", "cftr", "gucy2c")

##Add bars  to annotate which genes are transcription factors
anno <- list(green=intersect(bo4.tfs, genes),
             blue=intersect(bo4.ph, genes))

##Plot best4+ cell gene expression cascade as a heatmap - Figure S10
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

##Fit gene cascades to impulse plots for lysosome-rich enterocytes - Figure 7I
cascade.lre <- readRDS("~/Box/zfext/results/06b-URD/cascades/intestine/obj_cascades/intestine_LREs_cascade.rds")
geneCascadeImpulsePlots(cascade.lre, genes = c("atoh1b", "tfeb", "amn", "dab2", "lgmn", "cubn"))
geneCascadeImpulsePlots(cascade.lre, genes = c("atoh1b", "foxd2", "cdx1a", "mafa", "mafbb", "tfeb","satb2","cubn", "amn", "lgmn", "dab2", "skilb"))

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

##Plot full LRE cascade as a heatmap - not shown
dpi <- 300
pdf("post_LRE_gene_cascade.pdf", w = 28, h = 22)
geneCascadeHeatmap(cascade.lre, annotation.list = anno, row.font.size = 0.002*gene.num)
dev.off()







