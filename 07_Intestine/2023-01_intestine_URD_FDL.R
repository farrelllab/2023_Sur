# This code is designed to generate a force-directed layout for the URD branching tree calculated on the intestinal cells. This code was used to generate Figures 7D-E

##Load libraries
library(Seurat)
library(URD)
library(rgl)

##For saving, pre-load the paths
save.path <- "~/Box/zfext/results/06b-URD/obj_subsets/2022-11_intestine_trajectory/obj/"
plot.path <- "~/Box/zfext/results/06b-URD/plot_subsets/plots/"

# Generate a force-directed layout
## Choose cells that were visited more robustly
urd.tree <- obj.tree
# Data frame to measure cell visitation
visitation <- data.frame(
  cell=rownames(urd.tree@diff.data),
  seg=urd.tree@diff.data$segment,
  stringsAsFactors=F, row.names=rownames(urd.tree@diff.data)
)

visitation$visit <- log10(apply(visitation, 1, function(cr) urd.tree@diff.data[as.character(cr['cell']), paste0("visitfreq.raw.", as.character(cr['seg']))])+1)
# Choose those cells that were well visited
robustly.visited.cells <- visitation[visitation$visit >= 0.5, "cell"]
# Since some tips of the tree were combined in their entirety, get the terminal segments to use as the tips of the force-directed layout.
final.tips <- segTerminal(urd.tree)

## Choose the number of nearest neighbors
#Here, we generate the force-directed layout by varying the number of nearest neighbors (num.nn).

par(mfrow=c(3,2))
for (k in c(80:90)){
  urd.tree <- treeForceDirectedLayout(urd.tree, num.nn = 87, cells.to.do = robustly.visited.cells, cut.outlier.cells = NULL,
                                      cut.outlier.edges = NULL,
                                      cut.unconnected.segments = 2, min.final.neighbors = 4,
                                      tips = final.tips, verbose = T)
  
  plotTreeForce(urd.tree, "stage.group", alpha=1)
}

##Optimal num.nn = 87


#Change the number of unconnected segments to find the optimal number of unconnected segments
for (j in c(2,4,6,8)){
  urd.tree <- treeForceDirectedLayout(urd.tree, num.nn = 87, method = "fr", cells.to.do = robustly.visited.cells,
                                      cut.outlier.cells = NULL,
                                      cut.outlier.edges = NULL, max.pseudotime.diff = NULL,
                                      cut.unconnected.segments = j, min.final.neighbors = 4,
                                      tips = final.tips, coords = "auto", start.temp = NULL,
                                      density.neighbors = 10, plot.outlier.cuts = F, verbose = T)
  
  plotTreeForce(urd.tree, "segment", alpha=1)
}

## optimal cut.unconnected.segments = 2

## Calculate the layout that is used in the paper.

# use the optimal parameters
urd.tree <- treeForceDirectedLayout(urd.tree, num.nn = 87, cells.to.do = robustly.visited.cells, cut.outlier.cells = NULL,
                                    cut.outlier.edges = NULL,
                                    cut.unconnected.segments = 2, min.final.neighbors = 4,
                                    tips = final.tips, verbose = T)
plotTreeForce(urd.tree, "segment", alpha= 1)


## Rotate the tree and save the view
urd.build <- plotTreeForceStore3DView(urd.tree, "View2")
saveRDS(urd.build, file=paste0(save.path, sample, "_urd_treeView2.rds"))
# Load previous saved object
urd.build <- readRDS(file=paste0(save.path, sample, "_urd_treeView2.rds"))
plotTreeForce(urd.build, "stage.group", alpha=1, alpha.fade=0.08, size=4, density.alpha=T, label.tips=T, view = "View2")
rglwidget()

#Hand tuning the tree for publication
########********* The idea is to have the segment 15 along with all its children to the left and segment 14 along with all its children on the right ********** 
#SEGMENT 15 - progenitor-2 branch and its derivatives
#First let's move the "B4O2, "EECs", and "goblet-cells" branches outwards and towards the left so that the 3 branches of the enterocytes are together on the right making more space in the center.
for(throw.out in c(0, 100, 250, 500, 1000)) {
  obj.build <- treeForceRotateCoords(urd.build, seg = "13", angle = -3.3, axis = "x", around.cell = 10, throw.out.cells = throw.out, pseudotime = "pseudotime")
  plotTreeForce(obj.build, "segment", alpha=1, alpha.fade=0.08, size=4, density.alpha=T, label.tips=T, view = "View2")
}
obj.build <- treeForceRotateCoords(obj.build, seg = "16", angle = 3.3, axis = "x", around.cell = 10, throw.out.cells = 0, pseudotime = "pseudotime")
plotTreeForce(obj.build, "segment", alpha=1, alpha.fade=0.08, size=4, density.alpha=T, label.tips=T, view = "View2")
obj.build <- plotTreeForceStore3DView(obj.build, "View2")
saveRDS(obj.build, file=paste0(save.path, sample, "_urd_tree_FDL_modified.rds"))
obj.build <- readRDS(file=paste0(save.path, sample, "_urd_tree_FDL_modified.rds"))

##Now we have brought the entire segment 15 to the front, now we have to move that to the left. I tried different angles but settled on -45 degrees. I also tried different number of "throw.out" cells to find at which range I get the best results. 
for(throw.out in c(0, 100, 250, 500, 1000)) {
  obj.build <- treeForceRotateCoords(obj.build, seg = "15", angle = -45, axis = "x", around.cell = 10, throw.out.cells = throw.out, pseudotime = "pseudotime")
  plotTreeForce(obj.build, "segment", alpha=1, alpha.fade=0.08, size=4, density.alpha=T, label.tips=T, view = "View2")
}

obj.build <- treeForceRotateCoords(obj.build, seg = "14", angle = 3.5, axis = "z", around.cell = 2, throw.out.cells = 0, pseudotime = "pseudotime")

##Try more number of throw out cells
for (to in seq(0,360,length.out = 10)) {
  obj.build <- treeForceRotateCoords(obj.build, seg = "15", angle = -4.5, axis = "x", around.cell = NULL, throw.out.cells = to, pseudotime = "pseudotime", around.by = "x")
}
plotTreeForce(obj.build, "segment", alpha=1, alpha.fade=0.08, size=4, density.alpha=T, label.tips=T, view = "View2")
plotTreeForce(obj.build, "pseudotime", alpha=1, alpha.fade=0.08, size=4, density.alpha=T, label.tips=T, view = "View2", seg.show = "45")

##Save the FDL in the current view
obj.build <- plotTreeForceStore3DView(obj.build, "View2")
saveRDS(obj.build, file=paste0(save.path, sample, "_urd_treeView2_modified_final.rds"))
obj.build <- readRDS(file=paste0(save.path, sample, "_urd_treeView2_modified_final.rds"))


##Plot some genes
fire.with.grey <- c("#CECECE", "#DDC998", RColorBrewer::brewer.pal(9, "YlOrRd")[3:9])
geneList <- c("tnfrsf11a", "atoh1b", "ascl1a", "sox4a", "sox4b", "pou2f3", "sox8b", "best4", "otop2", "cftr", "amn", "cubn", "muc2.1", "agr2")
geneList <- c("penka", "nmbb", "pax6b", "adcyap1b")
for (gene in geneList) {
  plotTreeForce(obj.build, gene, alpha=0.2, alpha.fade=0.1, size= 5, density.alpha=T, label.tips=F, view = "View2",
                colors = fire.with.grey)
  
}

plotTreeForce(obj.build, "stage.group", alpha=0.2, alpha.fade=0.08, size=4, density.alpha=T, label.tips=F, view = "View2")

##Use the stage color scheme as used in Figure 1
stage.colors.new <- c(
  colorRampPalette(c("#8b5500"))(1), # 10
  rep(RColorBrewer::brewer.pal(6, "Oranges"), each = 2), # 24, 26 / 28, 30 / 32, 34 / 36, 38 / 40, 42 / 44, 46
  colorRampPalette(c("#FFC0CB", "#FFB6C1", "#FF69B4", "#DA70D6", "#BA55D3", "#C71585", "#FF1493"))(7), # 11, 12, 14, 16, 18, 21`
  rep(RColorBrewer::brewer.pal(6, "Greens"), each = 2), # 48, 50 / 52, 54 / 56, 58 / 60, 62 / 64, 66 / 68, 70
  rep(RColorBrewer::brewer.pal(6, "Blues"), each = 2),  # 72, 74 / 76, 78 / 80, 82 / 84, 86 / 88, 90 / 92, 94
  rep(RColorBrewer::brewer.pal(7, "Purples"), each = 2) # 96, 98 / 100, 102 / 104, 106 / 108, 110 / 112, 114 / 116, 118 / 120
)

##Plot the force directed layout  colored by stage - Figure 7D
dpi <- 300
png(file = "endoderm_FDL_stage.nice.png", width = 5 * dpi, height = 5 * dpi)
plotTreeForce(obj.build, "stage.nice", alpha=0.2, alpha.fade=0.08, size=4, density.alpha=T, label.tips=F, view = "View2", colors = stage.colors.new)
dev.off()

plotTreeForce(obj.build, "segment", alpha=0.2, alpha.fade=0.08, size=4, density.alpha=T, label.tips=F, view = "View2")


##Plot differentially expressed TFs for EECs and endo_panc
##Common markers
pond.with.grey <- c("#CECECE", "#CBDAC2", RColorBrewer::brewer.pal(9, "YlGnBu")[3:9])
genes.common <- c("pax6b", "pax4", "nkx2.2a", "neurod1", "insm1b")

##Plot genes on the FDL URD trajectory - Figure 7E
geneList <- c("tnfrsf11a", "atoh1b", "ascl1a", "sox4a", "sox4b", "pou2f3", "sox8b", "best4", "otop2", "cftr", "amn", "cubn", "muc2.1", "agr2")
for (gene in geneList) {
  plotTreeForce(obj.build, gene, alpha=1.2, alpha.fade=0.2, size= 8, density.alpha=T, label.tips=F, view = "View2",
                colors = pond.with.grey)
}

library(rgl)

##Save using the rgl package
rgl.snapshot(filename = "endoderm_FDL_tnfrsf11a.png", fmt = "png")

