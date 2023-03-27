#PART-3
#Determine genes enriched in individual trajectories to particular cell types
#We took each major cell type and its corresponding branch and compared them against each other pairwise to look for differentially expressed genes

#Get the parent segment of each clade to consider as a group
combined.tips <- c("11", "12") 
combined.tips <- c("8", "9", "1")
combined.tips <- c("5", "6", "3", "7")

#Get the cells in that segment and all child segments
cells.combined.tips <- lapply(combined.tips, function(n) whichCells(obj.tree, "segment", value = segChildrenAll(obj.tree, n, include.self = T)))
names(cells.combined.tips) <- combined.tips

#Loop through each of the branches and look for differentially expressed genes

combined.markers <- lapply(combined.tips, function(tip) {
  #Find all other tips
  opposing.tips <- setdiff(combined.tips, tip)
  #Perform pairwise comparisons to each other
  mo <- lapply(opposing.tips, function(other.tip){
    message(paste0(Sys.time(), ":Comparing tip", tip, "to", other.tip, "."))
    #Find differentially expressed genes pairwise
    markers <- markersAUCPR(object = obj.tree, cells.1 = cells.combined.tips[[tip]], cells.2 = cells.combined.tips[[other.tip]], effect.size = 0.4, auc.factor = 1.1)
    #In order to facilitate combining all of the results later, add columns about which two clades were compared and a duplicate entry of each gene that's recovered
    markers$gene <- rownames(markers)
    markers$tip1 <- tip
    markers$tip2 <- other.tip
    return(markers)
  })
  names(mo) <- opposing.tips
  return(mo)
})
names(combined.markers) <- combined.tips
markers.11.12.4 <- combined.markers

#Require that genes are markers against 1 other clades
combined.markers.good <- lapply(combined.markers, function(t) {
  names(which(table(unlist(lapply(t, rownames))) >= 1))
})

#Since genes might be a marker in comparison to several other clades, combine the results into a single table, where each gene is listed only once with the info from the previous comparison where it had the strongest differential expression

combined.markers.best <- lapply(1:length(combined.markers.good), function(i){
  cm <- do.call("rbind", combined.markers[[i]])
  cm <- cm[cm$gene %in% combined.markers.good[[i]], ]
  cmb <- do.call("rbind", lapply(combined.markers.good[[i]], function(g) {
    cmr <- cm[cm$gene == g, ]
    return(cmr[which.max(cmr$AUCPR.ratio), ])
  }))
  
  rownames(cmb) <- cmb$gene
  cmb <- cmb[order(cmb$AUCPR.ratio, decreasing = T), ]
  cmb$exp.global <- apply(obj.tree@logupx.data[rownames(cmb), unlist(obj.tree@tree$cells.in.segment)], 1, mean.of.logs)
  cmb$exp.global.fc <- cmb$nTrans_1 - cmb$exp.global
  return(cmb)
})
names(combined.markers.best) <- combined.tips

#AUCPR along tree to ask for genes that are differential markers of a lineage using URD's tree structure to make a comparison at each branchpoint from a particular cell type up to the root

#Get all of the tips from the tree
tips.in.tree <- as.character(obj.tree@tree$tips)

#Tree segments to use as root for each particular cell population
roots <- rep("11", length(tips.in.tree))
names(roots) <- tips.in.tree
roots[c("5", "6")] <- "8"
roots[c("3", "7")] <- "9"
roots[c("1", "8", "9")] <- "11"

#Perform loop of tests with each tip
markers <- lapply(tips.in.tree, function(t){
  this.root <- roots[t]
  message(paste0(Sys.time(), ":Starting tip ", t, "and root ", this.root))
  these.markers <- aucprTestAlongTree(obj.tree, pseudotime = "pseudotime", tips = as.character(t), genes.use = NULL, must.beat.sibs = 0.6, report.debug = F, root = this.root, auc.factor = 1.1, log.effect.size = 0.4)
  these.markers$gene <- rownames(these.markers)
  these.markers$tip <- t
  return(these.markers)
})
names(markers) <- tips.in.tree

#PART-4: Functions for curating differential expression results
#threshold.tree.markers
#Function to threshold markers from a markersAUCPRAlongTree test with additional criteria
#markers: list of results from markersAUCPRAlongTree tests
#tip: which tip (or element of the list to pursue)
#global.fc: fold.change that gene must have along the trajectory pursued vs. rest of the data
#aucpr.ratio.all: classifier score that gene must exhibit along trajectory test vs. rest of the data
#branch.fc: fold.change that gene must have (in best case) vs. the opposing branch at any branchpoint along the trajectory.
#Returns markers with only a subset of rows retained.

threshold.tree.markers <- function(markers, tip, global.fc = 0.1, branch.fc = 0.4, aucpr.ratio.all = 1.05){
  m <- markers[[tip]]
  #Firstly, lose all genes with global.fc < x
  bye.global.fc <- rownames(m)[m$expfc.all < global.fc]
  #Secondly, lose all genes with branch.fc < x
  bye.branch.fc <- rownames(m)[m$expfc.maxBranch < branch.fc]
  #Thirdly, get rid of stuff essentially worse than random classification on global level
  bye.badglobalaucpr <- rownames(m)[m$AUCPR.ratio.all < aucpr.ratio.all]
  bye.all <- unique(c(bye.global.fc, bye.branch.fc, bye.badglobalaucpr))
  m.return <- m[setdiff(rownames(m), bye.all), ]
  return(m.return)
}

#threshold.clade.markers
#Function to threshold markers of particular clades (see “Combined major branch families”) using
#additional criteria
#markers: result of markersAUCPR
#global.fc: fold.change that gene must have along the trajectory pursued vs. rest of the data
#(during testing, branches were compared pairwise. This compares one branch to all others together.)
#Returns markers with a subset of rows retained

threshold.clade.markers <- function(markers, global.fc = 0.1) {
  m <- markers
  # First off -- lose global FC < x
  bye.globalfc <- rownames(m)[m$exp.global.fc < global.fc]
  m.return <- m[setdiff(rownames(m), bye.globalfc), ]
  return(m.return)
}

#divide.branches
#Function to compare genes between two branches. Use this on a compiled list of markers to do
#a final selection of genes that are specific to one branch or another or markers of both (i.e. when making exocrine pancreas heatmap)
#object: An URD object
#genes: (Character vector) Genes to test
#clust.1: (Character) Cluster 1
#clust.2: (Character) Cluster 2
#clustering: (Character) Clustering to pull from
#exp.fc: (Numeric) Minimum expression fold-change between branches to consider different
#exp.thresh: (Numeric) Minimum fraction of cells in order to consider gene expressed in a branch
#exp.diff: (Numeric) Minimum difference in fraction of cells expressing to consider gene differential
#Returns list of gene names (“specific.1” = specific to clust.1, “specific.2” = specific to clust.2, “markers” = all genes tested)

divide.branches <- function(object, genes, clust.1, clust.2, clustering = "segment", exp.fc = 0.4, exp.thresh = 0.1, exp.diff = 0.1) {
  #Double check which markers are unique to one or the other population
  mcomp <- markersAUCPR(object, clust.1 = clust.1, clust.2 = clust.2, clustering = clustering, effect.size = -Inf, auc.factor = 0, genes.use = genes, frac.min.diff = 0, frac.must.express = 0)
  specific.b <- rownames(mcomp)[abs(mcomp$exp.fc) > exp.fc & mcomp[, 4] < exp.thresh & mcomp[, 5] > pmin((mcomp[, 4] + exp.diff), 1)]
  specific.a <- rownames(mcomp)[abs(mcomp$exp.fc) > exp.fc & mcomp[, 5] < exp.thresh & mcomp[, 4] > pmin((mcomp[, 5] + exp.diff), 1)]
  r <- list(specific.a, specific.b, mcomp)
  names(r) <- c("specific.1", "specific.2", "markers")
  return(r)
}

#Functions for heatmap generation
# 1. Color scale
cols <- (scales::gradient_n_pal(RColorBrewer::brewer.pal(9, "YlOrRd")))(seq(0, 1, length.out = 50))
cols <- c("#CECECE", "#CBDAC2", RColorBrewer::brewer.pal(9, "YlGnBu")[3:9])

# 2. determine timing
#Determines order to plot genes in heatmap. “Expression” is defined as 20% higher expression than the minimum observed value. 
#“Peak” expression is defined as 50% higher expression than minimum observed value.
#The two longest stretches of “peak” expression are found, and then the later one is used. The onset time of the stretch of expression that contains that peak is also determined. 
#Genes are then ordered by the pseudotime at which they enter “peak” expression, leave “peak” expression, start “expression”, and leave “expression”.
#s: result from geneSmoothFit
#genes: genes to order; default is all genes that were fit.
#Returns s but with an additional list entry ($timing) of the order to plot genes

determine.timing <- function(s, genes = rownames(s$mean.expression)) {
  s$timing <- as.data.frame(do.call("rbind", lapply(genes, function(g) {
    sv <- as.numeric(s$scaled.smooth[g, ])
    pt <- as.numeric(colnames(s$scaled.smooth))
    # Figure out baseline expression & threshold for finding peaks
    min.val <- max(min(sv), 0)
    peak.val <- ((1 - min.val)/2) + min.val
    exp.val <- ((1 - min.val)/5) + min.val
    # Run-length encoding of above/below the peak-threshold
    peak.rle <- rle(sv >= peak.val)
    
    peak.rle <- data.frame(lengths = peak.rle$lengths, values = peak.rle$values)
    peak.rle$end <- cumsum(peak.rle$lengths)
    peak.rle$start <- head(c(0, peak.rle$end) + 1, -1)
    # Run-length encoding of above/below the expressed-threshold
    exp.rle <- rle(sv >= exp.val)
    exp.rle <- data.frame(lengths = exp.rle$lengths, values = exp.rle$values)
    exp.rle$end <- cumsum(exp.rle$lengths)
    exp.rle$start <- head(c(0, exp.rle$end) + 1, -1)
    # Take top-two longest peak RLE & select later one. Find stretches that are above peak value
    peak <- which(peak.rle$values)
    # Order by length and take 1 or 2 longest ones
    peak <- peak[order(peak.rle[peak, "lengths"], decreasing = T)][1:min(2, length(peak))]
    # Order by start and take latest one.
    peak <- peak[order(peak.rle[peak, "start"], decreasing = T)][1]
    # Identify the actual peak value within that stretch
    peak <- which.max(sv[peak.rle[peak,"start"] : peak.rle[peak,"end"]]) + peak.rle[peak, "start"] - 1
    # Identify the start and stop of the expressed stretch that contains the peak
    exp.start <- exp.rle[which(exp.rle$end >= peak & exp.rle$start <= peak), "start"]
    exp.end <- exp.rle[which(exp.rle$end >= peak & exp.rle$start <= peak), "end"]
    # Identify values of expression at start and stop
    smooth.start <- sv[exp.start]
    smooth.end <- sv[exp.end]
    # Convert to pseudotime?
    exp.start <- pt[exp.start]
    exp.end <- pt[exp.end]
    peak <- pt[peak]
    # Return a vector
    v <- c(exp.start, peak, exp.end, smooth.start, smooth.end)
    names(v) <- c("pt.start", "pt.peak", "pt.end", "exp.start", "exp.end")
    return(v)
  })))
  rownames(s$timing) <- genes
  
  # Decide on ordering of genes
  s$gene.order <- rownames(s$timing)[order(s$timing$pt.peak, s$timing$pt.start,
                                           s$timing$pt.end, s$timing$exp.end, decreasing = c(F, F, F, T), method = "radix")]
  return(s)
}

#Filter heatmap genes
#Removes undesired (mitochondrial, ribosomal, tandem duplicated genes) from heatmaps for presentation purposes
#Genes (Character vector) genes to check
#Returns genes with undesired genes removed

filter.heatmap.genes <- function(genes) {
  mt.genes <- grep("^mt-", ignore.case = T, genes, value = T)
  many.genes <- grep("\\(1 of many\\)", ignore.case = T, genes, value = T)
  cox.genes <- grep("^cox", ignore.case = T, genes, value = T)
  return(setdiff(genes, c(mt.genes, many.genes, cox.genes)))
}

#PART-5: Generating Heatmaps of gene cascades
#Preparing cascade

#vSMCs_aretry: Seg 9
#get markers using the two different approaches

#Lineage markers from the combined clades
t9 <- threshold.clade.markers(combined.markers.best[["9"]], global.fc = 0.05)
#exocrine pancreas markers from aucprTestAlongTree
m7 <- threshold.tree.markers(markers, "7", global.fc = 0.6)
#liver markers from aucprTestAlongTree
m3 <- threshold.tree.markers(markers, "3", global.fc = 0.6)
exp.markers <- unique(c(rownames(t9), rownames(m7), rownames(m3))) 

# Make a duplicate of the pseudotime measurement (pseudotime.3.7)
obj.tree@pseudotime$pseudotime.3.7 <- obj.tree@pseudotime$pseudotime
# Grab pseudotime of branchpoint
pt.start.3.7 <- as.numeric(obj.tree@tree$segment.pseudotime.limits["3", "start"])
# Figure out lengths (and ratio) of the two branches in pseudotime
pt.end.3.7 <- as.numeric(obj.tree@tree$segment.pseudotime.limits[c("3", "7"), "end"]) -
  pt.start.3.7
pt.ratio.3.7 <- pt.end.3.7[1]/pt.end.3.7[2]
# For cells in the shorter branch (12), subtract the starting pseudotime,
# multiply by the ratio of branch lengths, then add the starting pseudotime back
# in order to stretch the branch.
obj.tree@pseudotime[cellsInCluster(obj.tree, "segment", "7"), "pseudotime.3.7"] <- (obj.tree@pseudotime[cellsInCluster(obj.tree, "segment", "7"), "pseudotime.3.7"] - pt.start.3.7) * pt.ratio.3.7 + pt.start.3.7

#Calculate spline curves using segments 9, 3, and 7. 
spline.7 <- geneSmoothFit(obj.tree, pseudotime = "pseudotime.3.7", cells = cellsInCluster(obj.tree, "segment", c("11", "9", "7")), genes = exp.markers, method = "spline", moving.window = 5, cells.per.window = 8, pseudotime.per.window = 0.005, spar = 0.5, verbose = T)

#Finding dynamic genes who change their expression values considerably
#Which genes change their actual mean expression value by at least 0.75
change.real <- apply(spline.7$mean.smooth, 1, function(x) diff(range(x)))
genes.mean <- names(which(change.real >= 0.75))

spline.7 <- geneSmoothFit(obj.tree, pseudotime = "pseudotime", cells = cellsInCluster(obj.tree, "segment", c("11", "7", "9")), genes = genes.mean, method = "spline", moving.window = 5, cells.per.window = 8, pseudotime.per.window = 0.005, spar = 0.5, verbose = T)

spline.3 <- geneSmoothFit(obj.tree, pseudotime = "pseudotime", cells = cellsInCluster(obj.tree, "segment", c("3")), genes = setdiff(exp.markers, c("cnn1a", "fbp2")), method = "spline", moving.window = 5, cells.per.window = 8, pseudotime.per.window = 0.005, spar = 0.5, verbose = T)

#Finding dynamic genes who change their expression values considerably
#Which genes change their actual mean expression value by at least 0.75
change.real <- apply(spline.7$mean.smooth, 1, function(x) diff(range(x)))
genes.mean.7 <- names(which(change.real >= 0.5))

spline.3 <- geneSmoothFit(obj.tree, pseudotime = "pseudotime.3.7", cells = cellsInCluster(obj.tree, "segment", c("11", "3", "9")), genes = genes.mean, method = "spline", moving.window = 5, cells.per.window = 8, pseudotime.per.window = 0.005, spar = 0.5, verbose = T)

spline.3.7 <- geneSmoothFit(obj.tree, pseudotime = "pseudotime.3.7", cells = cellsInCluster(obj.tree, "segment", c("11", "9", "3", "7")), genes = exp.markers, method = "spline", moving.window = 5, cells.per.window = 8, pseudotime.per.window = 0.005, spar = 0.5, verbose = T)

#Finding dynamic genes who change their expression values considerably
#Which genes change their actual mean expression value by at least 0.75
change.real <- apply(spline.3.7$mean.smooth, 1, function(x) diff(range(x)))
genes.mean.3.7 <- names(which(change.real >= 0.75))

spline.3.7 <- geneSmoothFit(obj.tree, pseudotime = "pseudotime.3.7", cells = cellsInCluster(obj.tree, "segment", c("11", "9", "3", "7")), genes =  genes.mean, method = "spline", moving.window = 5, cells.per.window = 8, pseudotime.per.window = 0.005, spar = 0.5, verbose = T)

#Want to plot a heatmap that shows expression in hepatopancreatic progenitors and then each branch (i.e. liver and exocrine pancreas) as separate columns. Hence we need to crop each spline fit to the correct pseudotime range and then combine them into a single one that can be plotted as a 3-column heatmap

pt.3v7 <- obj.tree@tree$segment.pseudotime.limits["7", "start"] #pseudotime where the crop needs to happen
splines.lp <- cropSmoothFit(spline.7, pt.min = pt.3v7, pt.max = Inf)
names(splines.lp) <- c("vSMC-artery_1")
splines.lp.hm <- combineSmoothFit(splines.lp) #Combine into single one

#Calculate gene expression timing for ordering rows
spline.3.7 <- determine.timing(s = spline.3.7)
spline.7 <- determine.timing(s = spline.7)
spline.3 <- determine.timing(s = spline.3)

d3v7 <- divide.branches(obj.tree, genes.mean, clust.1 = "3", clust.2 = "7", exp.fc = 0.4,
                        exp.thresh = 0.2, exp.diff = 0.1)
order.3.7 <- filter.heatmap.genes(setdiff(spline.3.7$gene.order, c(d3v7$specific.1, d3v7$specific.2)))
order.3 <- filter.heatmap.genes(intersect(spline.3$gene.order, d3v7$specific.1))
order.7 <- filter.heatmap.genes(intersect(spline.7$gene.order, d3v7$specific.2))
gene.order <- c(order.3.7, order.3, order.7)

#Output the gene table
base.path <- ("~/Box/zfext/annotations_celltype_curated_newMama/mural_cells/obj_cascades/")
table.save <- data.frame(gene = order.7, stringsAsFactors = F)
table.save$clade.AUCPR.ratio <- t9[table.save$gene, "AUCPR.ratio"]
table.save$clade.exp.fc <- t9[table.save$gene, "exp.fc"]
table.save$clade.exp.fc.global <- t9[table.save$gene, "exp.global.fc"]
table.save$vSMC.AUCPR.ratio.all <- m3[table.save$gene, "AUCPR.ratio.all"]
table.save$vSMC.AUCPR.ratio.maxBranch <- m3[table.save$gene, "AUCPR.ratio.maxBranch"]
table.save$vSMC.exp.fc.all <- m3[table.save$gene, "expfc.all"]
table.save$vSMC.exp.fc.best <- m3[table.save$gene, "expfc.maxBranch"]
write.csv(table.save, quote = F, file = paste0(base.path, "vSMCs_artery.csv"))

# Make sure any values <0 in the spline curves get set to 0 so that the heatmap
# scale doesn't get messed up.
splines.lp.hm$scaled.smooth[splines.lp.hm$scaled.smooth < 0] <- 0

# Determine where to place column separators (i.e. how many columns will each
# cell type occupy in the heatmap )
colsep <- cumsum(as.numeric(head(unlist(lapply(splines.lp, function(x) ncol(x$scaled.smooth))),
                                 -1)))
# Determine where to place row separators (i.e. how many common markers, and
# markers are specific to each cell type)
rowsep <- cumsum(c(length(order.3.7), length(order.3)))
# Open a PDF and generate the heatmap 
pdf(paste0('vSMC_artery_cascade_heatmap.pdf'), width=12, height=16)
gplots::heatmap.2(x = as.matrix(splines.lp.hm$scaled.smooth[gene.order, ]), Rowv = F, Colv = F,
                  dendrogram = "none", col = cols, trace = "none", density.info = "none", key = F,
                  cexCol = 0.8, cexRow = 0.35, margins = c(8, 8), lwid = c(0.3, 4), lhei = c(0.5,
                                                                                             6), labCol = NA, colsep = colsep, rowsep = rowsep, sepwidth = c(0.1, 0.2), labRow = genes.to.plot)
title(main = "vSMC_artery")
title(main = "common_prog", line = -71, adj = 0)
title(main = "vSMC_art_1", line = -71, adj = 0.45)
title(main = "vSMC_art_2", line = -71, adj = 0.96)

dev.off()


##Plotting some markers on the heatmap

genes.to.plot <- c("acta2", "tagln", "lmod1b", "pbx3b", "lmo4b", "meis1b", "ndufa4l2a", "pdgfra", "pdgfrb", "foxf2a", "isl2b", "klf2a")
rownames.to.plot <- gene.order
rtp <- rownames.to.plot %in% genes.to.plot
rownames.to.plot[!rtp] <- ""
rownames.to.plot[rtp] <- paste0("-", rownames.to.plot[rtp])
# Open a PDF and generate the heatmap pdf(paste0(base.path,
# '/heatmaps/retina-rgc-mainfig.pdf'), width=6, height=10)
pdf(paste0('vSMC_artery_cascade_heatmap_labeled.pdf'), width=12, height=16)
gplots::heatmap.2(as.matrix(splines.lp.hm$scaled.smooth[gene.order, ]), Rowv = F, Colv = F,
                  dendrogram = "none", col = cols, trace = "none", density.info = "none", key = F,
                    cexCol = 0.9, cexRow = 1.8, margins = c(8, 10), lwid = c(0.3, 4), lhei = c(0.3,
                                                                                             4), labCol = NA, labRow = rownames.to.plot)
title(main = "vSMC_artery")
title(main = "common_prog", line = -71, adj = 0)
title(main = "vSMC_art_1", line = -71, adj = 0.45)
title(main = "vSMC_art_2", line = -71, adj = 0.96)

dev.off()

##Now look at the spline curves to understand when and how these gene expression changes
pericyte.genes <- c("ndufa4l2a", "pdgfra", "pdgfrb", "cntfr", "etv4")
smc.genes <- c("acta2", "tagln", "lmod1b", "myl9a", "ndufa4l2a")
plotSmoothFit(smoothed.fit = splines.lp.hm, genes = genes.to.plot, scaled = T, multiplot = F, alpha.data = 0.2, alpha.smooth = 1)
plotSmoothFitMultiCascade(smoothed.fit = splines.lp, genes = smc.genes, scaled = T, alpha.data = 0.2, alpha.smooth = 1)

##Save the splines
saveRDS(splines.lp, "~/Box/zfext/annotations_celltype_curated_newMama/mural_cells/obj_cascades/vSMC_artery_mesoderm_spline.rds")
saveRDS(splines.lp.hm, "~/Box/zfext/annotations_celltype_curated_newMama/mural_cells/obj_cascades/vSMC_artery_mesoderm_combined_heatmap.rds")

splines.lp <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/mural/obj_cascades/vSMC_artery_mesoderm_spline.rds")
splines.lp.hm <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/mural/obj_cascades/vSMC_artery_mesoderm_combined_heatmap.rds")

gridExtra::grid.arrange(grobs = lapply(c("acta2", "ndufa4l2a", "lmod1b", "pdgfrb"), plotTree,
                                       object = obj.tree, label.x = F, plot.cells = T), ncol = 2)

##There appears to be a few acta2+/ndufa4l2a+ cells in the vSMC_artery_1 branch

##Now let's assess the branchpoint between the two visceral SMC populations
#Preparing cascade

#vSMCs_aretry: Seg 9
#get markers using the two different approaches

#Lineage markers from the combined clades
t8 <- threshold.clade.markers(combined.markers.best[["8"]], global.fc = 0.05)
#exocrine pancreas markers from aucprTestAlongTree
m5 <- threshold.tree.markers(markers, "5", global.fc = 0.6)
#liver markers from aucprTestAlongTree
m6 <- threshold.tree.markers(markers, "6", global.fc = 0.6)
exp.markers <- unique(c(rownames(t8), rownames(m5), rownames(m6))) 

# Make a duplicate of the pseudotime measurement (pseudotime.5.6)
obj.tree@pseudotime$pseudotime.5.6 <- obj.tree@pseudotime$pseudotime
# Grab pseudotime of branchpoint
pt.start.5.6 <- as.numeric(obj.tree@tree$segment.pseudotime.limits["5", "start"])
# Figure out lengths (and ratio) of the two branches in pseudotime
pt.end.5.6 <- as.numeric(obj.tree@tree$segment.pseudotime.limits[c("5", "6"), "end"]) -
  pt.start.5.6
pt.ratio.5.6 <- pt.end.5.6[1]/pt.end.5.6[2]
# For cells in the shorter branch (12), subtract the starting pseudotime,
# multiply by the ratio of branch lengths, then add the starting pseudotime back
# in order to stretch the branch.
cells.seg.6 <- cellsInCluster(obj.tree, "segment", "6")
cells.seg.5 <- cellsInCluster(obj.tree, "segment", "5")
cells.seg.8 <- cellsInCluster(obj.tree, "segment", "8")
seg.8.pt <- obj.tree@pseudotime[cells.seg.8, ]
cells.seg.8.add <- rownames(seg.8.pt)[which(seg.8.pt$pseudotime >=0.55)]
cells.seg.6.total <- unlist(unique(list(cells.seg.6, cells.seg.8.add)))
cells.seg.5.total <- unlist(unique(list(cells.seg.5, cells.seg.8.add)))

obj.tree@pseudotime[cells.seg.6.total, "pseudotime.5.6"] <- (obj.tree@pseudotime[cells.seg.6.total, "pseudotime.5.6"] - pt.start.5.6) * pt.ratio.5.6 + pt.start.5.6
obj.tree@pseudotime[cells.seg.5.total, "pseudotime.5.6"] <- (obj.tree@pseudotime[cells.seg.5.total, "pseudotime.5.6"] - pt.start.5.6) * pt.ratio.5.6 + pt.start.5.6

#Calculate spline curves using segments 9, 3, and 7. 
spline.6 <- geneSmoothFit(obj.tree, pseudotime = "pseudotime.5.6", cells = cellsInCluster(obj.tree, "segment", c("11", "6", "8")), genes = exp.markers, method = "spline", moving.window = 5, cells.per.window = 8, pseudotime.per.window = 0.005, spar = 0.5, verbose = T)
spline.5 <- geneSmoothFit(obj.tree, pseudotime = "pseudotime.5.6", cells = cellsInCluster(obj.tree, "segment", c("11", "5", "8")), genes = exp.markers, method = "spline", moving.window = 5, cells.per.window = 8, pseudotime.per.window = 0.005, spar = 0.5, verbose = T)

#Finding dynamic genes who change their expression values considerably
#Which genes change their actual mean expression value by at least 0.75
change.real <- apply(spline.5$mean.smooth, 1, function(x) diff(range(x)))
genes.mean <- names(which(change.real >= 0.75))

spline.6 <- geneSmoothFit(obj.tree, pseudotime = "pseudotime.5.6", cells = cellsInCluster(obj.tree, "segment", c("11", "6", "8")), genes = genes.mean, method = "spline", moving.window = 5, cells.per.window = 8, pseudotime.per.window = 0.005, spar = 0.5, verbose = T)

spline.5 <- geneSmoothFit(obj.tree, pseudotime = "pseudotime.5.6", cells = cellsInCluster(obj.tree, "segment", c("11", "5", "8")), genes = genes.mean, method = "spline", moving.window = 5, cells.per.window = 8, pseudotime.per.window = 0.005, spar = 0.5, verbose = T)

spline.5.6 <- geneSmoothFit(obj.tree, pseudotime = "pseudotime.5.6", cells = cellsInCluster(obj.tree, "segment", c("11", "8", "5", "6")), genes = exp.markers, method = "spline", moving.window = 5, cells.per.window = 8, pseudotime.per.window = 0.005, spar = 0.5, verbose = T)

#Want to plot a heatmap that shows expression in intestinal SMC progenitors and then each branch (i.e. circular and longitudinal SMC) as separate columns. Hence we need to crop each spline fit to the correct pseudotime range and then combine them into a single one that can be plotted as a 3-column heatmap

pt.5v6 <- 0.57 #pseudotime where the crop needs to happen
pt.8.end <- obj.tree@tree$segment.pseudotime.limits["8", "end"]
splines.lp <- list(cropSmoothFit(spline.5.6, pt.min = -Inf, pt.max = pt.8.end), cropSmoothFit(spline.5, pt.min = pt.5v6, pt.max = Inf), cropSmoothFit(spline.6, pt.min = pt.5v6, pt.max = Inf))
names(splines.lp) <- c("viSMC_prog","int_SMCs_long", "int_SMCs_circ")
splines.lp.hm <- combineSmoothFit(splines.lp) #Combine into single one

#Calculate gene expression timing for ordering rows
spline.5.6 <- determine.timing(s = spline.5.6)
spline.5 <- determine.timing(s = spline.5)
spline.6 <- determine.timing(s = spline.6)

d5v6 <- divide.branches(obj.tree, genes.mean, clust.1 = "5", clust.2 = "6", exp.fc = 0.4,
                        exp.thresh = 0.4, exp.diff = 0.1)
order.5.6 <- filter.heatmap.genes(setdiff(spline.5.6$gene.order, c(d5v6$specific.1, d5v6$specific.2)))
order.5 <- filter.heatmap.genes(intersect(spline.5$gene.order, d5v6$specific.1))
order.6 <- filter.heatmap.genes(intersect(spline.6$gene.order, d5v6$specific.2))
gene.order <- c(order.5.6, order.5, order.6)

#Output the gene table
base.path <- ("~/Box/zfext/annotations_celltype_curated_newMama/mural/obj_cascades/")
table.save <- data.frame(gene = order.5, stringsAsFactors = F)
table.save$clade.AUCPR.ratio <- t8[table.save$gene, "AUCPR.ratio"]
table.save$clade.exp.fc <- t8[table.save$gene, "exp.fc"]
table.save$clade.exp.fc.global <- t8[table.save$gene, "exp.global.fc"]
table.save$viSMC.AUCPR.ratio.all <- m5[table.save$gene, "AUCPR.ratio.all"]
table.save$viSMC.AUCPR.ratio.maxBranch <- m5[table.save$gene, "AUCPR.ratio.maxBranch"]
table.save$viSMC.exp.fc.all <- m5[table.save$gene, "expfc.all"]
table.save$viSMC.exp.fc.best <- m5[table.save$gene, "expfc.maxBranch"]
write.csv(table.save, quote = F, file = paste0(base.path, "viSMCs_intestine_circ_vs_long.csv"))

# Make sure any values <0 in the spline curves get set to 0 so that the heatmap
# scale doesn't get messed up.
splines.lp.hm$scaled.smooth[splines.lp.hm$scaled.smooth < 0] <- 0

# Determine where to place column separators (i.e. how many columns will each
# cell type occupy in the heatmap )
colsep <- cumsum(as.numeric(head(unlist(lapply(splines.lp, function(x) ncol(x$scaled.smooth))),
                                 -1)))
# Determine where to place row separators (i.e. how many common markers, and
# markers are specific to each cell type)
rowsep <- cumsum(c(length(order.5.6), length(order.5)))
# Open a PDF and generate the heatmap 
pdf(paste0('viSMC_intestine_vs_sb_cascade_heatmap.pdf'), width=12, height=16)
gplots::heatmap.2(x = as.matrix(splines.lp.hm$scaled.smooth[gene.order, ]), Rowv = F, Colv = F,
                  dendrogram = "none", col = cols, trace = "none", density.info = "none", key = F,
                  cexCol = 0.8, cexRow = 0.35, margins = c(8, 8), lwid = c(0.3, 4), lhei = c(0.5,
                                                                                             6), labCol = NA, colsep = colsep, rowsep = rowsep, sepwidth = c(0.1, 0.2))
title(main = "visceral_SMC")
title(main = "common_prog", line = -71, adj = 0)
title(main = "sb_SMCs", line = -71, adj = 0.65)
title(main = "int_SMCs", line = -71, adj = 0.96)

dev.off()


##Plotting some markers on the heatmap

genes.to.plot <- c("acta2", "tagln", "lmod1b", "myocd", "smtna", "smtnb", "nid2a", "foxq1b", "foxf2a", "desmb", "plxna3", "cabp1a")
rownames.to.plot <- gene.order
rtp <- rownames.to.plot %in% genes.to.plot
rownames.to.plot[!rtp] <- ""
rownames.to.plot[rtp] <- paste0("-", rownames.to.plot[rtp])
# Open a PDF and generate the heatmap pdf(paste0(base.path,
# '/heatmaps/retina-rgc-mainfig.pdf'), width=6, height=10)
pdf(paste0('viSMC_intestine_sb_cascade_heatmap_labeled.pdf'), width=12, height=16)
gplots::heatmap.2(as.matrix(splines.lp.hm$scaled.smooth[gene.order, ]), Rowv = F, Colv = F,
                  dendrogram = "none", col = cols, trace = "none", density.info = "none", key = F,
                  cexCol = 0.9, cexRow = 1.8, margins = c(8, 10), lwid = c(0.3, 4), lhei = c(0.3,
                                                                                             4), labCol = NA, labRow = rownames.to.plot)
title(main = "visceral_SMC")
title(main = "common_prog", line = -71, adj = 0)
title(main = "sb_SMCs", line = -71, adj = 0.65)
title(main = "int_SMCs", line = -71, adj = 0.96)
dev.off()

gridExtra::grid.arrange(grobs = lapply(c("acta2", "ndufa4l2a", "bgna", "tdrd9"), plotTree,
                                       object = urd.build, label.x = F, plot.cells = F), ncol = 2)

genes.to.plot <- c("smtnb", "smtna", "kcnk18", "desmb", "il13ra2", "rgs2", "sfrp1b", "foxf2a", "gucy1a1")
genes.to.plot <- c("foxq1a", "tcf21", "cremb", "itpr1a", "tead3")

##Now look at the spline curves to understand when and how these gene expression changes
plotSmoothFit(smoothed.fit = splines.lp.hm, genes = genes.to.plot, scaled = T, multiplot = F, alpha.data = 0.2, alpha.smooth = 1)
plotSmoothFitMultiCascade(smoothed.fit = splines.lp, genes = genes.to.plot, scaled = T, alpha.data = 0.2, alpha.smooth = 1, colors = c("#CC00CC", "#08519c", "#8B4513"))

saveRDS(splines.lp, "~/Box/zfext/annotations_celltype_curated_newMama/mural_cells/obj_cascades/viSMC_intestine_sb_spline.rds")
saveRDS(splines.lp.hm, "~/Box/zfext/annotations_celltype_curated_newMama/mural/obj_cascades/viSMC_intestine_sb_combined_heatmap_pt_adj.rds")

splines.lp <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/mural/obj_cascades/viSMC_intestine_sb_spline.rds")
splines.lp.hm <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/mural/obj_cascades/viSMC_intestine_sb_combined_heatmap_pt_adj.rds")

##Specify the path to save the markers
save.path <- "~/Box/zfext/results/06b-URD/cascades/mural"

segs.to.name <- c("5", "6", "3", "7", "1")
new.seg.names <- c("int_SMCs-circular", "int_SMCs-long", "vSMC_artery_2", "vSMC_artery_1", "SMC-visceral?")
obj.tree <- nameSegments(obj.tree, segments = segs.to.name, segment.names = new.seg.names)
plotTree(obj.tree, "stage.group", label.segments = T, cell.size = 1.2)

#Determine tips to run DE for
tips.to.run <- as.character(obj.tree@tree$segment.names)

genes.use <- NULL #Calculate for all genes

#Calculate the markers for each of the populations
gene.markers <- list()
for(tipn in 1:length(tips.to.run)) {
  tip <- tips.to.run[[tipn]]
  print(paste0(Sys.time(), ":", tip))
  markers <- aucprTestAlongTree(obj.tree, pseudotime = "pseudotime", tips = tip, log.effect.size = 0.4, auc.factor = 0.4, max.auc.threshold = 0.85, frac.must.express = 0.1, frac.min.diff = 0, genes.use = genes.use, root = "11", only.return.global = F, must.beat.sibs = 0.6, report.debug = T)
  saveRDS(markers, paste0(save.path, "/aucpr/", tip, ".rds"))
  gene.markers[[tip]] <- markers
}

saveRDS(gene.markers, paste0(save.path, "/aucpr/", "mesoderm_gene_markers_allTips.rds"))
gene.markers <- readRDS(paste0(save.path, "/aucpr_markers/", "mesoderm_gene_markers_allTips.rds"))

# Separate actual marker lists from the stats lists
gene.markers.de <- lapply(gene.markers, function(x) x[[1]])
gene.markers.stats <- lapply(gene.markers[1:5], function(x) x[[2]])
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

casc <- geneCascadeProcess(object = obj.tree, pseudotime='pseudotime', cells = cellsAlongLineage(obj.tree, c("9", "7"), remove.root=F), 
                           genes = rownames(gene.markers.de$vSMC_artery_1), 
                           moving.window=5, cells.per.window=18, 
                           pseudotime.per.window = 0.01, verbose = T, verbose.genes = T)

saveRDS(casc, "~/Box/zfext/annotations_celltype_curated_newMama/mural_cells/obj_cascades/vSMC-artery_pericyte-3_cascade.rds")

spline.artery <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/mural/obj_cascades/vSMC_artery_mesoderm_spline.rds")

##Just do cascade between vSMC-artery_1 and pericyte-3
#Prepare cascade for vSMC + pericyte-3 (Seg 3)
t9 <- threshold.clade.markers(combined.markers.best[["9"]], global.fc = 0.05)
#Progenitor markers from AUCPRTestalongTree
m7 <- threshold.tree.markers(markers, "7", global.fc = 0.6) # vSMC markers from aucprTestAlongTree
exp.markers <- unique(c(rownames(t9), rownames(m7))) 

#Calculate spline curves using segments 3 and 9
spline.7 <- geneSmoothFit(obj.tree, pseudotime = "pseudotime", cells = cellsInCluster(obj.tree,
                                                                                      "segment", c("11", "9", "7")), genes = exp.markers, method = "spline", moving.window = 5,
                          cells.per.window = 25, pseudotime.per.window = 0.005, spar = 0.5, verbose = T)

#Calculate gene expression timing for ordering rows
#spline.7 <- determine.timing(s = spline.7)
order.7 <- filter.heatmap.genes(spline.7$gene.order) 

# Output gene table
table.save <- data.frame(gene = order.7, stringsAsFactors = F)
table.save$smc.AUCPR.ratio.all <- m7[table.save$gene, "AUCPR.ratio.all"]
table.save$smc.AUCPR.ratio.maxBranch <- m7[table.save$gene, "AUCPR.ratio.maxBranch"]
table.save$smc.exp.fc.all <- m7[table.save$gene, "expfc.all"]
table.save$smc.exp.fc.best <- m7[table.save$gene, "expfc.maxBranch"]
write.csv(table.save, quote = F, file = paste0(base.path, "vSMC_pericytes_transition.csv"))

# Make sure any values <0 in the spline curves get set to 0 so that the heatmap
# scale doesn't get messed up.
spline.7$scaled.smooth[spline.7$scaled.smooth < 0] <- 0
# Open a PDF and generate the heatmap pdf(paste0(base.path,
# '/heatmaps/retina-smc.pdf'), width=6, height=10)
pdf("~/Box/zfext/annotations_celltype_curated_newMama/mural/obj_cascades/vSMC_pericyte_transition_2.pdf", height = 32, width = 26)
gplots::heatmap.2(x = as.matrix(spline.7$scaled.smooth[order.7, ]), Rowv = F, Colv = F,
                  dendrogram = "none", col = cols, trace = "none", density.info = "none", key = F,
                  cexCol = 0.8, cexRow = 0.15, margins = c(8, 8), lwid = c(0.3, 4), lhei = c(0.3,
                                                                                             4), labCol = NA)
title(main = "vSMC_pericyte-3_continuum")

dev.off()


##plot some spline curves
##Now look at the spline curves to understand when and how these gene expression changes
smc.genes.to.plot <- c("acta2", "tagln", "lmod1b", "myl9a", "mylkb", "tnfaip6", "barx1", "rcn3", "alcama", "fn1a")
peri.genes.plot  <- c("ndufa4l2a", "klf2a", "lmo4b", "pgk1", "igf1rb", "plk2a")
peri.genes.late <- c("ndufa4l2a", "vwde", "bgna", "dcn", "col2a1b", "loxl5b", "ldha", "pdgfra", "pdgfrb")
genes.plot <- c("cxcl12a", "acta2", "tagln", "myl9a", "lmod1b", "klf2a", "ndufa4l2a", "bgna", "dcn", )
genes.plot.early <- c("fgf10", "cxcl12a", "rbpms2a", "cdh2", "snrpf", "foxp4")
genes.plot.mes <- c("foxq1a", "isl2b", "barx1", "rcn3", "foxf2a")
plotSmoothFit(smoothed.fit = spline.7, genes = genes.plot.early, scaled = T, multiplot = F, alpha.data = 0.4, alpha.smooth = 1.2, lwd.smooth = 1.4)
plotSmoothFit(smoothed.fit = spline.7, genes = smc.genes.to.plot, scaled = T, multiplot = F, alpha.data = 0.4, alpha.smooth = 1.2, lwd.smooth = 1.4)
plotSmoothFit(smoothed.fit = spline.7, genes = peri.genes.late, scaled = T, multiplot = F, alpha.data = 0.4, alpha.smooth = 1.2, lwd.smooth = 1.4)
plotSmoothFit(smoothed.fit = spline.7, genes = genes.plot.mes, scaled = T, multiplot = F, alpha.data = 0.4, alpha.smooth = 1.2, lwd.smooth = 1.4)
plotSmoothFitMultiCascade(smoothed.fit = spline.7, genes = genes.plot, scaled = T, alpha.data = 0.2, alpha.smooth = 1)

saveRDS(spline.7, "~/Box/zfext/annotations_celltype_curated_newMama/mural/obj_cascades/vSMC_pericyte_transition_spline.rds")
spline.7 <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/mural/obj_cascades/vSMC_pericyte_transition_spline.rds")
spline <- spline.3

## Let's look at some of the spline curves of the post-LRE genes
#Calculate smoothed spline fits of their expression
spline <- geneSmoothFit(obj.tree, method = "spline", pseudotime = "pseudotime", cells = cellsAlongLineage(obj.tree, "3"), genes = rownames(gene.markers.de$vSMC_artery_1), moving.window = 2, cells.per.window = 5, spar = 0.9)

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
only.spline <- cropSmoothFit(spline, pt.min = pt.crop)

splines <- list(only.spline)
names(splines) <- "vSMC-artery_1_pericyte-3"

#Plotting gene expression
genes.plot <- c("acta2", "ndufa4l2a", "pdgfra", "pdgfrb")

plotSmoothFit(smoothed.fit = spline, genes = genes.plot, scaled = T, multiplot = F, alpha.data = 0.2, alpha.smooth = 1)
plotSmoothFitMultiCascade(smoothed.fit = splines, genes = genes.plot, scaled = T, alpha.data = 0.2, alpha.smooth = 1)

