##This script is to generate a URD tree of pericytes and vascular SMCs from the transient/quiescent states identified in our cell atlas

library(Seurat)
library(URD)

obj1 <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/mural/obj_seurat/mural_seurat.rds")
markers.mura <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/mural/obj_seurat/mural_markers.rds")

DimPlot(obj1, label = T)
DimPlot(obj1, group.by = "stage.group")

obj.crop <- subset(obj1, idents = c("7", "22", "11", "20", "14", "5", "6"))

#Trial-1: URD trajectory using seuratToURD2
seuratToURD2 <- function(seurat.object) {
  if (requireNamespace("Seurat", quietly = TRUE)) {
    # Create an empty URD object
    ds <- new("URD")
    
    # Copy over data
    ds@logupx.data <- as(as.matrix(seurat.object@assays$RNA@data), "dgCMatrix")
    if(!any(dim(seurat.object@assays$RNA@counts) == 0)) ds@count.data <- as(as.matrix(seurat.object@assays$RNA@counts[rownames(seurat.object@assays$RNA@data), colnames(seurat.object@assays$RNA@data)]), "dgCMatrix")
    
    # Copy over metadata
    ## TO DO - grab info
    get.data <- NULL
    if (.hasSlot(seurat.object, "data.info")) { 
      get.data <- as.data.frame(seurat.object@assays$RNA@data.info)
    } else if (.hasSlot(seurat.object, "meta.data")) { 
      get.data <- as.data.frame(seurat.object@meta.data) 
    }
    if(!is.null(get.data)) {
      di <- colnames(get.data)
      m <- grep("res|cluster|Res|Cluster", di, value=T, invert = T) # Put as metadata if it's not the result of a clustering.
      discrete <- apply(get.data, 2, function(x) length(unique(x)) / length(x))
      gi <- di[which(discrete <= 0.015)]
      ds@meta <- get.data[,m,drop=F]
      ds@group.ids <- get.data[,gi,drop=F]
    }
    
    # Copy over var.genes
    if(length(seurat.object@assays$RNA@var.features > 0)) ds@var.genes <- seurat.object@assays$RNA@var.features
    
    # Move over tSNE projection
    if (.hasSlot(seurat.object, "tsne.rot")) {
      if(!any(dim(seurat.object@tsne.rot) == 0)) {
        ds@tsne.y <- as.data.frame(seurat.object@tsne.rot)
        colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
      }
    } else if (.hasSlot(seurat.object, "reductions")) {
      if(("tsne" %in% names(seurat.object@reductions)) && !any(dim(seurat.object@reductions$tsne) == 0)) {
        ds@tsne.y <- as.data.frame(seurat.object@reductions$tsne@cell.embeddings)
        colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
      }
    }
    
    # Move over PCA results
    if (.hasSlot(seurat.object, "pca.x")) {
      if(!any(dim(seurat.object@pca.x) == 0)) {
        ds@pca.load <- seurat.object@pca.x
        ds@pca.scores <- seurat.object@pca.rot
        warning("Need to set which PCs are significant in @pca.sig")
      }
      ## TO DO: Convert SVD to sdev
    } else if (.hasSlot(seurat.object, "reductions")) {
      if(("pca" %in% names(seurat.object@reductions)) && !any(dim(Loadings(seurat.object, reduction = "pca")) == 0)) {
        ds@pca.load <- as.data.frame(Loadings(seurat.object, reduction = "pca"))
        ds@pca.scores <- as.data.frame(seurat.object@reductions$pca@cell.embeddings)
        ds@pca.sdev <- seurat.object@reductions$pca@stdev
        ds@pca.sig <- pcaMarchenkoPastur(M=dim(ds@pca.scores)[1], N=dim(ds@pca.load)[1], pca.sdev=ds@pca.sdev)
      }
    }
    return(ds)
  } else {
    stop("Package Seurat is required for this function. To install: install.packages('Seurat')\n")
  }
}

obj <- seuratToURD2(obj.crop)

counts <- obj.crop@assays$RNA@counts
metadata <- obj.crop@meta.data
obj <- createURD(count.data = counts, meta = metadata, min.cells = 3, min.genes = 200, verbose = T)

obj <- obj.urd

# Get variable genes for each group of stages
stages <- sort(unique(obj@meta$stage.group))
cells.each.stage <- lapply(stages, function(stage) rownames(obj@meta)[which(obj@meta$stage.group == stage)])
var.genes.by.stage <- lapply(1:length(stages), function(n) findVariableGenes(obj, cells.fit = cells.each.stage[[n]], set.object.var.genes = F,
                                                                             diffCV.cutoff = 0.4, mean.min = 0.005, mean.max = 100, main.use = stages[[n]], do.plot = T))

#Take union of variable genes across all stages
var.genes <- sort(unique(unlist(var.genes.by.stage)))

#Set variable genes in object
obj@var.genes <- var.genes

obj <- calcPCA(obj)
vg.dups <- duplicated(as.data.frame(as.matrix(t(obj@logupx.data[obj@var.genes,
]))))
if (length(which(vg.dups)) > 0) {
  print(paste("Removing", length(which(vg.dups)), "cell(s) with duplicated variable gene expression."))
  not.dup.cells <- colnames(obj@logupx.data)[!vg.dups]
  obj <- urdSubset(obj, not.dup.cells)
}
obj <- calcTsne(obj, perplexity = 30, theta = 0.5)
obj <- graphClustering(obj, dim.use = "pca", num.nn = c(10, 15, 20, 30), do.jaccard = T, method = "Louvain")

plotDim(obj, "stage.group", legend = T, plot.title = "Developmental Stage", alpha = 2.4)

umap.coord <- obj1@reductions$umap@cell.embeddings
head(umap.coord)
head(obj@tsne.y)
colnames(umap.coord) <- c("tSNE1", "tSNE2")
class(umap.coord)
class(obj@tsne.y)

obj@tsne.y <- as.data.frame(umap.coord)
plotDim(obj, "tfa")
plotDim(obj, "stage.group", label.clusters = T)
plotTree(obj, label.segments = T)

cells.to.exclude <- whichCells(obj, "Louvain-15", "17")
cells.to.include <- setdiff(WhichCells(obj.seurat), cells.to.exclude)
obj <- urdSubset(obj, cells.keep = cells.to.include)

# Save object
message(paste0(Sys.time(), ": Saving"))
saveRDS(obj, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD.rds")
obj <- readRDS(file="~/Box Sync/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD.rds")

obj <- calcKNN(obj)
outliers <- knnOutliers(obj, nn.1 = 1, nn.2 = 10, x.max = 25, slope.r = 1.1, int.r = 4, slope.b = 1.85, int.b = 8, title = "Identifying outliers by k-NN distance")                                                                             
gridExtra::grid.arrange(grobs=list(#Plot some apoptotic markers
  plotDim(obj, "isg15", alpha = 0.4, point.size = 0.5), 
  plotDim(obj, "foxo3b", alpha = 0.4, point.size = 0.5),
  plotDim(obj, "gadd45aa", alpha = 0.4, point.size = 0.5),
  #Figure out which clusters correspond to these cells
  plotDimHighlight(obj, clustering = "Louvain-15", cluster = "17", legend = F)))

apoptotic.like.cells <- cellsInCluster(obj, "Louvain-15", c(24, 32, 4, 1))
plotDim(obj, "Louvain-15", legend = T, plot.title = "Louvain_Jaccard Graph-based Clustering (15 NNs)", alpha = 1)
plotDim(obj, "tfa")
markersAUCPR(obj, "1", clustering = "Louvain-15", exp.thresh = 1.0, frac.must.express = 0.6)
plotDim(obj, "neurod1")
plotDimHighlight(obj, clustering = "Louvain-15", cluster = "29", legend = F)

#Subset object to eliminate outliers
cells.keep <- setdiff(colnames(obj@logupx.data), c(outliers))
obj <- urdSubset(obj, cells.keep = cells.keep)

cells.exclude <- WhichCells(obj.seurat, idents = c("1", "40"))
cells.all <- setdiff(WhichCells(obj.seurat), cells.exclude)
obj <- urdSubset(obj.all, cells.all)

saveRDS(obj, "~/Box Sync/zfext/results/06b-URD/obj_subsets/endoderm_superset/sculptedKNN/endoderm_URD_trimmed_v4.rds")
obj <- readRDS("~/Box Sync/zfext/results/06b-URD/obj_subsets/endoderm_superset/endoderm_URD_trimmed_v4.rds")
obj <- calcDM(obj, knn = 100, sigma.use = 12)

saveRDS(obj, file="~/Box Sync/zfext/results/06b-URD/obj_subsets/endoderm_superset/sculptedKNN/endoderm_URD_withDM_v4.rds")
obj <- readRDS(file="~/Box Sync/zfext/results/06b-URD/obj_subsets/endoderm_superset/sculptedKNN/endoderm_URD_withDM_v4.rds")
obj <- readRDS(file="~/Box Sync/zfext/results/06a-URD/2021-03 KNN Modification/endoderm_KNNsculpted_dm9.3.rds")

stage.colors <- c("darkblue", "#FFCCCC", "#99CC00", "#33CC00", "cyan3",
                  "gold", "goldenrod", "darkorange", "indianred1", "indianred2", "plum", "deepskyblue2",
                  "lightgrey")
plotDimArray(obj, reduction.use = "dm", dims.to.plot = 1:18, outer.title = "Diffusion Map (Sigma 14, 100 NNs): Stage", label="stage.group", alpha = 0.4, plot.title="", legend=T, discrete.colors = stage.colors)
plotDim(obj, "stage.group", transitions.plot = 10000, plot.title="Developmental stage (with transitions)")

saveRDS(obj, "~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD_DM.rds")

#Calculate pseudotime
message(paste0(Sys.time(), ": Defining earliest stage as the root.cells"))
root.cells.1 <- rownames(obj@meta)[obj@meta$stage.group == " 36-46"]
root.cells.2 <- rownames(obj@meta)[obj@meta$stage.group == " 48-58"]
root.cells <- unlist(unique(list(root.cells.1, root.cells.2)))
plotDimHighlight(obj, "stage.group", " 36-46", plot.title = "Root is 36 hpf cells")

#Do the flood
message(paste0(Sys.time(), ": Do the pseudotime_flood"))
flood.result <- floodPseudotime(obj, root.cells = root.cells, n = 100, minimum.cells.flooded = 2, verbose = T)

#Save the result
message(paste0(Sys.time(), ": Saving the flood object"))
saveRDS(flood.result, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD_flood.rds")
flood.result <- readRDS(file="~/Box Sync/zfext/results/06a-URD/obj_subsets/endoderm_superset/seuratToURD_endoderm/endoderm_ds_flood_sculptedKNN.rds")

#Process pseudotime floods
message(paste0(Sys.time(), ": Processing pseudotime_flood"))
obj <- floodPseudotimeProcess(obj, flood.result, floods.name = "pseudotime")

pseudotimePlotStabilityOverall(obj)

#Inspect pseudotime in Diffusion Plots
plotDimArray(object = obj, reduction.use = "dm", dims.to.plot = 1:14,
             label = "pseudotime", plot.title = "", outer.title = "Diffusion Map labeled by pseudotime",
             legend = F, alpha = 0.4)

#Inspect pseudotime using tSNE plot
message(paste0(Sys.time(), ": Plotting Pseudotime on tSNE plot"))
pdf(file=paste0(plot.path, sample, "_pseudotime.pdf"), width = 8, height = 8)
plotDim(obj, "pseudotime", plot.title = "Pseudotime")

#Plotting the distribution of pseudotime from each stage group
#Create a dataframe that contains both pseudotime and stage information
gg.data <- cbind(obj@pseudotime, obj@meta[rownames(obj@pseudotime), ])
#Plot
ggplot(gg.data, aes(x=pseudotime, color = stage.group, fill = stage.group)) + geom_density(alpha = 0.4)

#Plotting distances
message(paste0(Sys.time(), ": Plotting Distances"))
pdf(file=paste0(plot.path, sample, "_pseudotime.pdf"), width = 8, height = 8)
plotDists(obj, "pseudotime", "stage.group", plot.title="Pseudotime by stage")

#Save DM and PT information
saveRDS(obj, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_withDMandPT.rds")
obj <- readRDS(file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_withDMandPT.rds")


#Part 3 - Determining Tips
message(paste0(Sys.time(), ": Cropping the cells from the final stage_group"))
cells_120h <- grep("120", colnames(obj@logupx.data), value = T)
cells_120h <- colnames(obj@logupx.data)[obj@meta$stage.group == "120"]
cells_108h <- colnames(obj@logupx.data)[obj@meta$stage.group == " 108-118"]
cells.keep <- unlist(unique(list(cells_120h, cells_108h)))

obj_120h <- urdSubset(obj, cells.keep = cells.keep)

#Perform PCA/tSNE on final stage cells
message(paste0(Sys.time(), ": Load the variable genes specific to this stage"))
var.genes.120h <- scan("120_var.txt")
obj_120h@var.genes <- obj@var.genes

#Calculate PCA
message(paste0(Sys.time(), ": Calculating PCA"))
obj_120h <- calcPCA(obj_120h)

#Calculate tSNE
message(paste0(Sys.time(), ": Calculating tSNE"))
set.seed(18)
obj_120h <- calcTsne(obj_120h, perplexity = 20, theta = 0.5)

obj_120h <- graphClustering(obj_120h, num.nn = 50, do.jaccard = T, method = "Louvain")
obj_120h <- graphClustering(obj_120h, num.nn = c(5, 8, 10, 15, 30, 40), do.jaccard = T, method = "Louvain")
obj_120h <- graphClustering(obj_120h, num.nn = c(10, 15, 20, 30, 40, 50), do.jaccard = T, method = "Infomap")

clusterings <- c(paste0("Infomap-", c(10, 15, 20, 30, 40, 50)), paste0("Infomap-", c(10, 15, 20, 30, 40, 50)))
clusterings <- c(paste0("Louvain-", c(50, 40, 30, 20, 15, 10)))

for (c in clusterings) {
  plot(plotDim(obj_120h, c, legend = T, label.clusters = T))
}

#Looking at batch information
pdf(file=paste0(plot.path, sample, "_tSNE_batch.pdf"), width = 8, height = 8)
plotDim(obj_120h, "Louvain-15", plot.title = "Louvain-15_graph", label.clusters = T)
plotDim(obj_120h, "stage.group", plot.title = "Louvain-15_graph")
plotDimHighlight(obj_120h, clustering = "Louvain-15", cluster = "11", legend = F)

cells.bil <- whichCells(obj, "Louvain-15", "22")
clusters <- sort(unique(obj_120h@group.ids$`Louvain-15`))
pr.markers <- lapply(clusters, function(c) markersAUCPR(obj_120h, clust.1 = c, clustering = "Louvain-15", genes.use = obj_120h@var.genes))
names(pr.markers) <- clusters

plotDim(obj_120h, "kcnk18", plot.title="FOXA2 (Pharyngeal endoderm marker)")

#Make a set of data.frames to keep track during cluster assignment
I30.n <- length(unique(obj_120h@group.ids$`Louvain-15`))
I30.cluster.assignments <- data.frame(cluster = 1:I30.n, name = rep(NA, I30.n, name = rep(NA, I30.n), tip = rep(NA, I30.n)), row.names = 1:I30.n)

I20.n <- length(unique(obj_120h@group.ids$`Louvain-15`))
I20.cluster.assignments <- data.frame(cluster = 1:I20.n, name = rep(NA, I20.n, name = rep(NA, I20.n), tip = rep(NA, I20.n)), row.names = 1:I20.n)

plotDot(obj_120h, genes = c("ndufa4l2a", "nr4a1", "atp1a1b", "mmrn2a", "igfbp7", "atp1b4", "prrx1b", "kcnj8", "kcne4", "abcc9", "wfdc2", "fbln5", "pitx2", "lmx1bb",
                            "loxa", "gdf10a", "tgfb2", "ptgdsb.2", "olfm1b", "ccl20b", "myl9a", "acta2", "tagln", "ifitm1", "tlx1", "lama5", "tinagl1", "spns2", "rgs5b", "trdn", "rcn3", "tnfrsf9a", "golim4b", "rflnb"), clustering = "Louvain-15")

#Markers for early exocrine pancreas
plotDot(obj_120h, genes = c("c9", "pcsk6", "pdx1", "surf4", "scgn", "prox1a", "isl1", "muc5.1", "muc5.2", "muc5.3", "agr2", "cldnh", "malb", "chs1", "slc6a8", "muc2.1", "nkx2.1", "foxe1", "elovl2",
                            "bmp4", "anxa13", "faxdc2", "slc7a9", "chia.2", "ace2", "gcga", "lgals2a", "lyz", "sdc4", "cftr", "anxa3a", "anxa4", "anxa2b", "dck",
                            "srgn", "msna", "cd63"), clustering = "Louvain-15")

saveRDS(obj, file=paste0(obj.path, sample, "_withDMandPT.rds"))

I30.cluster.assignments["1", "name"] <- "stellate-cells" #Use as tip
I30.cluster.assignments["2", "name"] <- "progenitors" #Use as tips
I30.cluster.assignments["3", "name"] <- "vSMC-artery_1" #Use as tip
I30.cluster.assignments["4", "name"] <- "progenitors" #Use as tip
I30.cluster.assignments["5", "name"] <- "int_SMCs" #Use as tip
I30.cluster.assignments["6", "name"] <- "sb_SMCs" #Use as tips
I30.cluster.assignments["7", "name"] <- "vSMC_artery_2" #Use as tip
#I30.cluster.assignments["8", "name"] <- "fibroblasts" #Use as tip

#Generate final clusterings
message(paste0(Sys.time(), ": Combine clustering assignments from two clusterings"))
I30.cluster.assignments$clustering <- "Louvain-15"
I20.cluster.assignments$clustering <- "Louvain-15"
cluster.assignments <- rbind(I30.cluster.assignments)

#Remove any clusters that weren't assigned an identity
message(paste0(Sys.time(), ": Removing clusters without an identity"))
cluster.assignments <- cluster.assignments[!is.na(cluster.assignments$name), ]

#Renumber clusters
cluster.assignments$cluster.new <- 1:nrow(cluster.assignments)

#Create blank clusterings for obj_120h
obj_120h@group.ids$clusters.120h.name <- NA
obj_120h@group.ids$clusters.120h.num <- NA

#Copy cell identities over for each cluster
for (i in 1:nrow(cluster.assignments)) {
  cells <- cellsInCluster(obj_120h, clustering = cluster.assignments[i, "clustering"],
                          cluster = cluster.assignments[i, "cluster"])
  obj_120h@group.ids[cells, "clusters.120h.name"] <- cluster.assignments[i, "name"]
  obj_120h@group.ids[cells, "clusters.120h.num"] <- as.character(cluster.assignments[i, "cluster.new"])
}

#Transfer clusterings to main object
obj@group.ids$`tip.clusters` <- NA
obj@group.ids[rownames(obj_120h@group.ids), "tip.clusters"] <- obj_120h@group.ids$`Louvain-15`

obj@group.ids$`Louvain-15` <- NA
obj@group.ids[rownames(obj_120h@group.ids), "Louvain-15"] <- obj_120h@group.ids$`Louvain-15`

obj@group.ids$`Cluster` <- NA
obj@group.ids[rownames(obj_120h@group.ids), "Cluster"] <- obj_120h@group.ids$clusters.120h.name

obj@group.ids$`Cluster-Num` <- NA
obj@group.ids[rownames(obj_120h@group.ids), "Cluster-Num"] <- obj_120h@group.ids$clusters.120h.num

#Save objects
message(paste0(Sys.time(), ": Saving the 120h_seurat object"))
saveRDS(obj_120h, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_120h_Louvain-15.rds")
obj_120h <- readRDS(file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_120h_Louvain-15.rds")

message(paste0(Sys.time(), ": Saving the full object with 120h clustering added"))
saveRDS(obj, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD_Louvain-15.rds")
obj <- readRDS(file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD_Louvain-15.rds")

message(paste0(Sys.time(), ": Saving the data.frame with tips"))
write.csv(cluster.assignments, file = "~/Box Sync/zfext/results/06b-URD/obj_subsets/endoderm_superset/sculptedKNN/endoderm_tips-use_Louvain-15_scKNN_v4.csv")

#Plot tips in diffusion map
obj@group.ids$pop <- NA
obj@group.ids[cellsInCluster(obj, "Cluster", "Liver"), "pop"] <- "6"
plotDim(obj, label = "pop", plot.title = "Pancreas DCs 13 vs 14", reduction.use = "dm", dim.x = 4, dim.y = 5, 
        legend = F, alpha = 0.35, colors = stage.colors)

#PART-4
#Biased random walks

message(paste0(Sys.time(), ": Load previous saved object"))
obj <- readRDS(file="~/Box Sync/zfext/results/06b-URD/obj_subsets/endoderm_superset/seuratToURD/endoderm_URD_Louvain-15_scKNN_2.rds")

obj@group.ids[rownames(obj_120h@group.ids), "tip.clusters"] <- obj_120h@group.ids$`Louvain-15`

#Define parameters of logistic function to bias transition probablities
diffusion.logistic <- pseudotimeDetermineLogistic(obj, "pseudotime", optimal.cells.forward = 20, max.cells.back = 40, pseudotime.direction = "<",
                                                  do.plot = T, print.values = T)

#Create biased transition matrix
message(paste0(Sys.time(), ": Creating biased transition matrix"))
biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(obj, pseudotime = "pseudotime",
                                                        logistic.params = diffusion.logistic, pseudotime.direction = "<"))

#Define the root cells
message(paste0(Sys.time(), ": Defining root-cells"))
root.cells <- rownames(obj@meta)[obj@meta$stage.group == " 24-34"]
plotDimHighlight(obj, "stage.group", " 24-34")

#Define tip cells
message(paste0(Sys.time(), ": Defining tip-cells"))
tips <- setdiff(unique(obj@group.ids[, clustering]), NA)
this.tip <- tips[tip.to.walk]

#Simulate the biased random walks from each tip
message(paste0(Sys.time(), ": Simulating random walks from each tip"))
tip.walks <- simulateRandomWalksFromTips(obj, tip.group.id = "tip.clusters", root.cells = root.cells, transition.matrix = biased.tm, n.per.tip = 25000, root.visits = 1,
                                         max.steps = 5000, verbose = T)

biased.tm.good <- intersect(rownames(biased.tm), tip.cells)

walks <- lapply(rownames(table(obj@group.ids$`Cluster-Num`)), function(c) {
  # Exclude any tip cells that for whatever reason didn't end up in the
  # biased TM (e.g. maybe not assigned a pseudotime).
  tip.cells <- rownames(obj@group.ids)[which(obj@group.ids$tip.clusters == c)]
  tip.cells.good <- intersect(tip.cells, rownames(biased.tm))
  # Perform the random walk simulation
  this.walk <- simulateRandomWalk(start.cells = tip.cells.good, transition.matrix = biased.tm,
                                  end.cells = root.cells, n = 50000, end.visits = 1, verbose.freq = 1000,
                                  max.steps = 5000)
  return(this.walk)
})
names(tip.walks) <- rownames(table(obj@group.ids$`Cluster-Num`))

saveRDS(tip.walks, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_TipWalks_Louvain-15.rds")
tip.walks <- readRDS(file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_TipWalks_Louvain-15.rds")

#Process the biased random walks into visitation frequencies
message(paste0(Sys.time(), ": Processing the biased random walks"))
obj <- processRandomWalksFromTips(obj, tip.walks, verbose = T)

#Visualize visitation of cells from each tip
plotDim(obj, "tip.clusters", plot.title = "Cells in each tip")

plotDim(obj, "visitfreq.log.2", plot.title = "Visitation frequency from tip-1 (log10)", transitions.plot = 10000)
plotDim(obj, "visitfreq.log.1", plot.title = "Visitation frequency from tip-2 (log10)", transitions.plot = 10000)

saveRDS(obj, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_withWalks_Louvain-15.rds")
obj <- readRDS(file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_withWalks_Louvain-15.rds")

tips.walked <- setdiff(unique(obj@group.ids$`Cluster-Num`), NA)

#PART 5 - Building URD Tree
library(URD)
library(rgl)

#Set up knitr to capture rgl output
rgl::setupKnitr()

obj <- readRDS(file="~/Box Sync/zfext/results/06b-URD/obj_subsets/endoderm_superset/seuratToURD/endoderm_withWalks_Louvain-15_scKNN_3.rds")

#Load tip cells
message(paste0(Sys.time(), ": Loading tip cells"))
obj <- loadTipCells(obj, tips = "tip.clusters")

#Combine a few sets of tips where we walked from 2 groups of cells and are intermixed in the diffusion map
#Liver cells
obj <- combineTipVisitation(obj, "2", "3", "3")
#Exocrine Pancreas
obj <- combineTipVisitation(obj, "4", "5", "5")
#Anterior Intestine
obj <- combineTipVisitation(obj, "3", "15", "15")
#Posterior Intestine
obj <- combineTipVisitation(obj, "1", "9", "9")
#Biliary cells
obj <- combineTipVisitation(obj, "11", "19", "11")
#Digestive Gut
#obj <- combineTipVisitation(obj, "5", "10", "5")
#obj <- combineTipVisitation(obj, "30", "31", "30")


#Build the actual tree
message(paste0(Sys.time(), ": Decide on tips to use for tree construction"))
tips.to.exclude <- c("2", "4")
tips.to.use <- setdiff(as.character(1:7), tips.to.exclude)

#Build the tree
obj.built <- buildTree(object = obj, pseudotime = "pseudotime", divergence.method = "ks",
                       tips.use = 1:3, weighted.fusion = T, use.only.original.tips = T,
                       cells.per.pseudotime.bin = 25, bins.per.pseudotime.window = 8, minimum.visits = 1,
                       visit.threshold = 0.7, p.thresh = 0.1, save.breakpoint.plots = NULL, dendro.node.size = 100,
                       min.cells.per.segment = 10, min.pseudotime.per.segment = 0.0001, verbose = F)

obj.tree <- buildTree(obj, pseudotime = "pseudotime", tips.use = tips.to.use, divergence.method = "preference", cells.per.pseudotime.bin = 10, 
                      bins.per.pseudotime.window = 8, save.all.breakpoint.info = T, p.thresh = 0.001, verbose = F)

#Name the tips
tip.names <- unique(obj@group.ids[, c("Cluster", "Cluster-Num")])
tip.names <- tip.names[complete.cases(tip.names), ]
obj.tree <- nameSegments(obj.tree, segments = tip.names$`Cluster-Num`, segment.names = tip.names$`Cluster`)
obj.built <- nameSegments(obj.built, segments = tip.names$`Cluster-Num`, segment.names = tip.names$`Cluster`)

plotTree(obj.built, label.segments = T, "stage.group")
plotTree(obj.tree, label.segments = T, "stage.group")

saveRDS(obj.tree, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/vSMCs_artery_TREE.rds")
obj.tree <- readRDS(file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/vSMCs_artery_TREE.rds")

png("mural_cell_mesoderm_TREE_colored.png", width = 5*dpi, height = 5*dpi)
plotTree(obj.tree, "stage.nice", label.segments = F, plot.cells = T, legend = F, cell.alpha = 0.8, cell.size = 1.8, label.x = F, continuous.colors = stage.colors.new)
dev.off()

genes.plot <- c("hhex", "hmox1a", "glula")
for (gene in genes.plot) {
  plot(plotTree(obj.tree, gene))
}

gridExtra::grid.arrange(grobs = lapply(c("abcc9", "mcamb", "ednraa", "cd248a"), plotTree,
                                       object = obj.tree, label.x = F, plot.cells = T), ncol = 2)

gridExtra::grid.arrange(grobs = lapply(c("pdx1", "onecut1"), FeaturePlot,
                                       object = obj1), ncol = 2)

gridExtra::grid.arrange(grobs = lapply(c("acta2", "ndufa4l2a", "bgna", "pdgfrb"), plotTree,
                                       object = obj.tree, label.x = F, plot.cells = F, cell.alpha = 0.8, cell.size = 1.0), ncol = 2)

gridExtra::grid.arrange(grobs = lapply(common_genes.tfs, plotTree,
                                       object = obj.tree, label.x = F, plot.cells = T, cell.alpha = 0.8, cell.size = 1), ncol = 6)

#PART-2
#Create branchpoint plots to understand cellular dynamics between individual cell types
np.layout <- branchpointPreferenceLayout(obj.tree, pseudotime = "pseudotime", lineages.1 = "3", lineages.2 = "7", parent.of.lineages = c("9") , opposite.parent = 1, min.visit = 2)
plotBranchpoint(obj.tree, np.layout, label = "segment", point.alpha = 0.8, populations = c("aSMC-1", "aSMC-2"),
                pt.lim = c(0.7, 0.3), xlab = "", ylab = "", legend = T, axis.lines = F,
                fade.low = 0, title = "Stage.group")

#Examine gene expression in the branchpoint points
genes.plot <- c("acta2", "macrod2", "bgna", "ndufa4l2a", "emilin2a", "pthlha")
genes.plot <- c("amy2a", "ela3l", "ela2l", "cpa4", "cpa5", "prss59.1")
genes.plot <- c("gpx3", "gpx4a", "egr1")
branch.plots <- lapply(genes.plot, function(gene) plotBranchpoint(obj.tree, np.layout, label = gene, point.alpha = 1, populations = c("aSMC-1", "aSMC-2"),
                                                                  pt.lim = c(0.7, 0.3), xlab = "", ylab = "", title = gene, legend = F,
                                                                  axis.lines = F, fade.low = 0.66))
gridExtra::grid.arrange(grobs = branch.plots, ncol = 3)

#Manual refinement
message(paste0(Sys.time(), ": Descriptive names to be used on the dendrogram"))
new.seg.names <- c("int_SMCs", "sb_SMCs", "vSMC-artery_1", "vSMC-artery_2", "stellate-cells")
segs.to.name <- c("5", "6", "3", "7", "1")
obj.tree <- nameSegments(obj.tree, segments = segs.to.name, segment.names = new.seg.names)


##Ok, just want to try to find what these transition states are more likely to form - for that purpose, I am going to use the transition probablity calculations of URD

##Use the entire object to create an URD object
counts <- obj.smc@assays$RNA@counts
metadata <- obj.smc@meta.data
obj <- createURD(count.data = counts, meta = metadata, min.cells = 3, min.genes = 200, verbose = T)

# Get variable genes for each group of stages
stages <- sort(unique(obj@meta$stage.group))
cells.each.stage <- lapply(stages, function(stage) rownames(obj@meta)[which(obj@meta$stage.group == stage)])
var.genes.by.stage <- lapply(1:length(stages), function(n) findVariableGenes(obj, cells.fit = cells.each.stage[[n]], set.object.var.genes = F,
                                                                             diffCV.cutoff = 0.4, mean.min = 0.005, mean.max = 100, main.use = stages[[n]], do.plot = T))

#Take union of variable genes across all stages
var.genes <- sort(unique(unlist(var.genes.by.stage)))

#Set variable genes in object
obj@var.genes <- var.genes

obj <- calcPCA(obj)
vg.dups <- duplicated(as.data.frame(as.matrix(t(obj@logupx.data[obj@var.genes,
]))))
if (length(which(vg.dups)) > 0) {
  print(paste("Removing", length(which(vg.dups)), "cell(s) with duplicated variable gene expression."))
  not.dup.cells <- colnames(obj@logupx.data)[!vg.dups]
  obj <- urdSubset(obj, not.dup.cells)
}
obj <- calcTsne(obj, perplexity = 30, theta = 0.5)
obj <- graphClustering(obj, dim.use = "pca", num.nn = c(10, 15, 20, 30), do.jaccard = T, method = "Louvain")

plotDim(obj, "stage.group", legend = T, plot.title = "Developmental Stage", alpha = 2.4)

umap.coord <- obj.smc@reductions$umap@cell.embeddings
head(umap.coord)
head(obj@tsne.y)
colnames(umap.coord) <- c("tSNE1", "tSNE2")
class(umap.coord)
class(obj@tsne.y)

obj@tsne.y <- as.data.frame(umap.coord)
plotDim(obj, "ndufa4l2a")
plotDim(obj, "stage.group", label.clusters = T)
plotTree(obj, label.segments = T)

cells.to.exclude <- whichCells(obj, "Louvain-15", "17")
cells.to.include <- setdiff(WhichCells(obj.seurat), cells.to.exclude)
obj <- urdSubset(obj, cells.keep = cells.to.include)

# Save object
message(paste0(Sys.time(), ": Saving"))
saveRDS(obj, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_all_URD.rds")
obj <- readRDS(file="~/Box Sync/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD.rds")

obj <- calcKNN(obj)
outliers <- knnOutliers(obj, nn.1 = 1, nn.2 = 10, x.max = 75, slope.r = 1.2, int.r = 4, slope.b = 1.2, int.b = 8, title = "Identifying outliers by k-NN distance")                                                                             
gridExtra::grid.arrange(grobs=list(#Plot some apoptotic markers
  plotDim(obj, "isg15", alpha = 0.4, point.size = 0.5), 
  plotDim(obj, "foxo3b", alpha = 0.4, point.size = 0.5),
  plotDim(obj, "gadd45aa", alpha = 0.4, point.size = 0.5),
  #Figure out which clusters correspond to these cells
  plotDimHighlight(obj, clustering = "Louvain-15", cluster = "17", legend = F)))

apoptotic.like.cells <- cellsInCluster(obj, "Louvain-15", c(24, 32, 4, 1))
plotDim(obj, "Louvain-15", legend = T, plot.title = "Louvain_Jaccard Graph-based Clustering (15 NNs)", alpha = 1)
plotDim(obj, "tfa")
markersAUCPR(obj, "1", clustering = "Louvain-15", exp.thresh = 1.0, frac.must.express = 0.6)
plotDim(obj, "neurod1")
plotDimHighlight(obj, clustering = "Louvain-15", cluster = "29", legend = F)

#Subset object to eliminate outliers
cells.keep <- setdiff(colnames(obj@logupx.data), c(outliers))
obj <- urdSubset(obj, cells.keep = cells.keep)

cells.exclude <- WhichCells(obj.seurat, idents = c("1", "40"))
cells.all <- setdiff(WhichCells(obj.seurat), cells.exclude)
obj <- urdSubset(obj.all, cells.all)

saveRDS(obj, "~/Box Sync/zfext/results/06b-URD/obj_subsets/endoderm_superset/sculptedKNN/endoderm_URD_trimmed.rds")
obj <- readRDS("~/Box Sync/zfext/results/06b-URD/obj_subsets/endoderm_superset/endoderm_URD_trimmed_v4.rds")
obj <- calcDM(obj, knn = 100, sigma.use = NULL)

saveRDS(obj, file="~/Box Sync/zfext/results/06b-URD/obj_subsets/endoderm_superset/sculptedKNN/endoderm_URD_withDM_v4.rds")
obj <- readRDS(file="~/Box Sync/zfext/results/06b-URD/obj_subsets/endoderm_superset/sculptedKNN/endoderm_URD_withDM_v4.rds")
obj <- readRDS(file="~/Box Sync/zfext/results/06a-URD/2021-03 KNN Modification/endoderm_KNNsculpted_dm9.3.rds")

stage.colors <- c("darkblue", "#FFCCCC", "#99CC00", "#33CC00", "cyan3",
                  "gold", "goldenrod", "darkorange", "indianred1", "indianred2", "plum", "deepskyblue2",
                  "lightgrey")
plotDimArray(obj, reduction.use = "dm", dims.to.plot = 1:18, outer.title = "Diffusion Map (Sigma 14, 100 NNs): Stage", label="stage.group", alpha = 0.4, plot.title="", legend=T, discrete.colors = stage.colors)
plotDim(obj, "stage.group", transitions.plot = 1000, plot.title="Developmental stage (with transitions)")

saveRDS(obj, "~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD_all_DM.rds")

#Calculate pseudotime
message(paste0(Sys.time(), ": Defining earliest stage as the root.cells"))
root.cells.1 <- rownames(obj@meta)[obj@meta$stage.group == " 36-46"]
root.cells.2 <- rownames(obj@meta)[obj@meta$stage.group == " 48-58"]
root.cells <- unlist(unique(list(root.cells.1, root.cells.2)))
plotDimHighlight(obj, "stage.group", " 36-46", plot.title = "Root is 24 hpf cells")

#Do the flood
message(paste0(Sys.time(), ": Do the pseudotime_flood"))
flood.result <- floodPseudotime(obj, root.cells = root.cells, n = 100, minimum.cells.flooded = 2, verbose = T)

#Save the result
message(paste0(Sys.time(), ": Saving the flood object"))
saveRDS(flood.result, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD_flood.rds")
flood.result <- readRDS(file="~/Box Sync/zfext/results/06a-URD/obj_subsets/endoderm_superset/seuratToURD_endoderm/endoderm_ds_flood_sculptedKNN.rds")

#Process pseudotime floods
message(paste0(Sys.time(), ": Processing pseudotime_flood"))
obj <- floodPseudotimeProcess(obj, flood.result, floods.name = "pseudotime")

pseudotimePlotStabilityOverall(obj)

#Inspect pseudotime in Diffusion Plots
plotDimArray(object = obj, reduction.use = "dm", dims.to.plot = 1:14,
             label = "pseudotime", plot.title = "", outer.title = "Diffusion Map labeled by pseudotime",
             legend = F, alpha = 0.4)

#Inspect pseudotime using tSNE plot
message(paste0(Sys.time(), ": Plotting Pseudotime on tSNE plot"))
pdf(file=paste0(plot.path, sample, "_pseudotime.pdf"), width = 8, height = 8)
plotDim(obj, "pseudotime", plot.title = "Pseudotime")

#Plotting the distribution of pseudotime from each stage group
#Create a dataframe that contains both pseudotime and stage information
gg.data <- cbind(obj@pseudotime, obj@meta[rownames(obj@pseudotime), ])
#Plot
ggplot(gg.data, aes(x=pseudotime, color = stage.group, fill = stage.group)) + geom_density(alpha = 0.4)

#Plotting distances
message(paste0(Sys.time(), ": Plotting Distances"))
pdf(file=paste0(plot.path, sample, "_pseudotime.pdf"), width = 8, height = 8)
plotDists(obj, "pseudotime", "stage.group", plot.title="Pseudotime by stage")

#Differential expression with precision-recall along URD dendrogram
#Genes are considered differentially expressed if they are expressed in atleast 10% of cells in the trajectory segment under consideration and their mean expression was upregulated 1.5X compared to the sibling and the gene was 1.25x better than a random classifier for the population as defined by the area under a precision-recall curve

#Determine tips to run DE for
tips.to.run <- as.character(obj.tree@tree$segment.names)
genes.use <- NULL #Calculate for all genes

#Specify paths
save.path <- "~/Box/zfext/results/06b-URD/cascades/mural-cells"
plot.path <- "~/Box/zfext/results/06b-URD/cascades/mural-cells/spline_heatmaps/"

#Calculate the markers for each of the populations
gene.markers <- list()
for(tipn in 1:length(tips.to.run)) {
  tip <- tips.to.run[[tipn]]
  print(paste0(Sys.time(), ":", tip))
  markers <- aucprTestAlongTree(obj.tree, pseudotime = "pseudotime", tips = tip, log.effect.size = 0.4, auc.factor = 0.4, max.auc.threshold = 0.85, frac.must.express = 0.1, frac.min.diff = 0, genes.use = genes.use, root = "46", only.return.global = F, must.beat.sibs = 0.6, report.debug = T)
  saveRDS(markers, paste0(save.path, "/aucpr_markers/", tip, ".rds"))
  gene.markers[[tip]] <- markers
}

saveRDS(gene.markers, paste0(save.path, "/aucpr_markers/", "gene_markers_allTips.rds"))
gene.markers <- readRDS(paste0(save.path, "/aucpr_markers/", "gene_markers_allTips.rds"))

# Separate actual marker lists from the stats lists
gene.markers.de <- lapply(gene.markers, function(x) x[[1]])
gene.markers.stats <- lapply(gene.markers[1:18], function(x) x[[2]])
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
tips.to.run <- setdiff(tips.to.run, c("stellate-cells"))
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
saveRDS(gene.cascades, file = paste0(save.path, "/obj_cascades/", "mural_all.genes.cascades.rds"))

DimPlot(obj1, cells.highlight = cellsAlongLineage(obj.tree, "3"))


frac.exp <- rowSums(as.matrix(obj.tree@logupx.data > 0))/ncol(obj.tree@logupx.data)
#Consider all genes expressed in 1% of definitive erythropoietic cells
smc.genes <- names(frac.exp)[which(frac.exp > 0.01)]

## Let's look at some of the spline curves of the post-LRE genes
#Calculate smoothed spline fits of their expression
spline <- geneSmoothFit(obj.tree, method = "spline", pseudotime = "pseudotime", cells = cellsAlongLineage(obj.tree, "7"), genes = smc.genes, moving.window = 2, cells.per.window = 5, spar = 0.9)

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

only.spline.art2 <- only.spline
splines <- list(only.spline)
names(splines) <- c("vSMC-artery_1")

#Plotting gene expression
genes.plot <- c("fbp2", "cnn1a", "acta2", "ldha", "ndufa4l2a")

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
genes.plot <- rownames(spline$scaled.expression)

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
pdf(paste0(plot.path, "vSMC-artery_1_heatmap_spline.pdf"), width = 17, height = 22)
# Plot the heatmap
gplots::heatmap.2(as.matrix(ss[genes.plot[h.ss$order], ]), Rowv = F, RowSideColors = h.ss.clust.col[h.ss$order],
                  Colv = F, dendrogram = "none", col = cols, trace = "none", density.info = "none",
                  key = F, cexCol = 0.6, cexRow = 0.3, margins = c(8, 8), lwid = c(0.3, 4), lhei = c(0.3,
                                                                                                     4), labCol = NA)
dataset <- "vSMC-artery_1"  
# Put a title on it
title(dataset, line = -1, adj = 0.48, cex.main = 4)

dev.off()

common_genes <- intersect(rownames(gene.markers.de$`vSMC-artery_1`), rownames(gene.markers.de$vSMC_artery_2))
common_genes.tfs <- common_genes[which(common_genes %in% tf.list$Symbol)]

genes.artery1 <- setdiff(rownames(gene.markers.de$`vSMC-artery_1`), rownames(gene.markers.de$vSMC_artery_2))
genes.artery2 <- setdiff(rownames(gene.markers.de$`vSMC-artery_2`), rownames(gene.markers.de$vSMC_artery_1))






##Now try to resolve a trajectory for the other NCC-derived pericyte/SMC cluster
DimPlot(obj.smc, label = T)
obj.crop <- subset(obj.smc, idents = c("9", "8", "26", "27", "21", "13"))

obj <- seuratToURD2(obj.crop)

obj <- calcTsne(obj, perplexity = 30, theta = 0.5)
obj <- graphClustering(obj, dim.use = "pca", num.nn = c(10, 15, 20, 30), do.jaccard = T, method = "Louvain")

plotDim(obj, "stage.group", legend = T, plot.title = "Developmental Stage", alpha = 2.4)

umap.coord <- obj1@reductions$umap@cell.embeddings
head(umap.coord)
head(obj@tsne.y)
colnames(umap.coord) <- c("tSNE1", "tSNE2")
class(umap.coord)
class(obj@tsne.y)

obj@tsne.y <- as.data.frame(umap.coord)
plotDim(obj, "tfa")
plotDim(obj, "stage.group", label.clusters = T)
plotTree(obj, label.segments = T)

cells.to.exclude <- whichCells(obj, "Louvain-15", "17")
cells.to.include <- setdiff(WhichCells(obj.seurat), cells.to.exclude)
obj <- urdSubset(obj, cells.keep = cells.to.include)

# Save object
message(paste0(Sys.time(), ": Saving"))
saveRDS(obj, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD_NCC_derived.rds")
obj <- readRDS(file="~/Box Sync/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD.rds")

obj <- calcKNN(obj)
outliers <- knnOutliers(obj, nn.1 = 1, nn.2 = 10, x.max = 26, slope.r = 1.1, int.r = 4, slope.b = 1.85, int.b = 8, title = "Identifying outliers by k-NN distance")                                                                             
gridExtra::grid.arrange(grobs=list(#Plot some apoptotic markers
  plotDim(obj, "isg15", alpha = 0.4, point.size = 0.5), 
  plotDim(obj, "foxo3b", alpha = 0.4, point.size = 0.5),
  plotDim(obj, "gadd45aa", alpha = 0.4, point.size = 0.5),
  #Figure out which clusters correspond to these cells
  plotDimHighlight(obj, clustering = "Louvain-15", cluster = "17", legend = F)))

apoptotic.like.cells <- cellsInCluster(obj, "Louvain-15", c(24, 32, 4, 1))
plotDim(obj, "Louvain-15", legend = T, plot.title = "Louvain_Jaccard Graph-based Clustering (15 NNs)", alpha = 1)
plotDim(obj, "tfa")
markersAUCPR(obj, "1", clustering = "Louvain-15", exp.thresh = 1.0, frac.must.express = 0.6)
plotDim(obj, "neurod1")
plotDimHighlight(obj, clustering = "Louvain-15", cluster = "29", legend = F)

#Subset object to eliminate outliers
cells.keep <- setdiff(colnames(obj@logupx.data), c(outliers))
obj <- urdSubset(obj, cells.keep = cells.keep)

cells.exclude <- WhichCells(obj.seurat, idents = c("1", "40"))
cells.all <- setdiff(WhichCells(obj.seurat), cells.exclude)
obj <- urdSubset(obj.all, cells.all)

saveRDS(obj, "~/Box Sync/zfext/results/06b-URD/obj_subsets/endoderm_superset/sculptedKNN/endoderm_URD_trimmed_v4.rds")
obj <- readRDS("~/Box Sync/zfext/results/06b-URD/obj_subsets/endoderm_superset/endoderm_URD_trimmed_v4.rds")
obj <- calcDM(obj, knn = 100, sigma.use = 12)

saveRDS(obj, file="~/Box Sync/zfext/results/06b-URD/obj_subsets/endoderm_superset/sculptedKNN/endoderm_URD_withDM_v4.rds")
obj <- readRDS(file="~/Box Sync/zfext/results/06b-URD/obj_subsets/endoderm_superset/sculptedKNN/endoderm_URD_withDM_v4.rds")
obj <- readRDS(file="~/Box Sync/zfext/results/06a-URD/2021-03 KNN Modification/endoderm_KNNsculpted_dm9.3.rds")

stage.colors <- c("darkblue", "#FFCCCC", "#99CC00", "#33CC00", "cyan3",
                  "gold", "goldenrod", "darkorange", "indianred1", "indianred2", "plum", "deepskyblue2",
                  "lightgrey")
plotDimArray(obj, reduction.use = "dm", dims.to.plot = 1:18, outer.title = "Diffusion Map (Sigma 14, 100 NNs): Stage", label="stage.group", alpha = 0.4, plot.title="", legend=T, discrete.colors = stage.colors)
plotDim(obj, "stage.group", transitions.plot = 10000, plot.title="Developmental stage (with transitions)")

saveRDS(obj, "~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD_NCC_derived_DM.rds")

#Calculate pseudotime
message(paste0(Sys.time(), ": Defining earliest stage as the root.cells"))
root.cells.1 <- rownames(obj@meta)[obj@meta$stage.group == " 36-46"]
root.cells.2 <- rownames(obj@meta)[obj@meta$stage.group == " 48-58"]
root.cells <- WhichCells(obj.smc, idents = c("26"))
plotDimHighlight(obj, "stage.group", " 60-70", plot.title = "Root is 24 hpf cells")

#Do the flood
message(paste0(Sys.time(), ": Do the pseudotime_flood"))
flood.result <- floodPseudotime(obj, root.cells = root.cells, n = 100, minimum.cells.flooded = 2, verbose = T)

#Save the result
message(paste0(Sys.time(), ": Saving the flood object"))
saveRDS(flood.result, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD_NCC_derived_flood.rds")
flood.result <- readRDS(file="~/Box Sync/zfext/results/06a-URD/obj_subsets/endoderm_superset/seuratToURD_endoderm/endoderm_ds_flood_sculptedKNN.rds")

#Process pseudotime floods
message(paste0(Sys.time(), ": Processing pseudotime_flood"))
obj <- floodPseudotimeProcess(obj, flood.result, floods.name = "pseudotime")

pseudotimePlotStabilityOverall(obj)

#Inspect pseudotime in Diffusion Plots
plotDimArray(object = obj, reduction.use = "dm", dims.to.plot = 1:14,
             label = "pseudotime", plot.title = "", outer.title = "Diffusion Map labeled by pseudotime",
             legend = F, alpha = 0.4)

#Inspect pseudotime using tSNE plot
message(paste0(Sys.time(), ": Plotting Pseudotime on tSNE plot"))
pdf(file=paste0(plot.path, sample, "_pseudotime.pdf"), width = 8, height = 8)
plotDim(obj, "pseudotime", plot.title = "Pseudotime")

#Plotting the distribution of pseudotime from each stage group
#Create a dataframe that contains both pseudotime and stage information
gg.data <- cbind(obj@pseudotime, obj@meta[rownames(obj@pseudotime), ])
#Plot
ggplot(gg.data, aes(x=pseudotime, color = stage.group, fill = stage.group)) + geom_density(alpha = 0.4)

#Plotting distances
message(paste0(Sys.time(), ": Plotting Distances"))
pdf(file=paste0(plot.path, sample, "_pseudotime.pdf"), width = 8, height = 8)
plotDists(obj, "pseudotime", "stage.group", plot.title="Pseudotime by stage")

#Save DM and PT information
saveRDS(obj, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_NCC_derived_withDMandPT.rds")
obj <- readRDS(file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_withDMandPT.rds")
obj <- readRDS(file="~/Box Sync/zfext/results/06b-URD/obj_subsets/endoderm_superset/seuratToURD/endoderm_ds_withDMandPT_scKNN.rds")

#Part 3 - Determining Tips
message(paste0(Sys.time(), ": Cropping the cells from the final stage_group"))
cells_120h <- grep("120", colnames(obj@logupx.data), value = T)
cells_120h <- colnames(obj@logupx.data)[obj@meta$stage.group == "120"]
cells_108h <- colnames(obj@logupx.data)[obj@meta$stage.group == " 108-118"]
cells.keep <- unlist(unique(list(cells_120h, cells_108h)))

obj_120h <- urdSubset(obj, cells.keep = cells.keep)

#Perform PCA/tSNE on final stage cells
message(paste0(Sys.time(), ": Load the variable genes specific to this stage"))
var.genes.120h <- scan("120_var.txt")
obj_120h@var.genes <- obj@var.genes

#Calculate PCA
message(paste0(Sys.time(), ": Calculating PCA"))
obj_120h <- calcPCA(obj_120h)

#Calculate tSNE
message(paste0(Sys.time(), ": Calculating tSNE"))
set.seed(18)
obj_120h <- calcTsne(obj_120h, perplexity = 20, theta = 0.5)

obj_120h <- graphClustering(obj_120h, num.nn = 50, do.jaccard = T, method = "Louvain")
obj_120h <- graphClustering(obj_120h, num.nn = c(5, 8, 10, 15, 30, 40), do.jaccard = T, method = "Louvain")
obj_120h <- graphClustering(obj_120h, num.nn = c(10, 15, 20, 30, 40, 50), do.jaccard = T, method = "Infomap")

clusterings <- c(paste0("Infomap-", c(10, 15, 20, 30, 40, 50)), paste0("Infomap-", c(10, 15, 20, 30, 40, 50)))
clusterings <- c(paste0("Louvain-", c(50, 40, 30, 20, 15, 10)))

for (c in clusterings) {
  plot(plotDim(obj_120h, c, legend = T, label.clusters = T))
}

#Looking at batch information
pdf(file=paste0(plot.path, sample, "_tSNE_batch.pdf"), width = 8, height = 8)
plotDim(obj_120h, "Louvain-30", plot.title = "Louvain-15_graph", label.clusters = T)
plotDim(obj_120h, "stage.group", plot.title = "Louvain-15_graph")
plotDimHighlight(obj_120h, clustering = "Louvain-15", cluster = "11", legend = F)

cells.bil <- whichCells(obj, "Louvain-15", "22")
clusters <- sort(unique(obj_120h@group.ids$`Louvain-15`))
pr.markers <- lapply(clusters, function(c) markersAUCPR(obj_120h, clust.1 = c, clustering = "Louvain-15", genes.use = obj_120h@var.genes))
names(pr.markers) <- clusters

plotDim(obj_120h, "kcnk18", plot.title="FOXA2 (Pharyngeal endoderm marker)")

#Make a set of data.frames to keep track during cluster assignment
I30.n <- length(unique(obj_120h@group.ids$`Louvain-15`))
I30.cluster.assignments <- data.frame(cluster = 1:I30.n, name = rep(NA, I30.n, name = rep(NA, I30.n), tip = rep(NA, I30.n)), row.names = 1:I30.n)

I20.n <- length(unique(obj_120h@group.ids$`Louvain-15`))
I20.cluster.assignments <- data.frame(cluster = 1:I20.n, name = rep(NA, I20.n, name = rep(NA, I20.n), tip = rep(NA, I20.n)), row.names = 1:I20.n)

plotDot(obj_120h, genes = c("ndufa4l2a", "nr4a1", "atp1a1b", "mmrn2a", "igfbp7", "atp1b4", "prrx1b", "kcnj8", "kcne4", "abcc9", "wfdc2", "fbln5", "pitx2", "lmx1bb",
                            "loxa", "gdf10a", "tgfb2", "ptgdsb.2", "olfm1b", "ccl20b", "myl9a", "acta2", "tagln", "ifitm1", "tlx1", "lama5", "tinagl1", "spns2", "rgs5b", "trdn", "rcn3", "tnfrsf9a", "golim4b", "rflnb"), clustering = "Louvain-15")

#Markers for early exocrine pancreas
plotDot(obj_120h, genes = c("c9", "pcsk6", "pdx1", "surf4", "scgn", "prox1a", "isl1", "muc5.1", "muc5.2", "muc5.3", "agr2", "cldnh", "malb", "chs1", "slc6a8", "muc2.1", "nkx2.1", "foxe1", "elovl2",
                            "bmp4", "anxa13", "faxdc2", "slc7a9", "chia.2", "ace2", "gcga", "lgals2a", "lyz", "sdc4", "cftr", "anxa3a", "anxa4", "anxa2b", "dck",
                            "srgn", "msna", "cd63"), clustering = "Louvain-15")

saveRDS(obj, file=paste0(obj.path, sample, "_withDMandPT.rds"))

I30.cluster.assignments["1", "name"] <- "stellate-cells" #Use as tip
I30.cluster.assignments["2", "name"] <- "progenitors" #Use as tips
I30.cluster.assignments["3", "name"] <- "vSMC-artery_1" #Use as tip
I30.cluster.assignments["4", "name"] <- "progenitors" #Use as tip
I30.cluster.assignments["5", "name"] <- "int_SMCs" #Use as tip
I30.cluster.assignments["6", "name"] <- "sb_SMCs" #Use as tips
I30.cluster.assignments["7", "name"] <- "vSMC_artery_2" #Use as tip
#I30.cluster.assignments["8", "name"] <- "fibroblasts" #Use as tip

#Generate final clusterings
message(paste0(Sys.time(), ": Combine clustering assignments from two clusterings"))
I30.cluster.assignments$clustering <- "Louvain-15"
I20.cluster.assignments$clustering <- "Louvain-15"
cluster.assignments <- rbind(I30.cluster.assignments)

#Remove any clusters that weren't assigned an identity
message(paste0(Sys.time(), ": Removing clusters without an identity"))
cluster.assignments <- cluster.assignments[!is.na(cluster.assignments$name), ]

#Renumber clusters
cluster.assignments$cluster.new <- 1:nrow(cluster.assignments)

#Create blank clusterings for obj_120h
obj_120h@group.ids$clusters.120h.name <- NA
obj_120h@group.ids$clusters.120h.num <- NA

#Copy cell identities over for each cluster
for (i in 1:nrow(cluster.assignments)) {
  cells <- cellsInCluster(obj_120h, clustering = cluster.assignments[i, "clustering"],
                          cluster = cluster.assignments[i, "cluster"])
  obj_120h@group.ids[cells, "clusters.120h.name"] <- cluster.assignments[i, "name"]
  obj_120h@group.ids[cells, "clusters.120h.num"] <- as.character(cluster.assignments[i, "cluster.new"])
}

#Transfer clusterings to main object
obj@group.ids$`tip.clusters` <- NA
obj@group.ids[rownames(obj_120h@group.ids), "tip.clusters"] <- obj_120h@group.ids$`Louvain-15`

obj@group.ids$`Louvain-15` <- NA
obj@group.ids[rownames(obj_120h@group.ids), "Louvain-15"] <- obj_120h@group.ids$`Louvain-15`

obj@group.ids$`Cluster` <- NA
obj@group.ids[rownames(obj_120h@group.ids), "Cluster"] <- obj_120h@group.ids$clusters.120h.name

obj@group.ids$`Cluster-Num` <- NA
obj@group.ids[rownames(obj_120h@group.ids), "Cluster-Num"] <- obj_120h@group.ids$clusters.120h.num

#Save objects
message(paste0(Sys.time(), ": Saving the 120h_seurat object"))
saveRDS(obj_120h, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_120h_Louvain-15.rds")
obj_120h <- readRDS(file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_120h_Louvain-15.rds")

message(paste0(Sys.time(), ": Saving the full object with 120h clustering added"))
saveRDS(obj, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD_Louvain-15.rds")
obj <- readRDS(file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD_Louvain-15.rds")

message(paste0(Sys.time(), ": Saving the data.frame with tips"))
write.csv(cluster.assignments, file = "~/Box Sync/zfext/results/06b-URD/obj_subsets/endoderm_superset/sculptedKNN/endoderm_tips-use_Louvain-15_scKNN_v4.csv")

#Plot tips in diffusion map
obj@group.ids$pop <- NA
obj@group.ids[cellsInCluster(obj, "Cluster", "Liver"), "pop"] <- "6"
plotDim(obj, label = "pop", plot.title = "Pancreas DCs 13 vs 14", reduction.use = "dm", dim.x = 4, dim.y = 5, 
        legend = F, alpha = 0.35, colors = stage.colors)

#PART-4
#Biased random walks

message(paste0(Sys.time(), ": Load previous saved object"))
obj <- readRDS(file="~/Box Sync/zfext/results/06b-URD/obj_subsets/endoderm_superset/seuratToURD/endoderm_URD_Louvain-15_scKNN_2.rds")

obj@group.ids[rownames(obj_120h@group.ids), "tip.clusters"] <- obj_120h@group.ids$`Louvain-15`

#Define parameters of logistic function to bias transition probablities
diffusion.logistic <- pseudotimeDetermineLogistic(obj, "pseudotime", optimal.cells.forward = 20, max.cells.back = 40, pseudotime.direction = "<",
                                                  do.plot = T, print.values = T)

#Create biased transition matrix
message(paste0(Sys.time(), ": Creating biased transition matrix"))
biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(obj, pseudotime = "pseudotime",
                                                        logistic.params = diffusion.logistic, pseudotime.direction = "<"))

#Define the root cells
message(paste0(Sys.time(), ": Defining root-cells"))
root.cells <- rownames(obj@meta)[obj@meta$stage.group == " 24-34"]
plotDimHighlight(obj, "stage.group", " 24-34")

#Define tip cells
message(paste0(Sys.time(), ": Defining tip-cells"))
tips <- setdiff(unique(obj@group.ids[, clustering]), NA)
this.tip <- tips[tip.to.walk]

#Simulate the biased random walks from each tip
message(paste0(Sys.time(), ": Simulating random walks from each tip"))
tip.walks <- simulateRandomWalksFromTips(obj, tip.group.id = "tip.clusters", root.cells = root.cells, transition.matrix = biased.tm, n.per.tip = 25000, root.visits = 1,
                                         max.steps = 5000, verbose = T)

biased.tm.good <- intersect(rownames(biased.tm), tip.cells)

walks <- lapply(rownames(table(obj@group.ids$`Cluster-Num`)), function(c) {
  # Exclude any tip cells that for whatever reason didn't end up in the
  # biased TM (e.g. maybe not assigned a pseudotime).
  tip.cells <- rownames(obj@group.ids)[which(obj@group.ids$tip.clusters == c)]
  tip.cells.good <- intersect(tip.cells, rownames(biased.tm))
  # Perform the random walk simulation
  this.walk <- simulateRandomWalk(start.cells = tip.cells.good, transition.matrix = biased.tm,
                                  end.cells = root.cells, n = 50000, end.visits = 1, verbose.freq = 1000,
                                  max.steps = 5000)
  return(this.walk)
})
names(tip.walks) <- rownames(table(obj@group.ids$`Cluster-Num`))

saveRDS(tip.walks, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_TipWalks_Louvain-15.rds")
tip.walks <- readRDS(file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_TipWalks_Louvain-15.rds")

#Process the biased random walks into visitation frequencies
message(paste0(Sys.time(), ": Processing the biased random walks"))
obj <- processRandomWalksFromTips(obj, tip.walks, verbose = T)

#Visualize visitation of cells from each tip
plotDim(obj, "tip.clusters", plot.title = "Cells in each tip")

plotDim(obj, "visitfreq.log.2", plot.title = "Visitation frequency from tip-1 (log10)", transitions.plot = 10000)
plotDim(obj, "visitfreq.log.1", plot.title = "Visitation frequency from tip-2 (log10)", transitions.plot = 10000)

saveRDS(obj, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_withWalks_Louvain-15.rds")
obj <- readRDS(file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_withWalks_Louvain-15.rds")

tips.walked <- setdiff(unique(obj@group.ids$`Cluster-Num`), NA)

#PART 5 - Building URD Tree
library(URD)
library(rgl)

#Set up knitr to capture rgl output
rgl::setupKnitr()

obj <- readRDS(file="~/Box Sync/zfext/results/06b-URD/obj_subsets/endoderm_superset/seuratToURD/endoderm_withWalks_Louvain-15_scKNN_3.rds")

#Load tip cells
message(paste0(Sys.time(), ": Loading tip cells"))
obj <- loadTipCells(obj, tips = "tip.clusters")

#Combine a few sets of tips where we walked from 2 groups of cells and are intermixed in the diffusion map
#Liver cells
obj <- combineTipVisitation(obj, "2", "3", "3")
#Exocrine Pancreas
obj <- combineTipVisitation(obj, "4", "5", "5")
#Anterior Intestine
obj <- combineTipVisitation(obj, "3", "15", "15")
#Posterior Intestine
obj <- combineTipVisitation(obj, "1", "9", "9")
#Biliary cells
obj <- combineTipVisitation(obj, "11", "19", "11")
#Digestive Gut
#obj <- combineTipVisitation(obj, "5", "10", "5")
#obj <- combineTipVisitation(obj, "30", "31", "30")


#Build the actual tree
message(paste0(Sys.time(), ": Decide on tips to use for tree construction"))
tips.to.exclude <- c("2", "4")
tips.to.use <- setdiff(as.character(1:7), tips.to.exclude)

#Build the tree
obj.built <- buildTree(object = obj, pseudotime = "pseudotime", divergence.method = "ks",
                       tips.use = 1:3, weighted.fusion = T, use.only.original.tips = T,
                       cells.per.pseudotime.bin = 25, bins.per.pseudotime.window = 8, minimum.visits = 1,
                       visit.threshold = 0.7, p.thresh = 0.1, save.breakpoint.plots = NULL, dendro.node.size = 100,
                       min.cells.per.segment = 10, min.pseudotime.per.segment = 0.0001, verbose = F)

obj.tree <- buildTree(obj, pseudotime = "pseudotime", tips.use = tips.to.use, divergence.method = "preference", cells.per.pseudotime.bin = 10, 
                      bins.per.pseudotime.window = 8, save.all.breakpoint.info = T, p.thresh = 0.001, verbose = F)

#Name the tips
tip.names <- unique(obj@group.ids[, c("Cluster", "Cluster-Num")])
tip.names <- tip.names[complete.cases(tip.names), ]
obj.tree <- nameSegments(obj.tree, segments = tip.names$`Cluster-Num`, segment.names = tip.names$`Cluster`)
obj.built <- nameSegments(obj.built, segments = tip.names$`Cluster-Num`, segment.names = tip.names$`Cluster`)

plotTree(obj.built, label.segments = T, "stage.group")
plotTree(urd.build, label.segments = T, "stage.group", cell.alpha = 0.8, cell.size = 1.4)

saveRDS(obj.tree, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/vSMCs_artery_TREE.rds")
obj.tree <- readRDS(file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/vSMCs_artery_TREE.rds")

##Build a force-directed layout of the trajectory
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
robustly.visited.cells <- visitation[visitation$visit >= 0.6, "cell"]
# Since some tips of the tree were combined in their entirety, get the terminal segments to use as the tips of the force-directed layout.
final.tips <- segTerminal(urd.tree)

## Choose the number of nearest neighbors
#Here, we generate the force-directed layout by varying the number of nearest neighbors (num.nn).

par(mfrow=c(3,2))
for (k in c(110:120)){
  urd.tree <- treeForceDirectedLayout(urd.tree, num.nn = k, cells.to.do = robustly.visited.cells, cut.outlier.cells = NULL,
                                      cut.outlier.edges = NULL,
                                      cut.unconnected.segments = 2, min.final.neighbors = 4,
                                      tips = final.tips, verbose = T)
  
  plotTreeForce(urd.tree, "segment", alpha=1)
}

#OPTIMAL num.nn = 110

#Change the number of unconnected segments to find the optimal number of unconnected segments
for (j in c(2,4,6,8)){
  urd.tree <- treeForceDirectedLayout(urd.tree, num.nn = 110, method = "fr", cells.to.do = robustly.visited.cells,
                                      cut.outlier.cells = NULL,
                                      cut.outlier.edges = NULL, max.pseudotime.diff = NULL,
                                      cut.unconnected.segments = j, min.final.neighbors = 4,
                                      tips = final.tips, coords = "auto", start.temp = NULL,
                                      density.neighbors = 10, plot.outlier.cuts = F, verbose = T)
  
  plotTreeForce(urd.tree, "segment", alpha=1)
}

#OPTIMAL unconnected segments = 2

## Calculate the layout to be used for publication.
# use the optimal parameters
urd.tree <- treeForceDirectedLayout(urd.tree, num.nn = 110, cells.to.do = robustly.visited.cells, cut.outlier.cells = NULL,
                                    cut.outlier.edges = NULL,
                                    cut.unconnected.segments = 2, min.final.neighbors = 4,
                                    tips = final.tips, verbose = T)
plotTreeForce(urd.tree, "segment", alpha= 1)
plotTreeForce(urd.build, "stage.group", alpha = T)

## Rotate the tree and save the view
urd.build <- plotTreeForceStore3DView(urd.build, "View2")
saveRDS(urd.build, file=paste0("~/Box/zfext/annotations_celltype_curated_newMama/mural-cells/obj_urd/mural_cells_urd_treeView2.rds"))
# Load previous saved object
urd.build <- readRDS(file=paste0("~/Box/zfext/annotations_celltype_curated_newMama/mural/obj_urd/mural_cells_urd_treeView2.rds"))
plotTreeForce(urd.build, "segment", alpha=1, alpha.fade=0.08, size=4, density.alpha=T, label.tips=F, view = "View2")
rglwidget()


##Plot differentially expressed TFs for EECs and endo_panc
##Common markers
fire.with.grey <- c("#CECECE", "#DDC998", RColorBrewer::brewer.pal(9, "YlOrRd")[3:9])
genes.plot <- c("fsta")
for (gene in genes.plot) {
  plotTreeForce(urd.build, gene, alpha=1.8, alpha.fade=0.4, size= 8, density.alpha=T, label.tips=F, view = "View2",
                colors = pond.with.grey)
}

library(rgl)

plotTreeForce(urd.build, "stage.group", alpha=0.4, alpha.fade=0.08, size=6, density.alpha=T, label.tips=F, view = "View2", colors = stage.colors.new)
rgl.snapshot(filename = "mural-cells_FDL_fsta.png", fmt = "png")
dev.off()

fire.with.grey <- c("#CECECE", "#DDC998", RColorBrewer::brewer.pal(9, "YlOrRd")[3:9])
pond.with.grey <- c("#CECECE", "#CBDAC2", RColorBrewer::brewer.pal(9, "YlGnBu")[3:9])

genes.plot <- c("tead3a")

for (gene in genes.plot) {
  plotTreeForce(urd.build, gene, alpha=1.8, alpha.fade=0.4, size= 8, density.alpha=T, label.tips=F, view = "View2",
                colors = pond.with.grey)
}

rgl.snapshot(filename = paste0("mural_FDL_tree_", genes.plot, ".png"), fmt = "png")
dev.off()


