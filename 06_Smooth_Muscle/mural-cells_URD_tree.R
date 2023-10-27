##This script is to generate a URD tree of non-neural crest derived cells (foxc1a/foxc1b- and prrx1a/1b-) from this atlas that included the following clusters:
##C13: intestinal SMC progenitors
##C10: intestinal SMCs, longitudinal
##C8: intestinal SMCs, circular
##C11: vaSMCs, arterial 1
##C21: vaSMCs, arterial 2
##C3: myofibroblasts
##C16: viSMCs, col6a1+/rergla+/corin+

library(Seurat)
library(URD)

obj1 <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/mural/obj_seurat/mural_seurat.rds")
markers.mura <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/mural/obj_seurat/mural_markers.rds")

DimPlot(obj1, label = T)
DimPlot(obj1, group.by = "stage.group")

##Plot feature plots for foxc1a, foxc1b, prrx1a and prrx1b
FeaturePlot(obj1, features = c("foxc1a", "foxc1b", "prrx1a", "prrx1b"))



####################### CREATE URD OBJECT AND REMOVE OUTLIERS #################################
##Crop cells belonging to clusters pericyte or SMC clusters that do not express the above  markers - representing non-neural crest derived

##Get cells that belong to these clusters not expressing neural crest markers
cells.ident <- WhichCells(obj1, idents = c("13", "10", "8", "21", "11", "3"))

##Get cells that express neural crest markers
cells.foxc1b <- WhichCells(obj1, expression = foxc1b > 1)
cells.foxc1a <- WhichCells(obj1, expression = foxc1a > 1)
cells.prrx1b <- WhichCells(obj1, expression = prrx1b > 1)
cells.prrx1a<- WhichCells(obj1, expression = prrx1a > 1)

##Remove cells from the clusters above that may express neural crest genes
cells.keep <- setdiff(cells.ident, unlist(unique(list(cells.foxc1a, cells.foxc1b, cells.prrx1a, cells.prrx1b))))

##Crop the non-skeletal muscle seurat object to create a smaller object with just these non-neural crest derived cells
obj.crop <- subset(obj1, cells = cells.keep)

##Convert seurat object to URD object

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

##Convert cropped Seurat object to URD object
obj <- seuratToURD2(obj.crop)

##Calculate PCA using the seurat variable genes
#obj <- calcPCA(obj)
#vg.dups <- duplicated(as.data.frame(as.matrix(t(obj@logupx.data[obj@var.genes,
#]))))
#if (length(which(vg.dups)) > 0) {
#  print(paste("Removing", length(which(vg.dups)), "cell(s) with duplicated variable gene expression."))
#  not.dup.cells <- colnames(obj@logupx.data)[!vg.dups]
#  obj <- urdSubset(obj, not.dup.cells)
#}

##Calculate TSNE
obj <- calcTsne(obj, perplexity = 30, theta = 0.5)

##Perform Louvain clustering at different resolutions
obj <- graphClustering(obj, dim.use = "pca", num.nn = c(10, 15, 20, 30), do.jaccard = T, method = "Louvain")

##Plot graph in URD colored with stage groups
plotDim(obj, "stage.group", legend = T, plot.title = "Developmental Stage", alpha = 2.4)

##Transfer the UMAP coordinates from the Seurat object so that the UMAP plots remain same
##Get UMAP coordinates from the Seurat object
umap.coord <- obj1@reductions$umap@cell.embeddings
##Check coordinates
head(umap.coord)
#Check coordinates in the URD object
head(obj@tsne.y)
##Rename the colnames of the Seurat  UMAP coordinates dataframe to match that of the URD object
colnames(umap.coord) <- c("tSNE1", "tSNE2")
##Make sure the two coordinates are the same class
class(umap.coord)
class(obj@tsne.y)

##Replace  the TSNE coordinates with the Seurat UMAP cooridnates (after converting to a dataframe)
obj@tsne.y <- as.data.frame(umap.coord)
plotDim(obj, "stage.group", label.clusters = T)
plotTree(obj, label.segments = T)

# Save object
message(paste0(Sys.time(), ": Saving"))
saveRDS(obj, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD.rds")
obj <- readRDS(file="~/Box Sync/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD.rds")

##Calculate KNN on the URD object
obj <- calcKNN(obj)

##Remove outliers
outliers <- knnOutliers(obj, nn.1 = 1, nn.2 = 10, x.max = 25, slope.r = 1.1, int.r = 4, slope.b = 1.85, int.b = 8, title = "Identifying outliers by k-NN distance")                                                                             
gridExtra::grid.arrange(grobs=list(#Plot some apoptotic markers
  plotDim(obj, "isg15", alpha = 0.4, point.size = 0.5), 
  plotDim(obj, "foxo3b", alpha = 0.4, point.size = 0.5),
  plotDim(obj, "gadd45aa", alpha = 0.4, point.size = 0.5)))

#Subset object to eliminate outliers
cells.keep <- setdiff(colnames(obj@logupx.data), c(outliers))
obj <- urdSubset(obj, cells.keep = cells.keep)

saveRDS(obj, "~/Box Sync/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD_trimmed.rds")
obj <- readRDS("~/Box Sync/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD_trimmed.rds")



#################### CALCULATE DIFFUSION MAP AND PSEUDOTIME #######################################
## Calculate diffusion map
obj <- calcDM(obj, knn = 100, sigma.use = 12)

##Define stage colors for better visualization
stage.colors <- c("darkblue", "#FFCCCC", "#99CC00", "#33CC00", "cyan3",
                  "gold", "goldenrod", "darkorange", "indianred1", "indianred2", "plum", "deepskyblue2",
                  "lightgrey")

#@Plot  the diffusion map by stage using the newly defined stage colors
plotDimArray(obj, reduction.use = "dm", dims.to.plot = 1:18, outer.title = "Diffusion Map (Sigma 12, 100 NNs): Stage", label="stage.group", alpha = 0.4, plot.title="", legend=T, discrete.colors = stage.colors)

##Plot transition plot
plotDim(obj, "stage.group", transitions.plot = 10000, plot.title="Developmental stage (with transitions)")

##Save object
saveRDS(obj, "~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD_DM.rds")

#Calculate pseudotime 
## URD requires a starting point or "root" for determining pseudotime. Here we used all cells from the earliest timepoint represented in this ##atlas as the root

message(paste0(Sys.time(), ": Defining earliest stage as the root.cells (i.e., 36-58hpf"))
root.cells.1 <- rownames(obj@meta)[obj@meta$stage.group == " 36-46"]
root.cells.2 <- rownames(obj@meta)[obj@meta$stage.group == " 48-58"]
root.cells <- unlist(unique(list(root.cells.1, root.cells.2)))
plotDimHighlight(obj, "stage.group", " 36-46", plot.title = "Root is 36-58 hpf cells")

#Do the flood - run graph  search simulations to determine pseudotime
message(paste0(Sys.time(), ": Do the pseudotime_flood"))
flood.result <- floodPseudotime(obj, root.cells = root.cells, n = 100, minimum.cells.flooded = 2, verbose = T)

#Save the result
message(paste0(Sys.time(), ": Saving the flood object"))
saveRDS(flood.result, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD_flood.rds")

#Process the graph search simulations to determine the pseudotime value for each cell
message(paste0(Sys.time(), ": Processing pseudotime_flood"))
obj <- floodPseudotimeProcess(obj, flood.result, floods.name = "pseudotime")

##Inspect whether enough simulations have been run resulting in the overall change in pseudotime to reach an asymptote. If it doesn't then pseudotime needs to be recalculated with a higher n.  
pseudotimePlotStabilityOverall(obj)

#Inspect pseudotime in Diffusion Plots
plotDimArray(object = obj, reduction.use = "dm", dims.to.plot = 1:14,
             label = "pseudotime", plot.title = "", outer.title = "Diffusion Map labeled by pseudotime",
             legend = F, alpha = 0.4)

#Inspect pseudotime using tSNE plot or on the UMAP projections transferred from the seurat object
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

#Save the URD object with diffusion map and pseudotime calculated
saveRDS(obj, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_withDMandPT.rds")
obj <- readRDS(file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_withDMandPT.rds")



########################### DETERMINING TIPS ###################################################
#Part 3 - Determining Tips - We used cells from the final stage (i.e., 108-118 hpf and 120 hpf) as the tips for performing biased random walks
##Here we define cells belonging to each of those clusters

message(paste0(Sys.time(), ": Cropping the cells from the final stage_group"))
cells_120h <- colnames(obj@logupx.data)[obj@meta$stage.group == "120"]
cells_108h <- colnames(obj@logupx.data)[obj@meta$stage.group == " 108-118"]
cells.keep <- unlist(unique(list(cells_120h, cells_108h)))

##Subset the URD object to just include the  120 hpf cells
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

##Perform Louvain clustering at  various resolutions
obj_120h <- graphClustering(obj_120h, num.nn = 50, do.jaccard = T, method = "Louvain")
obj_120h <- graphClustering(obj_120h, num.nn = c(5, 8, 10, 15, 30, 40), do.jaccard = T, method = "Louvain")
obj_120h <- graphClustering(obj_120h, num.nn = c(10, 15, 20, 30, 40, 50), do.jaccard = T, method = "Infomap")

## plotting at  various clustering resolutions
clusterings <- c(paste0("Infomap-", c(10, 15, 20, 30, 40, 50)), paste0("Infomap-", c(10, 15, 20, 30, 40, 50)))
clusterings <- c(paste0("Louvain-", c(50, 40, 30, 20, 15, 10)))

for (c in clusterings) {
  plot(plotDim(obj_120h, c, legend = T, label.clusters = T))
}

##Final clustering used - Louvain-15
#Looking at batch information
pdf(file=paste0(plot.path, sample, "_tSNE_batch.pdf"), width = 8, height = 8)
plotDim(obj_120h, "Louvain-15", plot.title = "Louvain-15_graph", label.clusters = T)
plotDim(obj_120h, "stage.group", plot.title = "Louvain-15_graph")
plotDimHighlight(obj_120h, clustering = "Louvain-15", cluster = "11", legend = F)


##Calculate differentially expressed genes between the clusters of the 108-120 hpf object
pr.markers <- lapply(clusters, function(c) markersAUCPR(obj_120h, clust.1 = c, clustering = "Louvain-15", genes.use = obj_120h@var.genes))
names(pr.markers) <- clusters

##Plot a few markers
plotDim(obj_120h, "kcnk18", plot.title="KCNK18 (circular iSMC marker)")

#Make a set of data.frames to keep track during cluster assignment
I30.n <- length(unique(obj_120h@group.ids$`Louvain-15`))
I30.cluster.assignments <- data.frame(cluster = 1:I30.n, name = rep(NA, I30.n, name = rep(NA, I30.n), tip = rep(NA, I30.n)), row.names = 1:I30.n)

I20.n <- length(unique(obj_120h@group.ids$`Louvain-15`))
I20.cluster.assignments <- data.frame(cluster = 1:I20.n, name = rep(NA, I20.n, name = rep(NA, I20.n), tip = rep(NA, I20.n)), row.names = 1:I20.n)

##Plot a dotplot with differentially expressed genes to aid in cluster annotations of the 108-120 hpf object
plotDot(obj_120h, genes = c("ndufa4l2a", "nr4a1", "atp1a1b", "mmrn2a", "igfbp7", "atp1b4", "prrx1b", "kcnj8", "kcne4", "abcc9", "wfdc2", "fbln5", "pitx2", "lmx1bb",
                            "loxa", "gdf10a", "tgfb2", "ptgdsb.2", "olfm1b", "ccl20b", "myl9a", "acta2", "tagln", "ifitm1", "tlx1", "lama5", "tinagl1", "spns2", "rgs5b", "trdn", "rcn3", "tnfrsf9a", "golim4b", "rflnb"), clustering = "Louvain-15")


##Assign cluster identities based on markers
I30.cluster.assignments["1", "name"] <- "viSMCs_corin+/rergla+" #Use as tip
I30.cluster.assignments["2", "name"] <- "vSMC-artery_1" #Use as tips
I30.cluster.assignments["3", "name"] <- "vSMC-artery_1" #Use as tip
I30.cluster.assignments["4", "name"] <- "progenitors" #Use as tip
I30.cluster.assignments["5", "name"] <- "iSMCs_circular" #Use as tip
I30.cluster.assignments["6", "name"] <- "iSMCs_longitudinal" #Use as tips
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

#Save objects with the updated changes
message(paste0(Sys.time(), ": Saving the 120h_seurat object"))
saveRDS(obj_120h, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_120h_Louvain-15.rds")
obj_120h <- readRDS(file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_120h_Louvain-15.rds")

message(paste0(Sys.time(), ": Saving the full object with 120h clustering added"))
saveRDS(obj, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD_Louvain-15.rds")
obj <- readRDS(file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD_Louvain-15.rds")

message(paste0(Sys.time(), ": Saving the data.frame with tips"))
write.csv(cluster.assignments, file = "~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD_Louvain-15_clusters.csv")

#Plot tips in diffusion map
obj@group.ids$pop <- NA
obj@group.ids[cellsInCluster(obj, "Cluster", "iSMC_longitudinal"), "pop"] <- "6"
plotDim(obj, label = "pop", plot.title = "iSMCs DCs 5 vs 6", reduction.use = "dm", dim.x = 4, dim.y = 5, 
        legend = F, alpha = 0.35, colors = stage.colors)



################################# PERFORM BIASED RANDOM WALKS #######################################################

#PART-4
#Biased random walks
message(paste0(Sys.time(), ": Load previous saved object"))
obj <- readRDS(file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_URD_Louvain-15.rds")
obj_120h <- readRDS(file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_120h_Louvain-15.rds")

obj@group.ids[rownames(obj_120h@group.ids), "tip.clusters"] <- obj_120h@group.ids$`Louvain-15`

#Define parameters of logistic function to bias transition probablities
diffusion.logistic <- pseudotimeDetermineLogistic(obj, "pseudotime", optimal.cells.forward = 10, max.cells.back = 20, pseudotime.direction = "<",
                                                  do.plot = T, print.values = T)

#Create biased transition matrix
message(paste0(Sys.time(), ": Creating biased transition matrix"))
biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(obj, pseudotime = "pseudotime",
                                                        logistic.params = diffusion.logistic, pseudotime.direction = "<"))

#Define the root cells
message(paste0(Sys.time(), ": Defining earliest stage as the root.cells (i.e., 36-58hpf"))
root.cells.1 <- rownames(obj@meta)[obj@meta$stage.group == " 36-46"]
root.cells.2 <- rownames(obj@meta)[obj@meta$stage.group == " 48-58"]
root.cells <- unlist(unique(list(root.cells.1, root.cells.2)))

#Define tip cells
message(paste0(Sys.time(), ": Defining tip-cells"))
tips <- setdiff(unique(obj@group.ids[, clustering]), NA)
this.tip <- tips[tip.to.walk]

#Simulate the biased random walks from each tip
message(paste0(Sys.time(), ": Simulating random walks from each tip"))
tip.walks <- simulateRandomWalksFromTips(obj, tip.group.id = "tip.clusters", root.cells = root.cells, transition.matrix = biased.tm, n.per.tip = 25000, root.visits = 1,
                                         max.steps = 5000, verbose = T)

biased.tm.good <- intersect(rownames(biased.tm), tip.cells)

##Alternatively, the random walks can be run as a loop for each cluster
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

##Save the random walk output
saveRDS(tip.walks, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_TipWalks_Louvain-15.rds")
tip.walks <- readRDS(file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_TipWalks_Louvain-15.rds")

#Process the biased random walks into visitation frequencies
message(paste0(Sys.time(), ": Processing the biased random walks"))
obj <- processRandomWalksFromTips(obj, tip.walks, verbose = T)

#Visualize visitation of cells from each tip
plotDim(obj, "tip.clusters", plot.title = "Cells in each tip")

##Visualize visitation from user-defined tips
plotDim(obj, "visitfreq.log.2", plot.title = "Visitation frequency from tip-1 (log10)", transitions.plot = 10000)
plotDim(obj, "visitfreq.log.1", plot.title = "Visitation frequency from tip-2 (log10)", transitions.plot = 10000)

saveRDS(obj, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_withWalks_Louvain-15.rds")
obj <- readRDS(file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_withWalks_Louvain-15.rds")

tips.walked <- setdiff(unique(obj@group.ids$`Cluster-Num`), NA)


############################### BUILDING THE URD TREE #########################################################
#PART 5 - Building URD Tree
library(URD)
library(rgl)

#Set up knitr to capture rgl output
rgl::setupKnitr()

##Read in object with Walks calculated if starting here
obj <- readRDS(file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/mural-cells_withWalks_Louvain-15.rds")

##Tree building is destructive, so create a copy of the current URD object
object.tree <- obj

#Load tip cells
message(paste0(Sys.time(), ": Loading tip cells"))
obj <- loadTipCells(obj, tips = "tip.clusters")

#Combine a few sets of tips where we walked from 2 groups of cells and are intermixed in the diffusion map
#vaSMC, artery 1
obj <- combineTipVisitation(obj, "2", "3", "3")
#iSMC, progenitors
obj <- combineTipVisitation(obj, "4", "5", "5")

#Build the actual tree
message(paste0(Sys.time(), ": Decide on tips to use for tree construction"))
tips.to.exclude <- c("2", "4")
tips.to.use <- setdiff(as.character(1:7), tips.to.exclude)

#Build the tree
#obj.built <- buildTree(object = obj, pseudotime = "pseudotime", divergence.method = "ks",
                      # tips.use = 1:3, weighted.fusion = T, use.only.original.tips = T,
                      # cells.per.pseudotime.bin = 25, bins.per.pseudotime.window = 8, minimum.visits = 1,
                      # visit.threshold = 0.7, p.thresh = 0.1, save.breakpoint.plots = NULL, dendro.node.size = 100,
                      # min.cells.per.segment = 10, min.pseudotime.per.segment = 0.0001, verbose = F)

##Build the URD tree - "preference method"
obj.tree <- buildTree(obj, pseudotime = "pseudotime", tips.use = tips.to.use, divergence.method = "preference", cells.per.pseudotime.bin = 10, 
                      bins.per.pseudotime.window = 8, save.all.breakpoint.info = T, p.thresh = 0.001, verbose = F)

#Name the tips of the resultant URD tree
tip.names <- unique(obj@group.ids[, c("Cluster", "Cluster-Num")])
tip.names <- tip.names[complete.cases(tip.names), ]
obj.tree <- nameSegments(obj.tree, segments = tip.names$`Cluster-Num`, segment.names = tip.names$`Cluster`)
#obj.built <- nameSegments(obj.built, segments = tip.names$`Cluster-Num`, segment.names = tip.names$`Cluster`)

##Plot the URD tree
#plotTree(obj.built, label.segments = T, "stage.group")
plotTree(obj.tree, label.segments = T, "stage.group")

##Save the URD  tree
saveRDS(obj.tree, file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/vSMCs_artery_TREE.rds")
obj.tree <- readRDS(file="~/Box/zfext/results/06b-URD/obj_subsets/mural-cells/vSMCs_artery_TREE.rds")

##Save the branching tree structure to a separate file
png("mural_cell_mesoderm_TREE_colored.png", width = 5*dpi, height = 5*dpi)
plotTree(obj.tree, "stage.nice", label.segments = F, plot.cells = T, legend = F, cell.alpha = 0.8, cell.size = 1.8, label.x = F, continuous.colors = stage.colors.new)
dev.off()

##Plot some genes on the tree
genes.plot <- c("ndufa4l2a", "il13ra2", "desmb")
for (gene in genes.plot) {
  plot(plotTree(obj.tree, gene))
}

##Plotting genes on the URD tree
gridExtra::grid.arrange(grobs = lapply(c("abcc9", "mcamb", "ednraa", "cd248a"), plotTree,
                                       object = obj.tree, label.x = F, plot.cells = T), ncol = 2)

gridExtra::grid.arrange(grobs = lapply(c("acta2", "ndufa4l2a", "bgna", "pdgfrb"), plotTree,
                                       object = obj.tree, label.x = F, plot.cells = F, cell.alpha = 0.8, cell.size = 1.0), ncol = 2)


#PART-2
#Create branchpoint plots to understand cellular dynamics between individual cell types
##np.layout <- branchpointPreferenceLayout(obj.tree, pseudotime = "pseudotime", lineages.1 = "3", lineages.2 = "7", parent.of.lineages = c("9") , opposite.parent = 1, min.visit = 2)
#plotBranchpoint(obj.tree, np.layout, label = "segment", point.alpha = 0.8, populations = c("aSMC-1", "aSMC-2"),
 #               pt.lim = c(0.7, 0.3), xlab = "", ylab = "", legend = T, axis.lines = F,
##              fade.low = 0, title = "Stage.group")

#Examine gene expression in the branchpoint points
#genes.plot <- c("acta2", "macrod2", "bgna", "ndufa4l2a", "emilin2a", "pthlha")
#genes.plot <- c("gpx3", "gpx4a", "egr1")
#branch.plots <- lapply(genes.plot, function(gene) plotBranchpoint(obj.tree, np.layout, label = gene, point.alpha = 1, populations = c("aSMC-1", "aSMC-2"),
 #                                                                 pt.lim = c(0.7, 0.3), xlab = "", ylab = "", title = gene, legend = F,
#                                                                  axis.lines = F, fade.low = 0.66))
#gridExtra::grid.arrange(grobs = branch.plots, ncol = 3)

#Manual refinement - Rename segment names
message(paste0(Sys.time(), ": Descriptive names to be used on the dendrogram"))
new.seg.names <- c("iSMCs_longitudinal", "iSMCs_circular", "vSMC-artery_1", "vSMC-artery_2", "smooth muscle (visceral?)")
segs.to.name <- c("5", "6", "3", "7", "1")
obj.tree <- nameSegments(obj.tree, segments = segs.to.name, segment.names = new.seg.names)



############################## BUILDING A FORCE-DIRECTED LAYOUT OF THE URD TREE ############################################
##Build a force-directed layout of the trajectory

## Choose cells that were visited more robustly
urd.tree <- obj.tree
# Create a data frame to measure cell visitation
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
##Plot force directed layout by segment
plotTreeForce(urd.tree, "segment", alpha= 1)

##Plot force directed layout by stage groups - Figure 5H
plotTreeForce(urd.build, "stage.group", alpha = T)

## Rotate the tree and save the view
urd.build <- plotTreeForceStore3DView(urd.build, "View2")
saveRDS(urd.build, file=paste0("~/Box/zfext/annotations_celltype_curated_newMama/mural-cells/obj_urd/mural_cells_urd_treeView2.rds"))

# Load previous saved object
urd.build <- readRDS(file=paste0("~/Box/zfext/annotations_celltype_curated_newMama/mural/obj_urd/mural_cells_urd_treeView2.rds"))
plotTreeForce(urd.build, "segment", alpha=1, alpha.fade=0.08, size=4, density.alpha=T, label.tips=F, view = "View2")
rglwidget()


##Plot differentially expressed markers for each branch
##Common markers
## Plot genes on the FDL projection - shown in Figure 5I, 
pond.with.grey <- c("#CECECE", "#CBDAC2", RColorBrewer::brewer.pal(9, "YlGnBu")[3:9])
genes.plot <- c("fsta", "foxq1a", "foxq1b", "cremb", "itpr1a", "il13ra2")
for (gene in genes.plot) {
  plotTreeForce(urd.build, gene, alpha=1.8, alpha.fade=0.4, size= 8, density.alpha=T, label.tips=F, view = "View2",
                colors = pond.with.grey)
}

library(rgl)

plotTreeForce(urd.build, "stage.group", alpha=0.4, alpha.fade=0.08, size=6, density.alpha=T, label.tips=F, view = "View2", colors = stage.colors.new)
rgl.snapshot(filename = "mural-cells_FDL_fsta.png", fmt = "png")
dev.off()

##Define color schemes and try out different color schemes
fire.with.grey <- c("#CECECE", "#DDC998", RColorBrewer::brewer.pal(9, "YlOrRd")[3:9])
pond.with.grey <- c("#CECECE", "#CBDAC2", RColorBrewer::brewer.pal(9, "YlGnBu")[3:9])

##Plot differentially expressed genes for the two intestinal SMCs - Figure S6I-K
genes.plot <- c("tead3a", "foxf2a", "tcf21", "foxf1", "foxp4", "meis2a", "pbx3b")

for (gene in genes.plot) {
  plotTreeForce(urd.build, gene, alpha=1.8, alpha.fade=0.4, size= 8, density.alpha=T, label.tips=F, view = "View2",
                colors = pond.with.grey)
}

##Save FDL expression of genes using rgl package
rgl.snapshot(filename = paste0("mural_FDL_tree_", genes.plot, ".png"), fmt = "png")
dev.off()


