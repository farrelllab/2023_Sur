library(Seurat)
library(URD)

##Specify sample name
sample = "intestine"

##Specify paths
save.path <- "~/Box/zfext/results/06b-URD/obj_subsets/2022-11_intestine_trajectory/obj/"
plot.path <- "~/Box/zfext/results/06b-URD/plot_subsets/plots/"

##Load endoderm object - start with the endoderm object
obj1 <- readRDS(paste0(file = "~/Box/zfext/annotations_celltype_curated_newMama/", sample, "/obj_seurat/", sample, "_seurat.rds"))
markers <- readRDS(paste0(file = "~/Box/zfext/annotations_celltype_curated_newMama/", sample, "/obj_seurat/", sample, "_markers.rds"))
DimPlot(obj1, label = T)    
DimPlot(obj1, group.by = "stage.nice")

################################### CREATING URD OBJECT AND REMOVING OUTLIERS ##################################################

#Create function to transfer Seurat information to URD
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

#Creating URD object
obj <- seuratToURD2(obj1)

saveRDS(obj, paste0(save.path, sample, "_URD.rds"))
obj <- readRDS(paste0(save.path, sample, "_URD.rds"))

stages <- unique(obj@meta$stage.group)
cells.each.stage <- lapply(stages, function(stage) rownames(obj@meta)[which(obj@meta$stage.group == stage)])

#Calculate TSNE and graph clusterings
obj <- calcPCA(obj)
obj <- calcTsne(obj, perplexity = 30, theta = 0.5)
obj <- graphClustering(obj, dim.use = "pca", num.nn = c(15, 20, 30), do.jaccard = T, method = "Louvain")
obj <- graphClustering(obj, dim.use = "pca", num.nn = c(15, 20, 25, 30, 35, 40, 45), do.jaccard = T, method = "Infomap")

##Use the UMAP coordinates from the Seurat object
umap.coord <- obj1@reductions$umap@cell.embeddings
head(umap.coord)
head(obj@tsne.y)
colnames(umap.coord) <- c("tSNE1", "tSNE2")
class(umap.coord)
class(obj@tsne.y)

##Add them to the URD object
obj@tsne.y <- as.data.frame(umap.coord)
plotDim(obj, "cdx1b")

##Now let's subset the full endoderm object to just intestine cells
cells.intestine <- WhichCells(obj1, idents = c("9", "10", "11", "12", "13", "14", "15", "16", "17", "18"))
cells.cdx1b <- WhichCells(obj1, expression = cdx1b > 0.25)
cells.diff <- setdiff(cells.intestine, cells.cdx1b)
cells.rest <- cells.cdx1b[!(cells.cdx1b %in% cells.intestine)]

##Subset URD object to just intestine cells
obj <- urdSubset(obj, cells.keep = cells.intestine)

#Calculate PCA, TSNE and graph clusterings
obj <- calcPCA(obj)
obj <- calcTsne(obj, perplexity = 30, theta = 0.5)
obj <- graphClustering(obj, dim.use = "pca", num.nn = c(15, 20, 30, 40), do.jaccard = T, method = "Louvain")

sample = "intestine"
saveRDS(obj, paste0(save.path, sample, "_URD.rds"))
obj <- readRDS(paste0(save.path, sample, "_URD.rds"))

#Plotting the clustering and assessing tSNE plot
plotDim(obj, "stage.group", legend = T, plot.title = "Developmental Stage", alpha = 0.5)
plotDim(obj, "Louvain-20", legend = T, plot.title = "Louvain_Jaccard Graph-based Clustering (15 NNs)", alpha = 1)

#Calculating a k-nn graph
message(paste0(Sys.time(), ": Calculating a k-nearest neighbour graph"))
obj <- calcKNN(obj)


#Removing outliers
#Plot cells according to their distance to their nearest and 20th nearest neighbours, and identify those with unusually large distances.
outliers <- knnOutliers(obj, nn.1 = 1, nn.2 = 12, x.max = 23, slope.r = 1.2, int.r = 5, slope.b = 0.85, int.b = 7.5, title = "Identifying outliers by k-NN distance")       

gridExtra::grid.arrange(grobs=list(#Plot some apoptotic markers
  plotDim(obj, "isg15", alpha = 0.4, point.size = 0.5), 
  plotDim(obj, "foxo3b", alpha = 0.4, point.size = 0.5),
  plotDim(obj, "gadd45aa", alpha = 0.4, point.size = 0.5),
  #Figure out which clusters correspond to these cells
  plotDimHighlight(obj, clustering = "Louvain-20", cluster = "12", legend = F)))

#apoptotic.like.cells <- cellsInCluster(obj, "Louvain-20", "26")

#Subset object to eliminate outliers
cells.keep <- setdiff(colnames(obj@logupx.data), c(outliers))
obj <- urdSubset(obj, cells.keep = cells.keep)

##Saving trimmed object
saveRDS(obj, paste0(save.path, sample, "_URD_trimmed.rds"))

################################## DIFFUSION MAP AND PSEUDOTIME ###########################################

#Calculate diffusion map
obj <- calcDM(obj, knn = 100, sigma.use = NULL) ##sigma = 8.2

##Save the diffusion map under a variable
dm.8 <- obj@dm
stage.colors <- c("darkblue", "#FFCCCC", "#99CC00", "#33CC00", "cyan3",
                  "gold", "goldenrod", "darkorange", "indianred1", "indianred2", "plum", "deepskyblue2",
                  "lightgrey")

##Plot the diffusion map
plotDimArray(obj, reduction.use = "dm", dims.to.plot = 1:18, outer.title = "Diffusion Map (Sigma 8.2, 100 NNs): Stage", label="stage.group", alpha = 0.4, plot.title="", legend=T, discrete.colors = stage.colors)
plotDim(obj, "stage.group", transitions.plot = 10000, plot.title="Developmental stage (with transitions)")

saveRDS(obj, paste0(save.path, sample, "_URD_withDM.rds"))

#Calculate pseudotime
##Define the earliest timepoint cells as the "root". 
message(paste0(Sys.time(), ": Defining earliest stage as the root.cells"))
root.cells <- rownames(obj@meta)[obj@meta$stage.group == " 14-21"]

plotDim(obj, "cdx1b", plot.title="CDX1B (Intestine)")
plotDim(obj, "stage.group", legend = T, plot.title = "Developmental Stage", alpha = 0.5)
plotDimHighlight(obj, "stage.group", " 14-21", plot.title = "tip-cells")

#Do the flood
message(paste0(Sys.time(), ": Do the pseudotime_flood"))
flood.result <- floodPseudotime(obj, root.cells = root.cells, n = 100, minimum.cells.flooded = 2, verbose = T)

#Save the result
message(paste0(Sys.time(), ": Saving the flood object"))
saveRDS(flood.result, paste0(save.path, "intestine_URD_flood.rds"))
flood.result <- readRDS(paste0(save.path, sample, "_URD_flood.rds"))

#Process pseudotime floods
message(paste0(Sys.time(), ": Processing pseudotime_flood"))
obj <- floodPseudotimeProcess(obj, flood.result, floods.name = "pseudotime")
pseudotimePlotStabilityOverall(obj)

#Inspect pseudotime using tSNE plot
message(paste0(Sys.time(), ": Plotting Pseudotime on tSNE plot"))
pdf(file=paste0(plot.path, sample, "_pseudotime.pdf"), width = 8, height = 8)
plotDim(obj, "stage.group", legend = T, plot.title = "Developmental Stage", alpha = 0.5)
plotDim(obj, "pseudotime", plot.title = "Pseudotime")


#Plotting distances
message(paste0(Sys.time(), ": Plotting Distances"))
pdf(file=paste0(plot.path, sample, "_pseudotime.pdf"), width = 8, height = 8)
plotDists(obj, "pseudotime", "stage.group", plot.title="Pseudotime by stage")

gg.data <- cbind(obj@pseudotime, obj@meta[rownames(obj@pseudotime), ])
#Plot
ggplot(gg.data, aes(x=pseudotime, color = stage.group, fill = stage.group)) + geom_density(alpha = 0.4)

#Save object
message(paste0(Sys.time(), ": Saving object"))
saveRDS(obj, file=paste0(save.path, sample, "_URD_withDMandPT.rds"))
obj <- readRDS(paste0(save.path, sample, "_URD_withDMandPT.rds"))

################################# DETERMINING  TIPS ###############################################################
#Part 3 - Determining Tips
message(paste0(Sys.time(), ": Cropping the cells from the final stage_group"))
cells_120h <- rownames(obj@meta)[obj@meta$stage.group == " 120"]
cells.tuft <- WhichCells(obj1, idents = c("18")) ##Add the tuft-like cells to the tip clusters 
cells.total <- unlist(unique(list(cells_120h, cells.tuft)))
obj_120h <- urdSubset(obj, cells.keep = cells.total)

#Perform PCA/tSNE on final stage cells
message(paste0(Sys.time(), ": Load the variable genes specific to this stage"))
var.genes.120h <- scan(var.path, sample, "120_var.txt", what = "character")
obj_120h@var.genes <- obj@var.genes

#Calculate PCA
message(paste0(Sys.time(), ": Calculating PCA"))
obj_120h <- calcPCA(obj_120h)

#Calculate tSNE
message(paste0(Sys.time(), ": Calculating tSNE"))
set.seed(18)
obj_120h <- calcTsne(obj_120h, perplexity = 30, theta = 0.5)

##Perform clustering using different methods
obj_120h <- graphClustering(obj_120h, num.nn = 40, do.jaccard = T, method = "Louvain")
obj_120h <- graphClustering(obj_120h, num.nn = c(5, 8, 10, 15, 20, 30), do.jaccard = T, method = "Louvain")
obj_120h <- graphClustering(obj_120h, num.nn = c(10, 15, 20, 30, 40, 50), do.jaccard = T, method = "Infomap")

##Use the different clusterings to plot the  TSNE plot
clusterings <- c(paste0("Infomap-", c(10, 15, 20, 30, 40, 50)), paste0("Infomap-", c(10, 15, 20, 30, 40, 50)))
clusterings <- c(paste0("Louvain-", c(10, 15, 20, 30, 40)))

for (c in clusterings) {
  plot(plotDim(obj_120h, c, legend = T))
}

#Looking at batch information
pdf(file=paste0(plot.path, sample, "_tSNE_batch.pdf"), width = 8, height = 8)
plotDim(obj_120h, "Louvain-20", plot.title = "Louvain-20_graph", legend = T, label.clusters = T)

##The tuft-like cell cluster is really tiny and hence does not resolve as a tip cluster. So trying to force it to be a tip cluster
##Set tuft cells as a separate cluster
cells.3 <- whichCells(obj_120h, "Louvain-20", "3")
cells.tuft.in.3 <- cells.tuft[cells.tuft %in% cells.3]
obj_120h@group.ids[cells.tuft.in.3, "Louvain-20"] <- "9"
plotDim(obj_120h, "Louvain-20", plot.title = "Louvain-20_graph", legend = T, label.clusters = T)

##Calculate differential markers for clusters
clusters <- sort(unique(obj_120h@group.ids$`Louvain-20`))
pr.markers <- lapply(clusters, function(c) markersAUCPR(obj_120h, clust.1 = c, clustering = "Louvain-20", genes.use = obj_120h@var.genes))
names(pr.markers) <- clusters

##Plot few markers
plotDim(obj_120h, "foxa2", plot.title="FOXA2 (Pharyngeal endoderm marker)")

#Make a set of data.frames to keep track during cluster assignment
I30.n <- length(unique(obj_120h@group.ids$`Louvain-20`))
I30.cluster.assignments <- data.frame(cluster = 1:I30.n, name = rep(NA, I30.n, name = rep(NA, I30.n), tip = rep(NA, I30.n)), row.names = 1:I30.n)

I20.n <- length(unique(obj_120h@group.ids$`Louvain-20`))
I20.cluster.assignments <- data.frame(cluster = 1:I20.n, name = rep(NA, I20.n, name = rep(NA, I20.n), tip = rep(NA, I20.n)), row.names = 1:I20.n)

#plotDimHighlight(obj_120h, "Louvain-20", "33")
#plotDim(obj, "tip.clusters", plot.title = "Cells in each tip")
#plotDimHighlight(obj, "tip.clusters", "33")
#plotDim(obj, "sox9b")


##Assign cluster identities based on prior annotations
I30.cluster.assignments["1", "name"] <- "enterocyte-3" #Use as tip
I30.cluster.assignments["2", "name"] <- "enterocyte-1" #Use as tips
I30.cluster.assignments["3", "name"] <- "posterior_prog" #Don't use as tip
I30.cluster.assignments["4", "name"] <- "posterior-LREs" #Use as tip
I30.cluster.assignments["5", "name"] <- "best4+/otop2+" #Use as tip
I30.cluster.assignments["6", "name"] <- "enterocyte-2" #Use as tip
I30.cluster.assignments["7", "name"] <- "EECs" #Use as tip
I30.cluster.assignments["8", "name"] <- "goblet-cells" #Don't use as tip
I30.cluster.assignments["9", "name"] <- "tuft-like" #Use as tip

#Generate final clusterings
message(paste0(Sys.time(), ": Combine clustering assignments from two clusterings"))
I30.cluster.assignments$clustering <- "Louvain-20"
I20.cluster.assignments$clustering <- "Louvain-20"
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
obj@group.ids$`Louvain-20` <- NA
obj@group.ids[rownames(obj_120h@group.ids), "tip.clusters"] <- obj_120h@group.ids$`Louvain-20`

obj@group.ids$`Louvain-20` <- NA
obj@group.ids[rownames(obj_120h@group.ids), "Louvain-20"] <- obj_120h@group.ids$`Louvain-20`

obj@group.ids$`Cluster` <- NA
obj@group.ids[rownames(obj_120h@group.ids), "Cluster"] <- obj_120h@group.ids$clusters.120h.name

obj@group.ids$`Cluster-Num` <- NA
obj@group.ids[rownames(obj_120h@group.ids), "Cluster-Num"] <- obj_120h@group.ids$clusters.120h.num

#Save objects
message(paste0(Sys.time(), ": Saving the 120h_seurat object"))
saveRDS(obj_120h, file = paste0(save.path, "intestine_URD_120h.rds"))
saveRDS(obj_120h, file = paste0(save.path, "intestine_URD_120h_with_tuft.rds"))
obj_120h <- readRDS(file = paste0(save.path, "intestine_URD_120h.rds"))

message(paste0(Sys.time(), ": Saving the full object with 120h clustering added"))
saveRDS(obj, file = paste0(save.path, "intestine_URD_withTips.rds"))
saveRDS(obj, file = paste0(save.path, "intestine_URD_withTips_with_tuftCells.rds"))
obj <- readRDS(file= paste0(save.path, "intestine_URD_withTips_with_tuftCells.rds"))

##Plot the clustering with tip clusters
plotDim(obj, "tip.clusters", label.clusters = T)

message(paste0(Sys.time(), ": Saving the data.frame with tips"))
write.csv(cluster.assignments, file = paste0(save.path, "intestine_tips-use_Louvain-20.csv"))

#Plot tips in diffusion map
obj@group.ids$pop <- NA
obj@group.ids[cellsInCluster(obj, "Cluster", "tuft-like"), "pop"] <- "9"
plotDim(obj, label = "pop", plot.title = "Intestine DCs 1 vs 2", reduction.use = "dm", dim.x = 1, dim.y = 2, 
        legend = F, alpha = 0.35)



#################################### PERFORM BIASED RANDOM WALKS ########################################################
#PART-4
#Biased random walks
message(paste0(Sys.time(), ": Load previous saved object"))
obj <- readRDS(file= paste0(save.path, "endoderm_URD_withTips.rds"))
#Define parameters of logistic function to bias transition probablities
diffusion.logistic <- pseudotimeDetermineLogistic(obj, "pseudotime", optimal.cells.forward = 20, max.cells.back = 40, pseudotime.direction = "<",
                                                  do.plot = T, print.values = T)

#Create biased transition matrix
message(paste0(Sys.time(), ": Creating biased transition matrix"))
biased.tm <- pseudotimeWeightTransitionMatrix(obj, pseudotime = "pseudotime",
                                              logistic.params = diffusion.logistic, pseudotime.direction = "<")

#Define the root cells
message(paste0(Sys.time(), ": Defining root-cells"))
root.cells <- rownames(obj@meta)[obj@meta$stage.group == "14-21"]

#Define tip cells
message(paste0(Sys.time(), ": Defining tip-cells"))
clustering <- "Louvain-20"
tips <- setdiff(unique(obj@group.ids[, clustering]), NA)
this.tip <- tips[tip.to.walk]

#Simulate the biased random walks from each tip
message(paste0(Sys.time(), ": Simulating random walks from each tip"))
tip.walks <- simulateRandomWalksFromTips(obj, tip.group.id = "tip.clusters", root.cells = root.cells,
                                         transition.matrix = biased.tm, n.per.tip = 25000, root.visits = 1,
                                         max.steps = 5000, verbose = T)

##Alternatively, the random walks can be looped over individual clusters
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
names(walks) <- rownames(table(obj@group.ids$`Cluster-Num`))

saveRDS(walks, file = paste0(save.path, sample, "_walks.rds"))
saveRDS(tip.walks, file = paste0(save.path, "/", folder_name, "/", dataset, "_walks_more_var_genes.rds"))
#Process the biased random walks into visitation frequencies
message(paste0(Sys.time(), ": Processing the biased random walks"))
obj <- processRandomWalksFromTips(obj, walks, verbose = T)

#Visualize visitation of cells from each tip
plotDim(obj, "tip.clusters", plot.title = "Cells in each tip")

plotDim(obj, "visitfreq.log.1", plot.title = "Visitation frequency from tip-1 (log10)", transitions.plot = 10000)
plotDim(obj, "visitfreq.log.5", plot.title = "Visitation frequency from tip-2 (log10)", transitions.plot = 10000)

##Save the object
saveRDS(obj, file=paste0(save.path, sample, "_URD_withWalks.rds"))
saveRDS(obj, file=paste0(save.path, "/", folder_name, "/", dataset, "_URD_withWalks_withTuftCells.rds"))
obj <- readRDS(file=paste0(save.path, sample, "_URD_withWalks.rds"))
obj <- readRDS(file=paste0(save.path, sample, "_URD_withWalks_withTuftCells.rds"))


################################# BUILDING  THE URD TREE #####################################################
#PART 5 - Building URD Tree
library(URD)
library(rgl)

#Set up knitr to capture rgl output
rgl::setupKnitr()

#Load tip cells
message(paste0(Sys.time(), ": Loading tip cells"))
obj <- loadTipCells(obj, tips = "Cluster-Num")

#Build the actual tree
message(paste0(Sys.time(), ": Decide on tips to use for tree construction"))
tips.to.exclude <- c("3")
tips.to.use <- setdiff(as.character(1:9), tips.to.exclude)
tips.to.use <- c("1", "2", "3", "4", "5", "6", "7")

#Build the tree
##obj.built <- buildTree(object = obj, pseudotime = "pseudotime", divergence.method = "ks",
                     #  tips.use = tips.to.use, weighted.fusion = T, use.only.original.tips = T,
                      # cells.per.pseudotime.bin = 25, bins.per.pseudotime.window = 5, minimum.visits = 1,
                       #visit.threshold = 0.7, p.thresh = 0.001, save.breakpoint.plots = NULL, dendro.node.size = 100,
                       #min.cells.per.segment = 10, min.pseudotime.per.segment = 0.1, verbose = F)

obj.tree <- buildTree(obj, pseudotime = "pseudotime", tips.use = tips.to.use, divergence.method = "preference", cells.per.pseudotime.bin = 25,
                          bins.per.pseudotime.window = 5, save.all.breakpoint.info = T, p.thresh = 0.1, verbose = F)

#Name the tips
tip.names <- unique(obj@group.ids[, c("Cluster", "Cluster-Num")])
tip.names <- tip.names[complete.cases(tip.names), ]
obj.tree <- nameSegments(obj.tree, segments = tip.names$`Cluster-Num`, segment.names = tip.names$`Cluster`)
obj.built <- nameSegments(obj.built, segments = tip.names$`Cluster-Num`, segment.names = tip.names$`Cluster`)

##Finally, plot the tree
plotTree(obj.tree, "stage.group", label.segments = T, title = "Developmental_Stage_tree")
plotTree(obj.built, "stage.group", label.segments = T, title = "Developmental_Stage_tree")

genes.plot <- c("tnfrsf11a", "atoh1b", "ascl1a", "sox4b")
for (gene in genes.plot) {
  plot(plotTree(obj.tree, gene))
}

##Plot gene expression on the URD trajectory
gridExtra::grid.arrange(grobs = lapply(c("pou2f3", "sox8b", "tnfrsf11a", "atoh1b"), plotTree,
                                       object = obj.tree, label.x = F, plot.cells = T), ncol = 2)

gridExtra::grid.arrange(grobs = lapply(c("cdx1b", "ascl1a", "sox4a", "fabp1b.1"), plotTree,
                                       object = obj.tree, label.x = F, plot.cells = T), ncol = 2)

gridExtra::grid.arrange(grobs = lapply(c("blnk", "ehf", "cobll1a", "sall1a"), plotTree,
                                       object = obj.tree, label.x = F, plot.cells = T), ncol = 2)


##Check where all the cells from the Seurat object lie in the URD tree
#cells.15 <- WhichCells(obj1, idents = c("15"))

##Save the tree object
saveRDS(obj.tree, file=paste0(save.path, sample, "_TREE.rds"))
obj.tree <- readRDS(file=paste0(save.path, sample, "_TREE.rds"))

#Manual refinement - rename the segments manually
message(paste0(Sys.time(), ": Descriptive names to be used on the dendrogram"))
new.seg.names <- c("EC-1", "EC-2", "EC-3", "LREs", "EECs", "goblet", "B4O2", "tuft-like")
segs.to.name <- c("2", "6", "1", "4", "7", "8", "5", "9")
obj.tree <- nameSegments(obj.tree, segments = segs.to.name, segment.names = new.seg.names)

plotTree(obj.tree, "stage.group", title = "Developmental_stage_tree", label.segments = T)

diff.exp <- function(segment.1, segment.2) {
  markers.comp <- markersAUCPR(obj.tree, cells.1 = obj.tree@tree[["cells.in.segment"]][[segment.1]], cells.2 = obj.tree@tree[["cells.in.segment"]][[segment.2]], effect.size = 0.5, frac.min.diff = 0.1)
  
write.csv(markers.comp, file = paste0("~/Box/zfext/results/06a-URD/markers/Hematopoietic_superset_postprocess/", dataset, "_markers_", segment.1, "vs", segment.2, ".csv"))
}
diff.exp("2", "1")







