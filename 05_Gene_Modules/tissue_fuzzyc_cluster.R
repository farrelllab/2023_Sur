#!/bin/Rscript

# tissue_fuzzyc_cluster.R
# Jeff Farrell - 04-19-2022
# Testing a parameter sweep of fuzzy-c clustering for finding gene modules

# Arguments:
#       1: c (Number of clusters to generate)
#       2: m (Fuzziness parameter for Mfuzz: 1 is hard clustering, >1 allows multiple membership)
#       3: do.smooth (Should k=5 gene expression smoothing be performed prior to clustering?)
#       4: tissue (Character)

library(Seurat)
library(Mfuzz)
library(Matrix)
#library(heatmaply)

source("/data/CSD/zfext/LTA/scripts/11-Modules/knn_smooth.R")
source("/data/CSD/zfext/LTA/scripts/11-Modules/2022-04-13 Gene Module Functions.R")

## Get arguments -----------------------------

message(paste0(Sys.time(), ": Starting."))
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) stop("Please provide 3 arguments. (See R script for details.)")

n <- as.numeric(args[1])
m <- as.numeric(args[2])
do.smooth <- as.logical(args[3])
tissue <- as.character(args[4])

# Load data --------------

# Load object
obj.path <- paste0("/data/CSD/zfext/results/04f-SubsetsV6/obj_subsets_v6/", tissue, "_seurat.rds")
message(Sys.time(), ": Loading", obj.path)
obj <- readRDS(obj.path)

# Load markers for endoderm object 
marker.path <- paste0("/data/CSD/zfext/results/04f-SubsetsV6/obj_subsets_v6/", tissue, "_markers.rds")
message(Sys.time(), ": Loading", marker.path)
markers <- readRDS(marker.path)

# Load GO annotations for scoring modules for parameter estimation
#zfin.go <- read.table("~/Downloads/zfin.gaf.gz", header = F, sep = "\t", comment.char = "!", quote = "")
#colnames(zfin.go) <- c("database", "marker", "symbol", "qualifier", "go.id", "ref.id", "go.evidence", "inferred", "ontology", "name", "synonyms", "type", "taxon", "mod.date", "assigned.by", "annotation.extension", "gene.product.form.id")

# Load mama for inspecting whether modules are expressed in other cell types
message(Sys.time(), ": Loading mama.")
mama <- readRDS("/data/CSD/zfext/LTA/results/02-Clustering/merged_timepoints/obj/merged_mama_ds_slim_seurat_25Jan2022.rds")

# Calculate modules ------------------

# Genes to operate on: all variable genes used in the object, plus the best markers of each cluster
markers.roc.top <- top.markers(markers$roc, min.diff = 0.47, min.power = 0.2, min.pct.diff = 0.1) # 1.6-fold / calibrated on best4.
genes.use <- sort(unique(c(obj@assays$RNA@var.features, markers.roc.top$gene)))

# Assemble expression set from regular or smoothed data.
if (do.smooth) {
  d.use <- max(obj@commands$FindNeighbors.RNA.pca$dims)
  if (!is.numeric(d.use)) {
    d.use <- 75
    message("Could not locate dimensions used in previous analysis, so going with default of 75.")
  }
  message("Smoothing with ", d.use, " PCs.")
  data.smoothed <- knn_smoothing(obj@assays$RNA@data[genes.use,], k = 5, d = d.use) # This uses KNN smoothing from Itai's lab
  exp.data <- ExpressionSet(as.matrix(data.smoothed))
  smooth.text <- "_smooth"
} else {
  exp.data <- ExpressionSet(as.matrix(obj@assays$RNA@data[genes.use,]))
  smooth.text <- ""
}
exp.data.std <- Mfuzz::standardise(exp.data)

# Generate c-means fuzzy clustering
cl <- mfuzz(exp.data.std, c=n, m=m, verbose = T)       # 3497 unassigned genes, 1147 single,  160 multi // 28 empty modules,  21 - 238 genes/module

# Filter clusters
# Add clustering result to Seurat
obj <- fuzzy.clust.to.seurat(obj, cl = cl, min.a = 0.2, min.gene.per.mod = 5, module.combine.thresh = 0.95, do.plot = F, add.to.meta = T)

# Project modules onto mama object to look in other cell types
mama <- project.reduction.between.seurat(seurat.source = obj, seurat.dest = mama, reduction = "fc", add.to.meta = T)

# Do plots!
file.prefix <- paste0("c", n, "_m", m, smooth.text)
plot.path <- paste0("/data/CSD/zfext/LTA/results/11-Modules/fc_tissue_20220526/", tissue, "_plots/", file.prefix, "/")

dir.create(plot.path, recursive = T)
write.csv(obj@reductions$fc@feature.loadings, file=paste0(plot.path, file.prefix, "_feature_loadings.csv"))
saveRDS(mama@reductions$fc, file = paste0(plot.path, file.prefix, "_mamaFCreduc.rds"))
saveRDS(obj@reductions$fc, file = paste0(plot.path, file.prefix, "_endoFCreduc.rds"))
plot.all.modules(obj = obj, mama = mama, plot.path = plot.path)


