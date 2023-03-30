##Load libraries
library(Seurat)
library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(ggplot2)

#' Safely expand df or matrix to match rows
#' 
#' @param x (Matrix or data.frame)
#' @param margin (Numeric) Match Rows (1) or Columns (2)
#' @param match (Character) Rows/Columns that should be present in final object
#' @param fill (What to fill missing rows with, default is NA)
#' @return Matrix or data.frame, with rows/columns in same order as match
expand.matrix <- function(x, margin, match, fill = NA) {
  if (margin == 1) { # ROWS
    rows.keep <- intersect(rownames(x), match)
    rows.missing <- setdiff(match, rows.keep)
    add.matrix <- matrix(rep(fill, length(rows.missing) * ncol(x)), ncol = ncol(x))
    rownames(add.matrix) <- rows.missing
    colnames(add.matrix) <- colnames(x)
    y <- rbind(x, add.matrix)
    y <- y[match,]
    return(y)
  } else if (margin == 2) { # COLUMNS
    cols.keep <- intersect(colnames(x), match)
    cols.missing <- setdiff(match, cols.keep)
    add.matrix <- matrix(rep(NA, length(cols.missing) * nrow(x)), nrow = nrow(x))
    colnames(add.matrix) <- cols.missing
    rownames(add.matrix) <- rownames(x)
    y <- cbind(x, add.matrix)
    y <- y[,match]
    return(y)
  } else {
    stop ("Margin must be 1 (rows) or 2 (columns)")
  }
}

cosine.similarity <- function(mat, by.col = T, dissim = F) {
  if (!("matrix" %in% class(mat))) stop("mat must be a matrix.")
  if (by.col) mat <- t(mat)
  sim <- mat / sqrt(rowSums(mat * mat))
  sim <- sim %*% t(sim)
  if (dissim) return(as.dist(1 - sim)) else return(sim)
}

# Load global dataset with original data
message(Sys.time(), ": Loading ", mama.path)
mama <- readRDS("~/Box/zfext/02-Clustering/mama+ds_integrated/merged_mama_ds_slim_seurat_25Jan2022.rds")

# Load tissue-specific modules
base.path <- "~/Box/zfext/genemodule_celltype_curated/fc_tissue_20220526/"
path.df <- data.frame(
  tissues = c("axial_mesoderm", "basal_epidermis", "blastula", "endoderm", "eye", "fin", "gastrula", "glial_cells", "hematopoietic", "ionocytes_mucous-secreting", "mural-cells", "muscle_superset", "otic", "periderm", "pgc", "pigment-cells", "pronephros", "taste_olfactory"),
  tissue.short = c("axia", "epid", "blas", "endo", "eye", "fin", "gast", "glia", "hema", "iono", "mura", "musc", "otic", "peri", "pgc", "pigm", "pron", "tast"),
  c.num = c(25, 39, 7, 30, 30, 22, 8, 33, 26, 21, 24, 36, 22, 32, 6, 23, 17, 18),
  stringsAsFactors = F
)
paths <- paste0(base.path, path.df$tissues, "_plots/c", path.df$c.num, "_m1.04_smooth/c", path.df$c.num, "_m1.04_smooth_mamaFCreduc.rds")
mama.mod.reduc <- lapply(paths, readRDS)
names(mama.mod.reduc) <- path.df$tissue.short

# Get feature loadings out
mod.fl <- lapply(path.df$tissue.short, function(tissue) {
  fl <- mama.mod.reduc[[tissue]]@feature.loadings
  colnames(fl) <- paste0(tissue, "_", colnames(fl))
  return(fl)
})
names(mod.fl) <- path.df$tissue.short

# Get all genes
all.mod.genes <- sort(unique(unlist(lapply(mod.fl, rownames))))
mod.fl.expand <- lapply(mod.fl, expand.matrix, margin = 1, match = all.mod.genes, fill = 0)
all.fl <- Reduce("cbind", mod.fl.expand)
range(rowSums(all.fl))

# Cluster modules? First try correlating them
mod.cor <- cor(all.fl)
heatmaply::heatmaply_cor(mod.cor, main = "Module Correlation", file = "~/Box/zfext/genemodule_celltype_curated/fc_tissue_20220526/20220812-module-cor-heatmap_with_mama.html")

mod.cossim <- cosine.similarity(as.matrix(all.fl), by.col = T)
heatmaply::heatmaply(mod.cossim, main = "Module Cosine Similarity", file = "~/Box/zfext/genemodule_celltype_curated/fc_tissue_20220526/20220812-module-cossim-heatmap_with_mama.html")

# Identify modules with a correlation to another module of at least 0.25
mod.cor.thresh <- mod.cor > 0.25
mod.with.match <- names(which(rowSums(mod.cor.thresh) > 1))
mod.cor.match <- mod.cor[mod.with.match, mod.with.match]

# Use hclust to turn this into clusters
mod.cor.dist <- as.dist(1-mod.cor.match)
mod.cor.hclust <- hclust(d = mod.cor.dist, method = "ward.D2")
table(cutree(mod.cor.hclust, h = 0.75))

##Add correlaion with the modules calculated on the global dataset to the previous correlation matrix
mama.full.red <- readRDS("~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/c200_m1.04_mamaFCreduc.rds")

##Ok, now the question is how to match the two sets of modules - modules calculated on tissues and modules calculated on mama
##So to do that, we need to correlate the tissue specific modules with the mama modules
mama.mod.reduc$mama.full <- mama.full.red

##Add the mama modules to the full list
mama.mod.reduc$mama.full <- mama.full.red

##Add the mama slot to the path.df
path.df <- data.frame(
  tissues = c("axial_mesoderm", "basal_epidermis", "blastula", "endoderm", "eye", "fin", "gastrula", "glial_cells", "hematopoietic", "ionocytes_mucous-secreting", "mural-cells", "muscle_superset", "otic", "periderm", "pgc", "pigment-cells", "pronephros", "taste_olfactory", 
              "mama.full"),
  tissue.short = c("axia", "epid", "blas", "endo", "eye", "fin", "gast", "glia", "hema", "iono", "mura", "musc", "otic", "peri", "pgc", "pigm", "pron", "tast", "mama.full"),
  c.num = c(25, 39, 7, 30, 30, 22, 8, 33, 26, 21, 24, 36, 22, 32, 6, 23, 17, 18, 147),
  stringsAsFactors = F
)


# Get feature loadings out
mod.fl <- lapply(path.df$tissue.short, function(tissue) {
  fl <- mama.mod.reduc[[tissue]]@feature.loadings
  colnames(fl) <- paste0(tissue, "_", colnames(fl))
  return(fl)
})
names(mod.fl) <- path.df$tissue.short

# Get all genes
all.mod.genes <- sort(unique(unlist(lapply(mod.fl, rownames))))
mod.fl.expand <- lapply(total.fl, expand.matrix, margin = 1, match = all.mod.genes, fill = 0)
all.fl <- Reduce("cbind", mod.fl.expand)
range(rowSums(all.fl))

##Using this list, I will correlate modules found in each tissue to that to mama
for(sample in path.df$tissue.short){
tissue.fl <- mama.mod.reduc[[sample]]@feature.loadings
mama.fl <- mama.mod.reduc$mama.full@feature.loadings
total.fl <- list(tissue.fl, mama.fl)

# Get all genes
all.mod.genes <- sort(unique(unlist(lapply(total.fl, rownames))))
mod.fl.expand <- lapply(total.fl, expand.matrix, margin = 1, match = all.mod.genes, fill = 0)
all.fl <- Reduce("cbind", mod.fl.expand)
range(rowSums(all.fl))

saveRDS(all.fl, file = paste0("~/Box/zfext/genemodule_celltype_curated/fc_tissue_20220526/tissue_mama_comparisons/20220812-", sample, "_mama_comparisons.rds"))

# Cluster modules? First try correlating them
mod.cor <- cor(all.fl)
heatmaply::heatmaply_cor(mod.cor, main = paste0(sample, "_module_correlation"), file = paste0("~/Box/zfext/genemodule_celltype_curated/fc_tissue_20220526/tissue_mama_comparisons/20220812-", sample, "_module-cor-heatmap_with_mama.html"))
}

mod.cor.thresh <- mod.cor > 0.25
mod.with.match <- names(which(rowSums(mod.cor.thresh) > 1))
mod.cor.match <- mod.cor[mod.with.match, mod.with.match]

##Ok, find shared genes between the tissue-specific modules and mama modules and see if all genes were represented in mama
tissue.mod.genes <- all.mod.genes
mama.mod.genes <- rownames(mama.full.red@feature.loadings)

genes.not.in.mama <- setdiff(tissue.mod.genes, mama.mod.genes)

#Check which modules have these genes not represented in mama
fl.diff <- all.fl[genes.not.in.mama, ]
fl.diff.good <- fl.diff > 0
fl.match <- names(which(colSums(fl.diff.good) > 1))

df.clean <- all.fl[fl.match, ]
df.good <- df.clean > 0
mod.names <- names(which(colSums(df.good) > 1))

##So these modules contain genes that were not included in the mama module calculations.
##What I need to figure out is that whether these genes represent important modules that were missed in the mama analysis

##Find multifunctional gene sets that are enriched in multiple modules calculated on mama
mama.fl <- mama.full.red@feature.loadings
all.mod.genes <- rownames(mama.fl)
mama.cor <- cor(mama.fl)
heatmap::heatmap(mama.cor, main = "mama_module_correlation", file = "~/Box/zfext/genemodule_celltype_curated/fc_tissue_20220526/20220812-mama_module-cor-heatmap_with_mama.html")





