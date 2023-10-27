##This code was used to generate the shared gene module heatmap shown in Figure 2A and the subsequent dotplots shown in Figure 2B-D and Figure S4.
##Once the gene expression programs were calculated on the global dataset, they were then resolved onto the clustering shown in Figure 1B. 

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

# Load cell annotations to introduce for downstream purposes
cell.annot.path <- "~/Box/zfext/annotations_celltype_curated_newMama/celltype_annotations_df_full.tsv"
message(Sys.time(), ": Loading ", cell.annot.path)
cell.annot <- read.table(cell.annot.path, sep = " ", header = T, stringsAsFactors = F)
cell.annot$stage.group <- plyr::mapvalues(
  from = c("  3-4", "  5-6", "  7-9", " 10-12", " 14-21", " 24-34", " 36-46", 
           " 48-58", " 60-70", " 72-82", " 84-94", " 96-106", " 108-118", 
           "120", "5-6", "7-9", "10-12", "14-21", "24-34", "36-46", "48-58", 
           "60-70", "72-82", "84-94", "96-106", "108-118"),
  to =   c("  3-4", "  5-6", "  7-9", " 10-12", " 14-21", " 24-34", " 36-46", 
           " 48-58", " 60-70", " 72-82", " 84-94", " 96-106", " 108-118", 
           "120", "  5-6", "  7-9", " 10-12", " 14-21", " 24-34", " 36-46", " 48-58", 
           " 60-70", " 72-82", " 84-94", " 96-106", " 108-118"),
  x = cell.annot$stage.group
)
cell.annot$tissue <- unlist(lapply(strsplit(x = cell.annot$clust, split = "\\."), function(x) x[1]))
cell.annot$clust.sg <- paste0(cell.annot$clust, "_", gsub(" ", "", cell.annot$stage.group))
cell.annot$clust.stage <- paste0(cell.annot$clust, "_", gsub(" ", "", cell.annot$stage.nice))
cell.annot$tissue.sg <- paste0(cell.annot$tissue, "_", gsub(" ", "", cell.annot$stage.group))
cell.annot$tissue.stage <- paste0(cell.annot$tissue, "_", gsub(" ", "", cell.annot$stage.nice))

mama@meta.data <- cbind(
  mama@meta.data[,setdiff(colnames(mama@meta.data), colnames(cell.annot))],
  cell.annot[rownames(mama@meta.data),]
)


## CREATE ADDITIONAL META DATA IN MAMA --------------------

# Create any additional needed meta.data
mama@meta.data$hpf <- as.numeric(mama@meta.data$stage.nice)

##Replace all underscores in the tissue.subset slot
mama@meta.data$tissue.annot <- plyr::mapvalues(
  from = c("blastula", "periderm", "axial_mesoderm", "gastrula", "neural", "hematopoietic", "muscle_superset", "endoderm", "cephalic_mesoderm",
           "PGCs", "pronephros", "unknown", "mesenchyme", "basal_epidermis", "ionocytes_mucous-secreting", "taste_olfactory", "glial_cells", "eye",
           "otic", "pigment-cells", "fin", "non-skeletal_muscle"),
  to =   c("blastula", "periderm", "axial", "gastrula", "neural", "hematopoietic", "muscle", "endoderm", "cephalic",
           "PGCs", "pronephros", "unknown", "mesenchyme", "epidermis", "ionocytes", "taste", "glial", "eye",
           "otic", "pigment", "fin", "mural"),
  x = mama@meta.data$tissue.subsets
)

# Create stage.group / cluster
mama@meta.data$tissue.group <- URD::cutstring(delimiter = "-", field = 1, x = 
                                                gsub(" ", "", 
                                                     apply(mama@meta.data[,c("tissue.annot", "stage.group")], 1, paste0, collapse = "_")
                                                )
)

mama@meta.data$tissue.hpf <- gsub(" ", "", apply(mama@meta.data[,c("tissue.annot", "hpf")], 1, paste0, collapse = "_"))
mama@meta.data$tissue_stage.diff <- paste0(mama@meta.data$clust, "_", gsub(" ", "", mama@meta.data$stage.diff))

##Now load the cell-embeddings calculated on mama
all.emb <- mama.full.red@cell.embeddings
all.emb <- as.data.frame(all.emb)
all.emb$tissue.clust <- mama@meta.data$tissue.group

# Load cluster annotations 
is.empty <- function(x) {
  is.na(x) | is.null(x) | (sapply(trimws(x), nchar) == 0)
}
annot.confirmed <- readxl::read_xlsx("~/Box/zfext/02-Clustering/mama_tissue_clustering/annot/mama_annotations_full_extended_jeffedit_newtissues.xlsx")
annot.confirmed <- annot.confirmed[!is.empty(annot.confirmed$clust),3:ncol(annot.confirmed)]
rownames(annot.confirmed) <- annot.confirmed$daniocell.cluster
annot.confirmed$confirmed <- T

##Add annotations to mama metadata
mama@meta.data$tissue.annot <- NA
for(i in annot.confirmed$clust){
  mama@meta.data[which(mama@meta.data$clust == i), "tissue.annot"] <- annot.confirmed[which(annot.confirmed$clust == i), "identity.super"]
}


##Load curated clustering that is shown in Figure 1C
curated.clusters <- readRDS("~/Box/zfext/02-Clustering/mama+ds_integrated/tissue_curated_clusters_mama_figure.rds")
mama@meta.data$curated.clusters <- curated.clusters
DimPlot(mama, group.by = "curated.clusters")


##Resolve the module calculations to the various tissue and cluster resolutions calculated on the global dataset
module.exp.per.group <- aggregate(mama.full.red@cell.embeddings, by = list(mama@meta.data[rownames(mama.full.red@cell.embeddings), "tissue.group"]), FUN = mean)
module.exp.per.clust <- aggregate(mama.full.red@cell.embeddings, by = list(mama@meta.data[rownames(mama.full.red@cell.embeddings), "clust"]), FUN = mean)
rownames(module.exp.per.group) <- module.exp.per.group$Group.1
rownames(module.exp.per.clust) <- module.exp.per.clust$Group.1
module.exp.per.group.log <- log2(module.exp.per.group[,2:ncol(module.exp.per.group)] + 1)
module.exp.per.clust.log <- log2(module.exp.per.clust[,2:ncol(module.exp.per.clust)] + 1)
module.exp.per.group.log$tissue <- unlist(lapply(strsplit(rownames(module.exp.per.group.log), ".", fixed = T), function(x) x[1]))
module.exp.per.clust.log$clust <- unlist(lapply(strsplit(rownames(module.exp.per.clust.log), ".", fixed = T), function(x) x[2]))
module.exp.per.clust.log[is.na(module.exp.per.clust.log$clust), "clust"] <- 1

module.exp.per.annot <- aggregate(mama.full.red@cell.embeddings, by = list(mama@meta.data[rownames(mama.full.red@cell.embeddings), "identity.super"]), FUN = mean)
rownames(module.exp.per.annot) <- module.exp.per.annot$Group.1
module.exp.per.annot.log <- log2(module.exp.per.annot[,2:ncol(module.exp.per.annot)] + 1)

##Calculate module expression at the level of curated.clustering
module.exp.per.curated.cluster <- aggregate(mama.full.red@cell.embeddings, by = list(mama@meta.data[rownames(mama.full.red@cell.embeddings), "curated.clusters"]), FUN = mean)
rownames(module.exp.per.curated.cluster) <- module.exp.per.curated.cluster$Group.1
module.exp.per.curated.cluster.log <- log2(module.exp.per.curated.cluster[,2:ncol(module.exp.per.curated.cluster)] + 1)

##Calculate module expression at the level of Figure 1C - the identity.super catgories
module.exp.per.identity.super <- aggregate(mama.full.red@cell.embeddings, by = list(mama@meta.data[rownames(mama.full.red@cell.embeddings), "tissue.annot"]), FUN = mean)
rownames(module.exp.per.identity.super) <- module.exp.per.identity.super$Group.1
module.exp.per.identity.super.log <- log2(module.exp.per.identity.super[,2:ncol(module.exp.per.identity.super)] + 1)

##Make a heatmap for tissue vs modules
module.exp.per.tissue <- aggregate(mama.full.red@cell.embeddings, by = list(mama@meta.data[rownames(mama.full.red@cell.embeddings), "tissue"]), FUN = mean)
rownames(module.exp.per.tissue) <- module.exp.per.tissue$Group.1
module.exp.per.tissue.log <- log2(module.exp.per.tissue[,2:ncol(module.exp.per.tissue)] + 1)

##Saving the modules calculated on various categories (including log calculations)
saveRDS(module.exp.per.clust, "~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/data/module_expression_per_clust.rds")  
saveRDS(module.exp.per.clust.log, "~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/data/module_expression_per_clust_log.rds")  
saveRDS(module.exp.per.tissue, "~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/data/module_expression_per_tissue.rds")  
saveRDS(module.exp.per.tissue.log, "~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/data/module_expression_per_tissue_log.rds") 
saveRDS(module.exp.per.group, "~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/data/module_expression_per_tissue_stg_group.rds")
saveRDS(module.exp.per.annot, "~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/data/module_expression_per_identity.super.rds")
saveRDS(module.exp.per.identity.super, "~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/data/module_expression_per_identity_super_new_annot.rds")
saveRDS(module.exp.per.identity.super.log, "~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/data/module_expression_per_identity_super_new_annot_log.rds")

module.exp.per.tissue <- readRDS("~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/data/module_expression_per_tissue.rds")
module.exp.per.tissue.log <- readRDS("~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/data/module_expression_per_tissue_log.rds")
module.exp.per.annot <- readRDS("~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/data/module_expression_per_identity.super.rds")

##Plot a heatmap of expression of these modules across the tissue_stage.groups
library(ComplexHeatmap)
library(NMF)
library(viridis)
library(pheatmap)

##Generate a binary heatmap for the clustering shown in Figure 1C and their expression of shared modules - shown in Figure 3A
##Have the full identity.super shared module heatmap 

##Load relevant module calculations
module.exp.per.annot <- readRDS("~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/data/module_expression_per_identity.super.rds")
module.exp.per.annot.log <- log2(module.exp.per.annot[,2:ncol(module.exp.per.annot)] + 1)

##Subset matrix to have only shared modules
mat <- module.exp.per.annot.log[, shared.fc.num]

##Make a binary matrix
mat.norm <- sweep(mat, 2, apply(mat, 2, max), "/")
mat.norm.bin <- matrix(0, nrow = nrow(mat.norm), ncol = ncol(mat.norm))
mat.norm.bin[mat.norm > 0.5] <- 1

##Use the same column names
colnames(mat.norm.bin) <- colnames(mat.norm)
rownames(mat.norm.bin) <- rownames(mat.norm)

##Subset this matrix such that only columns that are present in more than one rows is kept
mat.norm.bin.filtered <- mat.norm.bin[, colSums(mat.norm.bin) > 1]

##Remove the module 10 that is a weird cluster shared with all cells
mat.norm.bin <- mat.norm.bin[, -7]

##Define a function to calculate Z-score
cal_z_score <- function(x){
  (x-mean(x))/sd(x)
}

data_subset_norm <- t(apply(mat.norm.bin.filtered, 1, cal_z_score))
pheatmap(data_subset_norm, annotation_colors = df$n.genes,
         cutree_cols = 4,
         cutree_rows = 4,
         legend = T)

##Plot a heatmap with modules vs tissue.groups - Figure 3A
library(pheatmap)
pdf(file = "~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/figures/module_expression_across_identity.super_shared_new_annot_v2.pdf", width = 72, height = 60)
pheatmap(mat.norm.bin.filtered, 
         annotation_colors = df,
         cutree_cols = 4,
         cutree_rows = 4,
         color = c('#e0e0e0','#f5f5f5','#5ab4ac'), legend = T)
dev.off()


##Find the number of genes in each module
genes.in.module <- names(which(sort(mama.full.red@feature.loadings, decreasing = T) > 0))
##Create a dataframe with the module names on a column
df <- as.data.frame(colnames(mama.full.red@feature.loadings))
colnames(df) <- "modules"

n.genes <- lapply(df$modules, function (i) length(names(which(sort(mama.full.red@feature.loadings[, i], decreasing = T) > 0))))
genes.each.module <- unlist(n.genes)        

df$n.genes <- genes.each.module
df <- df[which(df$modules %in% colnames(mat.norm.bin)), ]

##Find number of genes present in more than one module
genes.each.module <- lapply(colnames(mama.full.red@feature.loadings), function (i) names(which(sort(mama.full.red@feature.loadings[, i], decreasing = T) > 0)))
names(genes.each.module) <- colnames(mama.full.red@feature.loadings)

##Find all unique genes
genes.unique <- unlist(unique(genes.each.module))
genes.all <- unlist(genes.each.module)


########################## Plotting individual genes expressed within shared gene modules shared across 2 or more tissue ########################
##Plot the genes expressed in the epithelial matrix
epidermis.fc <-  mod.annot$Module[which(mod.annot$`Module description` == "epithelial")]
grep("epithelia", mod.annot$`Module description`, value = T)

##Create a dotplot for the ciliary genes - Figure S3A
##Ok get all the cell-types that shares cilia genes and a cluster that is not ciliated
##Clusters with ciliary genes - pron.11, pron.16, tast.17, otic.14, otic.17, otic.20, glia.16, glia.19
obj <- subset(mama, idents = c("pron.11", "pron.16", "tast.18", "otic.14", "otic.17", "otic.20", "glia.15", "glia.26", "endo.21"))
Idents(obj) <- factor(levels(Idents(obj)), levels = c("pron.11", "pron.16", "tast.18", "otic.14", "otic.17", "otic.20", "glia.15", "glia.26", "endo.21"))
genes.plot <- c("ccdc173", "ccdc169", "ccdc151", "ccdc114", "lrrc6", "dnah5", "dnah11", "dnah6", "dnah12", "dnaaf1", "dnali1", "dydc2", "dnai1.2", "spata17", "armc4", "spag1a", "spaca9", "cfap100", "cfap99", "cfap43", "cfap52", "tbata", "daw1", "foxj1a")
pdf(file = "~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/figures/ciliary_genes_ODA.pdf", width = 16, height = 20)
DotPlot(mama, features = genes.plot, idents = c("pron.11", "pron.16", "tast.18", "otic.14", "otic.17", "otic.20", "glia.15", "glia.26", "endo.21"), dot.scale = 18, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

genes.rsp <- c("tekt1", "tekt4", "rsph14", "rsph9", "rsph4a", "spag6", "spag8", "hydin", "capsla", "capslb", "lrrc23", "dnajb13", "enkur", "foxj1b")
pdf(file = "~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/figures/ciliary_genes_RSP.pdf", width = 16, height = 20)
DotPlot(mama, features = genes.rsp, idents = c("pron.11", "pron.16", "tast.18", "otic.14", "otic.17", "otic.20", "glia.15", "glia.26", "endo.21"), dot.scale = 18, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

genes.ift <- c("ift57", "ift74", "ift22", "ift46", "ift27", "ift52", "ift80", "dynlt1", "dync2li1", "dpcd", "dnal1", "cfap20", "ttc26", "rfx2")
pdf(file = "~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/figures/ciliary_genes_IFT.pdf", width = 16, height = 20)
DotPlot(mama, features = genes.ift, idents = c("pron.11", "pron.16", "tast.18", "otic.14", "otic.17", "otic.20", "glia.15", "glia.26", "endo.21"), dot.scale = 18, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

genes.mcc <- c("odf3b", "ccdc33", "crocc2", "dnah5l", "dnah1", "dnah2", "dnah6", "dnah7", "dnah3", "tekt2", "zmynd12", "spata18", "tex36", "lrrc18b", "gmnc")
pdf(file = "~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/figures/ciliary_genes_MCC.pdf", width = 16, height = 20)
DotPlot(mama, features = genes.mcc, idents = c("pron.11", "pron.16", "tast.18", "otic.14", "otic.17", "otic.20", "glia.15", "glia.26", "endo.21"), dot.scale = 18, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()


##Load curated clustering
curated.clusters <- readRDS("~/Box/zfext/02-Clustering/mama+ds_integrated/tissue_curated_clusters_2.rds")
mama@meta.data$curated.clusters <- curated.clusters
DimPlot(mama, group.by = "curated.clusters")

##For epithelial dotplot, we need to assign different idents to mama - Figure 3D
Idents(mama) <- mama@meta.data$curated.clusters
mama@meta.data[which(is.na(mama@meta.data$curated.clusters)), "curated.clusters"] <- "unknown"
##Order the idents
#Idents(mama) <- factor(levels(Idents(mama)), levels = c("EVL/periderm", "skin/epidermis", "ionocytes", "olfactory epithelium", "taste epithelium", "RPE", "otic/lateral-line", "lens" , "hatching gland", "mucous-secreting", "fin-epithelia", "fin", "pronephros", "liver", "gut", "muscle"))

epithelial.genes <- c("nfe2l2a", "mhc1zba", "irf2", "tmem59", "b2ml", "litaf", "enosf1", "wbp2nl", "tuft1a", "cish", "prdx5", "adipor2", "chmp5b")
pdf(file = "~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/figures/epithelial_module_shared_all_1_v2.pdf", width = 30, height = 20)
DotPlot(mama, assay = "RNA", features = epithelial.genes, idents = c("EVL/periderm", "skin/epidermis", "ionocytes", "olfactory epithelium", "taste epithelium", "RPE", "otic/lateral-line", "lens" , "hatching gland", "mucous-secreting", "fin-epithelia", "pronephros", "liver", "gut", "muscle"), dot.scale = 18, scale = F) + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()


epithelial.genes.1 <- c("gbgt1l4", "oclnb", "cldnb", "cldnf", "spint2", "chmp4c", "tmem30b", "gmds", "ildr1a", "tnks1bp1", "pkp1b", "macc1", "sh3yl1", "esrp1", "nectin4a")
pdf(file = "~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/figures/epithelial_module_shared_1_v2.pdf", width = 30, height = 20)
DotPlot(mama, assay = "RNA", features = epithelial.genes.1, idents = c("EVL/periderm", "skin/epidermis", "ionocytes", "olfactory epithelium", "taste epithelium", "RPE", "otic/lateral-line", "lens" , "hatching gland", "mucous-secreting", "fin-epithelia", "pronephros", "liver", "gut", "muscle"), dot.scale = 18, scale = F) + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()


epithelial.genes.2 <- c("cfl1l", "tmsb1", "capns1a", "oclna", "krt5", "krtt1c19e", "epcam", "anxa1a", "anxa2a", "arhgef5", "spint1a", "spint1b", "pycard", "foxq1a")
pdf(file = "~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/figures/epithelial_module_shared_2_v2.pdf", width = 30, height = 20)
DotPlot(mama, assay = "RNA", features = epithelial.genes.2, idents = c("EVL/periderm", "skin/epidermis", "ionocytes", "olfactory epithelium", "taste epithelium", "RPE", "otic/lateral-line", "lens" , "hatching gland", "mucous-secreting", "fin-epithelia", "pronephros", "liver", "gut", "muscle"), dot.scale = 18, scale = F) + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()


epithelial.genes.3 <- c("anxa1c", "anxa1b", "anxa3a", "lye", "agr1", "stx11b.1", "icn2", "sft2d2a", "elovl7a", "elovl7b", "elovl6l", "cldne", "capn9", "znf185")
pdf(file = "~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/figures/epithelial_module_shared_3_v2.pdf", width = 30, height = 20)
DotPlot(mama, assay = "RNA", features = epithelial.genes.3, idents = c("EVL/periderm", "skin/epidermis", "ionocytes", "olfactory epithelium", "taste epithelium", "RPE", "otic/lateral-line", "lens" , "hatching gland", "mucous-secreting", "fin-epithelia", "pronephros", "liver", "gut", "muscle"), dot.scale = 18, scale = F) + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()


epithelial.genes.4 <- c("cldn15la", "cldn15lb", "cldn15a", "cldnc", "anxa2b", "krt18a.2", "cdhr5b", "gpa33b", "anks4b", "eps8l3b", "pllp", "itga6l")
pdf(file = "~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/figures/epithelial_module_shared_4_v2.pdf", width = 30, height = 20)
DotPlot(mama, assay = "RNA", features = epithelial.genes.4, idents = c("EVL/periderm", "skin/epidermis", "ionocytes", "olfactory epithelium", "taste epithelium", "RPE", "otic/lateral-line", "lens" , "hatching gland", "mucous-secreting", "fin-epithelia", "pronephros", "liver", "gut", "muscle"), dot.scale = 18, scale = F) + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()


##Let's do the a figure for the lysosomal degradation pathway - Figure S3B
genes.lyz <- c("slc10a2", "ctsa", "ctsc", "ctsh", "ctsz", "ctsba", "lgmn", "dpp4", "dpp7", "hexa", "hexb", "man2b1", "man2b2", "manba", "pepd", "scpep1", "cpvl")
Idents(mama) <- mama@meta.data$clust
pdf(file = "~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/figures/lysosome_degradation_module_all.pdf", width = 30, height = 20)
DotPlot(mama, assay = "RNA", features = genes.lyz, idents = sort(c("pron.1", "pron.3", "pron.5", "pron.7", "pron.14", "endo.11",
                                                                   "hema.12", "hema.26", "hema.37", "pigm.14", "neur.14")), dot.scale = 18, scale = F) + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

##Also plot the members of the SLC pathway and the megalin expression pathway - Figure 3B
genes.slc <- c("cubn", "dab2", "amn", "lrp2a", "lrp2b", "slc5a11", "slc22a7b.3", "slc16a13", "slc6a19a.2", "slc6a19b", "slc5a12", "slc5a2", "slc22a4", "slc22a13b", "slc20a1a", "slc34a1a")
Idents(mama) <- mama@meta.data$clust
pdf(file = "~/Box/zfext/genemodule_celltype_curated/2022-05 modules on mama/figures/SLC_transporter_module_2.pdf", width = 30, height = 20)
DotPlot(mama, assay = "RNA", features = genes.slc, idents = sort(c("pron.1", "pron.3", "pron.5", "pron.7", "pron.14", "endo.13")), dot.scale = 18, scale = F) + coord_flip() + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

##Plot a dotplot for the mucin genes and the associated O-glycosylation enzymes  - Figure 3C
sample <- "ionocytes"
obj.iono <- readRDS(paste0(file = "~/Box/zfext/annotations_celltype_curated_newMama/", sample, "/obj_seurat/", sample, "_seurat.rds"))

sample <- "endoderm"
obj.endo <- readRDS(paste0(file = "~/Box/zfext/annotations_celltype_curated_newMama/", sample, "/obj_seurat/", sample, "_seurat.rds"))

DimPlot(obj.iono, label = T)
cells.iono.skin <- WhichCells(obj.iono, idents = c("1", "7", "11", "14", "18"))
cells.iono.muc <- unlist(unique(list(cells.iono.gill, cells.iono.skin)))
obj.iono.skin <- subset(obj.iono, cells = cells.iono.skin)

DimPlot(obj.endo, label = T)
cells.esophagus <- WhichCells(obj.endo, idents = c("30"))
cells.goblet <- WhichCells(obj.endo, idents = c("15"))
cells.endo <- unlist(unique(list(cells.esophagus, cells.goblet)))

obj.endo <- subset(obj.endo, cells = cells.endo)

obj.muc <- merge(obj.endo, obj.iono.skin)

obj.muc@meta.data$ident <- NA
obj.muc@meta.data[cells.iono.skin, "ident"] <- "skin"
obj.muc@meta.data[cells.esophagus, "ident"] <- "esophagus"
obj.muc@meta.data[cells.goblet, "ident"] <- "goblet-cells"

Idents(obj.muc) <- obj.muc@meta.data$ident


library(viridis)
pdf("mucin_dotplot.pdf", width = 8, height = 6)
genes.list <- c("muc5.1", "muc5.2", "muc5.3", "muc5f", "muc2.1", "muc13b")
DotPlot(obj.muc, features = genes.list, dot.scale = 16, scale = F) + coord_flip() + 
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

pdf("mucin_GALNT_dotplot.pdf", width = 8, height = 6)
genes.list <- c("galnt8b.1", "galnt12", "galnt6", "galnt7", "galnt2")
DotPlot(obj.muc, features = genes.list, dot.scale = 16, scale = F) + coord_flip() + 
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

pdf("mucin_B4GALNT_dotplot.pdf", width = 8, height = 6)
genes.list <- c("b3gnt7", "b3gnt3.4", "b3gnt5a", "b3gnt5b", "b4galnt3a", "gcnt4a")
DotPlot(obj.muc, features = genes.list, dot.scale = 16, scale = F) + coord_flip() + 
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

pdf("mucin_ST6galnac_dotplot.pdf", width = 8, height = 6)
genes.list <- c("st6galnac1.1", "st6gal2b", "st3gal7", "st3gal1l2", "st3gal1", "gal3st3")
DotPlot(obj.muc, features = genes.list, dot.scale = 16, scale = F) + coord_flip() + 
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

genes.mucus <- c("muc5.1", "muc5.2", "muc5.3", "muc5f", "muc2.1", "muc13b", "galnt8b.1", "galnt12", "galnt6", "galnt7", "galnt2", 
                 "b3gnt7", "b3gnt3.4", "b3gnt5a", "b3gnt5b", "b4galnt3a", "gcnt4a", "st6galnac1.1", "st6gal2b", "st3gal7", "st3gal1l2", "st3gal1", "gal3st3")
pdf("mucin_enzyme_diff_exp.pdf", width = 8, height = 6)
DotPlot(obj.muc, features = genes.mucus, dot.scale = 10, scale = F) + coord_flip() + 
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()


