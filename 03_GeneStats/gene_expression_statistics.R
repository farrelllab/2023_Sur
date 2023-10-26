library(Seurat)
library(ggplot2)

## GENERAL FUNCTIONS -------------

#' Cut string into fields
#'
#' Reworks \code{\link{strsplit}} to act like Unix cut. Instead of returning a list, returns a vector of the chosen field.
#' @param x (Character) A character vector
#' @param delimiter (Character) Delimiter character to use to split
#' @param field (Numeric) Which field to return; count starts at left with 1.
#' @return Character vector
cutstring <- function(x, delimiter, field) {
  unlist(lapply(strsplit(x, split = delimiter), function(y) y[field]))
}

is.empty <- function(x) {
  is.na(x) | is.null(x) | (sapply(trimws(x), nchar) == 0)
}

## LOAD DATA --------------
# Data and daniocell_load.R script available from Daniocell webiste (https://daniocell.nichd.nih.gov/)
source("daniocell_load.R")
daniocell@meta.data$celltype <- daniocell.annot[daniocell@meta.data$cluster, "clust.collapse"]
daniocell@meta.data$celltype.sg <- paste0(daniocell@meta.data$celltype, "_", gsub(" ", "", daniocell@meta.data$stage.group))

## GENERATE SUMMARY STATISTICS -------------------

# Load tables of expression data calculated per cluster/tissue/etc
num.expr.tissue <- readRDS("/Volumes/Daniocell/data-generated/mean-percent/num.exp.tissuegroup.rds")
mean.exp.stagegroup <- readRDS("/Volumes/Daniocell/data-generated/mean-percent/mean.exp.stagegroup.rds")
mean.exp.cluster <- readRDS("/Volumes/Daniocell/data-generated/mean-percent/mean.exp.cluster.rds")
mean.exp.celltype.sg <- readRDS("/Volumes/Daniocell/data-generated/mean-percent/mean.exp.celltype.sg.rds")

## EXPRESSED GENES ---------------
# Identify which genes are 'expressed' in enough cells to be considered 'captured'.

# 36250 genes in the annotation we used / in Seurat object.
num.expr <- rowSums(num.expr.tissue)
min.cells <- 22 # Number of cells in smallest cluster (thyroid)
genes.expr <- names(which(num.expr >= min.cells)) # 30928 genes expressed (>0 in 22 cells)
max.clust.mean <- apply(mean.exp.cluster, 1, max)
genes.expr.mean <- names(which(max.clust.mean >= 0.137)) # 24263 genes with mean of 0.1 counts/cell in at least 1 cluster
genes.expr.use <- intersect(genes.expr, genes.expr.mean) # 24216 genes (= 66.8% of genes in reference)

# Eliminate small categories 
# Limit categorizations to those categories with at least minimum number of cells
# This is important in cases where we've further split clusters by stage-group, etc.
# Using 22 cells as the threshold: smallest cluster (thyroid)
min.cells <- 22
keep.cat.min.n <- function(mean.matrix, meta.col, n, exclude = c("", " ", NA), verbose = T) {
  n.orig <- ncol(mean.matrix)
  n.cells <- table(daniocell@meta.data[,meta.col])
  pass.min <- names(which(n.cells >= n))
  pass.min <- setdiff(pass.min, exclude)
  mean.matrix <- mean.matrix[, intersect(colnames(mean.matrix), pass.min)]
  if (verbose) message("Retained ", ncol(mean.matrix), " of ", n.orig, " original categories.")
  return(mean.matrix)
}

mean.exp.cluster <- keep.cat.min.n(mean.exp.cluster, "cluster", n = min.cells)
mean.exp.celltype.sg <- keep.cat.min.n(mean.exp.celltype.sg, "celltype.sg", n = min.cells)

## DIVIDE GENES BY TYPE -------

# Define categories of genes
xloc.genes.all <- grep("^XLOC", rownames(daniocell@assays$RNA@data), value = T) # 1948 genes
loc.genes.all <- grep("^LOC", rownames(daniocell@assays$RNA@data), value = T) # 4779 genes
cdna.genes.all <- grep(":", rownames(daniocell@assays$RNA@data), value = T) # 4841 genes
contig.genes.all <- grep(
                  grep(".", rownames(daniocell@assays$RNA@data), fixed = T, value = T), 
                  pattern = "^AL|^CABZ|^BX|^CR|^CT|^CU|^CZ|^FO|^FP", ignore.case = F, value = T) # 4993 genes

xloc.genes.expr <- intersect(xloc.genes.all, genes.expr.use) # 723 genes
723/1984*100 # 37% expressed 
723/24216*100 # 3% of expressed genes
loc.genes.expr <- intersect(loc.genes.all, genes.expr.use) # 1159 genes (24%)
1159/4779*100 # 24% expressed
1159/24216*100 # 4.8% of expressed genes
cdna.genes.expr <- intersect(cdna.genes.all, genes.expr.use) # 3266 genes (67%)
3226/4841*100 # 67% expressed
3266/24216*100 # 13% of expressed genes
contig.genes.expr <- intersect(contig.genes.all, genes.expr.use) # 1618 genes (32%)
1618/4993*100 # 32% expressed
1618/24216*100 # 7% of expressed genes 
(24216-(723+1159+3266+1618))/24216*100 # Named genes are 72% of all expressed genes

# Calculate CV across populations -----------

# Calculate CV of mean expression across clusters, 'cell type' reduced clusters, and tissues
CV.exp.sg <- apply(mean.exp.stagegroup, 1, function(x) sd(x)/mean(x))
CV.exp.cluster <- apply(mean.exp.cluster, 1, function(x) sd(x)/mean(x))

# Calculate CV of mean expression across 'cell types' from each stage.group
all.stage.groups <- gsub(" ", "", sort(unique(daniocell@meta.data$stage.group)))
ct.per.sg <- lapply(all.stage.groups, function(sg) {
  return(colnames(mean.exp.celltype.sg)[which(cutstring(colnames(mean.exp.celltype.sg), "_", 2) == sg)])
})
names(ct.per.sg) <- all.stage.groups
CV.exp.ct.per.sg <- as.data.frame(lapply(ct.per.sg, function(clusts) {
  apply(mean.exp.celltype.sg[,clusts], 1, function(x) sd(x)/mean(x))
}))
colnames(CV.exp.ct.per.sg) <- all.stage.groups
is.exp.per.sg <- as.data.frame(lapply(ct.per.sg, function(clusts) {
  apply(mean.exp.celltype.sg[,clusts], 1, function(x) any(x > 0.1))
}))
colnames(is.exp.per.sg) <- all.stage.groups
CV.exp.ct.per.sg[!is.exp.per.sg] <- NA

# Plot relative distributions
df <- data.frame(
  gene = genes.expr.use,
  type = "Named",
  CV.exp.cluster = CV.exp.cluster[genes.expr.use],
  CV.exp.cluster.log = log(CV.exp.cluster[genes.expr.use]),
  CV.exp.sg = CV.exp.sg[genes.expr.use],
  CV.exp.sg.log = log(CV.exp.sg[genes.expr.use]),
  row.names = genes.expr.use,
  stringsAsFactors = F
)
# Add categories
df[xloc.genes.expr, "type"] <- "XLOC"        # "Lawson v4.3.2 (XLOC)"
df[loc.genes.expr, "type"] <- "Predicted"    # "Refseq Model (LOC)"
df[cdna.genes.expr, "type"] <- "cDNA"        # "cDNA clone (e.g. si:dkey/si:ch211)"
df[contig.genes.expr, "type"] <- "Predicted" # "Predicted (e.g. BX/CABZ)"

# Identify maximum variation among cell types PER stage.group -- (i.e. if time is held still, does the gene vary between cell types)
# Make NA values -Inf so they stop messing up max, but never get picked.
CV.exp.ct.per.sg.naRM <- CV.exp.ct.per.sg
CV.exp.ct.per.sg.naRM[is.na(CV.exp.ct.per.sg.naRM)] <- -Inf 
df$CV.exp.ct.per.sg <- apply(CV.exp.ct.per.sg.naRM[genes.expr.use,], 1, max)
df$CV.exp.ct.per.sg.log <- log(df$CV.exp.ct.per.sg)

# Calculate number of 'expressed' clusters
max.expression <- apply(mean.exp.cluster, 1, max) # Max expression per cluster
expression.threshold <- 0.15 * max.expression
expression.threshold[expression.threshold < 0.1] <- 0.1 # Minimum of 0.1, 20% of max expression.
exp.cluster.rel.threshold <- sweep(mean.exp.cluster, 1, expression.threshold, "/")
positive.clusters <- rowSums(exp.cluster.rel.threshold >= 1)
df$positive.clusters <- positive.clusters[rownames(df)]

# Define positive tissues based on positive clusters, rather than just expression level average within the tissue
tissue.of.clust <- unique(daniocell.annot[,c("clust", "tissue")])
tissue.of.clust <- tissue.of.clust[apply(tissue.of.clust, 1, function(x) !any(is.empty(x))),]
rownames(tissue.of.clust) <- tissue.of.clust$clust
df$positive.tissues <- unlist(apply(exp.cluster.rel.threshold[rownames(df),], 1, function(x) {
  length(setdiff(sort(unique(tissue.of.clust[names(which(x >= 1)), "tissue"])), c(NA, "")))
}))

# Define positive 'celltypes' based on positive clusters, rather than just expression level average within the tissue
celltype.of.clust <- unique(daniocell.annot[,c("clust", "clust.collapse")])
celltype.of.clust <- celltype.of.clust[apply(celltype.of.clust, 1, function(x) !any(is.empty(x))),]
rownames(celltype.of.clust) <- celltype.of.clust$clust
df$positive.ct <- unlist(apply(exp.cluster.rel.threshold[rownames(df),], 1, function(x) {
  length(setdiff(sort(unique(celltype.of.clust[names(which(x >= 1)), "clust.collapse"])), c(NA, "")))
}))

## CLASSIFY GENES BASED ON VARIATION -----------------

df$class.biorxiv <- "Unclassified"
df$class.devcell <- "Unclassified"

# Original biorxiv classification
df[df$CV.exp.cluster.log <= -0.85, "class.biorxiv"] <- "Ubiquitous"
df[df$CV.exp.cluster.log > -0.85 & df$CV.exp.cluster.log < 1.52, "class.biorxiv"] <- "Regulated"
df[df$CV.exp.cluster.log >= 1.52, "class.biorxiv"] <- "Specific"

# Dev Cell paper classifications
# "Housekeeping": Genes that are VERY consistent across time & tissue
genes.hk <- rownames(df)[df$CV.exp.cluster.log <= -0.85] # 863 (3.6%) genes - ubiquitous
df[genes.hk, "class.devcell"] <- "Housekeeping"
# "Ubiquitous": Genes that are in most tissues, but have a higher CV -- so exhibiting some variation
genes.ubi <- rownames(df)[df$positive.tissues >= 38 & df$CV.exp.cluster.log > -0.85]
df[genes.ubi, "class.devcell"] <- "Ubiquitous varying"
# "Cell-type specific" = high CV + expressed in a max of 3 cell types
genes.ctspec <- rownames(df)[df$CV.exp.cluster.log >= 1.5 & df$positive.ct <= 3] # 2706 genes
df[genes.ctspec, "class.devcell"] <- "Cell type-specific"
# "Cell-type restricted" = high CV + expressed in a max of 8 cell types
genes.ctrestr <- rownames(df)[df$CV.exp.cluster.log >= 1.5 & df$positive.ct > 3 & df$positive.ct <= 8] # 1502 genes
df[genes.ctrestr, "class.devcell"] <- "Cell type-restricted"
# "Tissue-specific" = moderately high CV, limited to 1 or 2 tissues, was not identified as cell-type specific.
genes.tissuespec <- setdiff(rownames(df)[df$CV.exp.cluster.log >= 1.25 & df$positive.tissues <= 2], genes.ctspec) # 748 genes
df[genes.tissuespec, "class.devcell"] <- "Tissue-specific"
# "Tissue-restricted" = moderately high CV, limited to 4 tissues, not identified as cell-type restricted.
genes.tissuerestr <- setdiff(rownames(df)[df$CV.exp.cluster.log >= 1.25 & df$positive.tissues > 2 & df$positive.tissues <= 4], c(genes.ctspec, genes.ctrestr)) # 600 genes
df[genes.tissuerestr, "class.devcell"] <- "Tissue restricted"


## PLOTS FOR THE PAPER -----------------------

# Category color palette
okabe.ito <- rgb(
  red = c(000, 000, 000, 086, 240, 230, 213, 204),
  green = c(000, 158, 114, 180, 228, 159, 094, 121),
  blue = c(000, 115, 178, 233, 066, 000, 000, 167),
  maxColorValue = 255
)
cell.class.colors.2 <- c(okabe.ito[c(4, 6, 7, 3, 2, 8)], "#CCCCCC")
names(cell.class.colors.2) <- c("Housekeeping", "Ubiquitous varying", "Cell type-specific", "Cell type-restricted", "Tissue-specific", "Tissue restricted", "Unclassified")

# Plot: CV of cluster mean expression vs. positive tissues for each gene, colored based on 'category'
pdf(width = 8, height = 7, file = "~/Desktop/Time_vs_CellTypes_ColoredPoints.pdf")
ggplot(data = df, aes(x = CV.exp.cluster.log, y = positive.ct, color = class.devcell)) + geom_point(alpha = 0.35, size = 0.75) +
  theme_bw() + scale_color_manual(values = cell.class.colors.2) +
  labs(x = "Variation between clusters (log CV of mean expression)", y = "Positive cell types", fill = "Number\nof genes") +
  geom_vline(xintercept = c(-0.85, 1.25, 1.5), lty = 'dashed', color = "#555555") + geom_hline(yintercept = c(3.5, 8.5), lty = 'dashed', color = "#555555")
dev.off()

# Plot: Variation between timepoints vs. Variation between cell types within a timepoint, colored based on 'category'.
pdf(width = 8, height = 7, file = "~/Desktop/Time_vs_CellTypes_PointsColored.pdf")
ggplot(data = df, aes(x = CV.exp.sg.log, y = CV.exp.ct.per.sg.log, color = class.devcell)) + geom_point(alpha = 0.5, size = 0.85) +
  theme_bw() + scale_color_manual(values = cell.class.colors.2) +
  labs(x = "Variation between timepoints (log CV of mean expression)", y = "Variation between cell types within a timepoint (max log CV of mean expression)")
dev.off()

# Plot: Distribution of cluster mean CV across all genes
ggplot(data = df, aes(x = CV.exp.cluster.log)) + geom_density(alpha = 0.2, fill = "#000000DD") + theme_bw() + labs(y = "Frequency", x = "CV of mean expression within clusters (log)") + geom_vline(xintercept = c(-0.85, 1.52), color = 'red', lty = 2)

# Plot: Distribution of cluster mean CV across each type of gene
plot.colors <- RColorBrewer::brewer.pal(n = 8, "Set1")[c(2,4,5,8)]
ggplot(data = df, aes(x = CV.exp.cluster.log, fill = type, color = type)) + geom_density(alpha = 0.35) + theme_bw() + labs(y = "Frequency", x = "CV of mean expression within clusters (log)", fill = "Gene type", color = "Gene type") + scale_fill_manual(breaks = c("Named", "cDNA", "Predicted", "XLOC"), values = plot.colors) + scale_color_manual(breaks = c("Named", "cDNA", "Predicted", "XLOC"), values = plot.colors) + geom_vline(xintercept = c(-0.85, 1.52), color = '#777777', lty = 2)

 