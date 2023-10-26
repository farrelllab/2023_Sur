##This code is designed to compare zebrafish intestinal cells to human intestinal cells by using correlation of a specific set of genes that are enriched in various intestinal cell types

library(Seurat)
library(pheatmap)

##Load human colonic dataset from Smillie et al., 2019
obj <- readRDS("~/Box/zfext/results/2022-11_human_intestine_data/human_colonic_dataset_Smillie_2019_seurat.rds")
obj1 <- obj

DimPlot(obj, label = T)

##Get cells that need to be excluded
cells.exclude <- WhichCells(obj1, idents = c("Goblet"), expression = AQP8 > 0.5)
cells.also.exclude <- WhichCells(obj1, idents = c("Goblet"), expression = GUCA2B > 0.5)
cells.also.exclude.1 <- WhichCells(obj1, idents = c("Goblet"), expression = GUCA2A > 0.5)
cells.also.exclude.2 <- WhichCells(obj1, idents = c("Goblet"), expression = SLC26A3 > 0.5)

cells.remove <- unlist(unique(list(cells.also.exclude, cells.also.exclude.1, cells.also.exclude.2)))

##Subset colon object without these doublet cells
cells.keep <- setdiff(WhichCells(obj), cells.remove)
obj <- subset(obj, cells = cells.keep)

##Scale expression per cluster
getZData <- function(matrix, genes=NULL) {
  # Z-scored data
  if (is.null(genes)) genes <- rownames(matrix)
  data.mean <- apply(matrix[genes,], 1, mean)
  data.sd <- apply(matrix[genes,], 1, sd)
  z.data <- as.matrix(sweep(sweep(matrix[genes,], 1, data.mean, "-"), 1, data.sd, "/"))
  z.data <- z.data[complete.cases(z.data),]
  return(z.data)
}

##Load the gene conversion table to convert between human and zebrafish gene names
gene.conversion.table <- read.table("~/Box/zfext/results/2022-11_human_intestine_data/HGNC_human_ZF_genenamesonly.tsv.gz")

##Load intestine seurat object
obj.intestine <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/endoderm/best4_markers/intestine_seurat.rds")
##Assign clusters
cluster.assignments <- Idents(obj.intestine)

# Find genes that we want to compare
gene.conversion.table.scaled <- gene.conversion.table
gene.conversion.table.scaled <- gene.conversion.table.scaled[gene.conversion.table.scaled$zebrafish_genes %in% intersect(rownames(obj.intestine@assays$RNA@data), gene.list.compare),]
gene.conversion.table.scaled <- gene.conversion.table.scaled[gene.conversion.table.scaled$human_genes != "",]
gene.conversion.table.scaled <- gene.conversion.table.scaled[gene.conversion.table.scaled$human_genes %in% rownames(obj@assays$RNA@data),]

##For zebrafish calculate average gene expression per cluster
clusters <- Idents(obj.intestine)
clusters.human.colon <- Idents(obj)

##Get genes to use for comparison
##Don't use the progenitor clusters
clusters.keep <- setdiff(unique(markers.int$cluster), c("int_progenitors-1", "int_progenitors-2"))
for(cluster in clusters.keep){
  ##Get genes that are quite specific to the cluster
  genes.use <- markers.int[which(markers.int$cluster == cluster & markers.int$avg_log2FC >= 0.25 & markers.int$pct.1 - markers.int$pct.2 >= 0.2), "gene"]
  genes.strong <- markers.int[which(markers.int$cluster == cluster & markers.int$avg_log2FC >= 2), "gene"]
  genes.keep <- unlist(unique(list(genes.use, genes.strong)))
  ##Now add that to the list 
  genes.keep <- list(genes.keep)
  names(genes.keep) <- cluster
  gene.list <- c(gene.list, genes.keep)
}

##Ok now curate the final marker list that will be used to calculate between the zebrafish intestine and human colon
gene.list.specific <- unlist(unique(gene.list))
gene.list.abs <- scan("~/Box/zfext/annotations_celltype_curated_newMama/endoderm/best4_markers/markers_zf_specific/markers_intestine_absorptive_only.txt", what = "character")
gene.list.sec <- scan("~/Box/zfext/annotations_celltype_curated_newMama/endoderm/best4_markers/markers_zf_specific/markers_intestine_secretory_only.txt", what = "character")

##Ok store FULL CURATED GENE LIST
gene.list.compare <- unlist(unique(list(gene.list.specific, gene.list.abs, gene.list.sec))) ##3317 genes 
##This above gene list will be used to compare zebrafish and human intestinal cell types

## Get zebrafish and human normalized data
zf.normalized.counts <- as.matrix(obj.intestine@assays$RNA@data[sort(unique(gene.conversion.table.scaled$zebrafish_genes)),])
human.colon.normalized.counts <- as.matrix(obj@assays$RNA@data[sort(unique(gene.conversion.table.scaled$human_genes)),])

## Sum expression when there are multiple orthologs
sum.of.logs <- function(x) {
  log(sum(exp(x) - 1)+1)
}

human.colon.norm.counts.sum.zf <- matrix(data = NA, nrow = nrow(zf.normalized.counts), ncol = ncol(human.colon.normalized.counts), dimnames = list(rownames(zf.normalized.counts), colnames(human.colon.normalized.counts)))

for (r in rownames(human.colon.norm.counts.sum.zf)) {
  hg <- gene.conversion.table.scaled[gene.conversion.table.scaled$zebrafish_genes == r, "human_genes"]
  if (length(hg) == 1) {
    human.colon.norm.counts.sum.zf[r,] <- human.colon.normalized.counts[hg,]
  } else {
    human.colon.norm.counts.sum.zf[r,] <- apply(human.colon.normalized.counts[hg,], 2, sum.of.logs)
  }
}

zf.norm.counts.sum.human <- matrix(data = NA, ncol = ncol(zf.normalized.counts), nrow = nrow(human.colon.normalized.counts), dimnames = list(rownames(human.colon.normalized.counts), colnames(zf.normalized.counts)))

for (r in rownames(zf.norm.counts.sum.human)) {
  zg <- gene.conversion.table.scaled[gene.conversion.table.scaled$human_genes == r, "zebrafish_genes"]
  if (length(zg) == 1) {
    zf.norm.counts.sum.human[r,] <- zf.normalized.counts[zg,]
  } else {
    zf.norm.counts.sum.human[r,] <- apply(zf.normalized.counts[zg,], 2, sum.of.logs)
  }
}

zf.norm.counts.sum.human.scaled <- getZData(zf.norm.counts.sum.human)
human.colon.norm.counts.sum.zf.scaled <- getZData(human.colon.norm.counts.sum.zf)
zf.scaled.counts <- getZData(zf.normalized.counts)
human.scaled.counts <- getZData(human.colon.normalized.counts)

# Average scaled data per cluster
agg.to.df <- function(x) {
  rownames(x) <- x$Group.1
  y <- as.data.frame(t(x[,2:ncol(x)]))
  return(y)
}

av.exp.zf.scaled <- aggregate(t(zf.scaled.counts), by = list(clusters[colnames(zf.scaled.counts)]), FUN = mean)
av.exp.human.scaled <- aggregate(t(human.scaled.counts), by = list(clusters.human[colnames(human.scaled.counts)]), FUN = mean)
av.exp.zf.to.human.scaled <- aggregate(t(zf.norm.counts.sum.human.scaled), by = list(clusters[colnames(zf.norm.counts.sum.human.scaled)]), FUN = mean)
av.exp.human.to.zf.scaled  <- aggregate(t(human.colon.norm.counts.sum.zf.scaled), by = list(clusters.human[colnames(human.colon.norm.counts.sum.zf.scaled)]), FUN = mean)

av.exp.zf.scaled <- agg.to.df(av.exp.zf.scaled)
av.exp.human.scaled <- agg.to.df(av.exp.human.scaled)
av.exp.zf.to.human.scaled <- agg.to.df(av.exp.zf.to.human.scaled)
av.exp.human.to.zf.scaled <- agg.to.df(av.exp.human.to.zf.scaled)

# Calculate correlations in each direction and also try mean of the two directions
genes.cor.zf <- intersect(rownames(av.exp.zf.scaled), rownames(av.exp.human.to.zf.scaled))
genes.cor.human <- intersect(rownames(av.exp.human.scaled), rownames(av.exp.zf.to.human.scaled))

cor.zf <- as.data.frame(cor(x = av.exp.zf.scaled[genes.cor.zf,], y = av.exp.human.to.zf.scaled[genes.cor.zf,]))
cor.human <- as.data.frame(cor(x = av.exp.zf.to.human.scaled[genes.cor.human,], y = av.exp.human.scaled[genes.cor.human,]))
cor.mean <- matrix(NA, nrow = nrow(cor.zf), ncol = ncol(cor.zf), dimnames = list(rownames(cor.zf), colnames(cor.zf)))
for (i in 1:nrow(cor.zf)) cor.mean[i,] <- colMeans(rbind(cor.zf[i,], cor.human[i,]))

# Heatmaps
heatmap.colors <- colorRampPalette(c("blue", "white", "red"))(99)
heatmap.colors <- c(heatmap.colors[1:49], rep("#FFFFFF", 20), heatmap.colors[51:99])
pheatmap(cor.mean, color = heatmap.colors)

##calculate correlation on the zebrafish mtrix
cor.zf.exp <- as.data.frame(cor(av.exp.zf.scaled))
pheatmap(cor.zf.exp, color = heatmap.colors)

##Create a duplicate of the rows for the zebrafish best4+ cell expression and converted human colon best4 cell expression
av.exp.zf.scaled$genes <- rownames(av.exp.zf.scaled)
av.exp.human.to.zf.scaled$genes <- rownames(av.exp.human.to.zf.scaled)

##Now subset the dataframe for the scaled expression for best4 cells for zebrafish and colon
exp.scaled.best4.colon <- av.exp.human.to.zf.scaled[, c("Best4+ Enterocytes", "genes")]
exp.scaled.best4.zf <- av.exp.zf.scaled[, c("best4+_enterocytes", "genes")]

##Do the same as above just by using the bioMart list
biomart_conversion_table <- read.delim("~/Box/zfext/results/2022-11_human_intestine_data/ZF_human_homologs/zebrafish-to-human_gene_conversion_biomaRt.tsv", sep = "\t")

##In the biomart_list find how many are 1:1 relationships compared to the full gene conversion table
zebrafish_match <- character()
human_match <- character()

df <- gene.conversion.table.scaled
colnames(df) <- c("zebrafish_genes", "human_genes")
# Iterate through each row
for (i in 1:nrow(df)) {
  # Check matches in column 1 (zebrafish_genes)
  matches_col1 <- sum(df$zebrafish_genes == df$zebrafish_genes[i])
  
  # Check matches in column 2 (human_genes)
  matches_col2 <- sum(df$human_genes == df$human_genes[i])
  
  # Check if it's a 1:1 match
  if (matches_col1 == 1 && matches_col2 == 1) {
    # Save the gene names in separate lists
    zebrafish_match <- c(zebrafish_match, df$zebrafish_genes[i])
    human_match <- c(human_match, df$human_genes[i])
  }
}

zebrafish_match <- unique(zebrafish_match)
human_match <- unique(human_match)

# Find genes that we want to compare - limit to just zebrafish and human gene names
gene.conversion.table.scaled <- biomart_conversion_table[, c("UniProtKB.Gene.Name.symbol", "UniProtKB.Gene.Name.symbol.1")]
colnames(gene.conversion.table.scaled) <- c("zebrafish_genes", "human_genes")

gene.conversion.table.scaled <- gene.conversion.table.scaled[gene.conversion.table.scaled$zebrafish_genes %in% intersect(rownames(obj.intestine@assays$RNA@data), gene.list.compare),]
gene.conversion.table.scaled <- gene.conversion.table.scaled[gene.conversion.table.scaled$human_genes != "",]
gene.conversion.table.scaled <- gene.conversion.table.scaled[gene.conversion.table.scaled$human_genes %in% rownames(obj@assays$RNA@data),]

## Get zebrafish and human normalized data
zf.normalized.counts <- as.matrix(obj.intestine@assays$RNA@data[sort(unique(gene.conversion.table.scaled$zebrafish_genes)),])
human.colon.normalized.counts <- as.matrix(obj@assays$RNA@data[sort(unique(gene.conversion.table.scaled$human_genes)),])

## Sum expression when there are multiple orthologs
sum.of.logs <- function(x) {
  log(sum(exp(x) - 1)+1)
}

human.colon.norm.counts.sum.zf <- matrix(data = NA, nrow = nrow(zf.normalized.counts), ncol = ncol(human.colon.normalized.counts), dimnames = list(rownames(zf.normalized.counts), colnames(human.colon.normalized.counts)))

for (r in rownames(human.colon.norm.counts.sum.zf)) {
  hg <- gene.conversion.table.scaled[gene.conversion.table.scaled$zebrafish_genes == r, "human_genes"]
  if (length(hg) == 1) {
    human.colon.norm.counts.sum.zf[r,] <- human.colon.normalized.counts[hg,]
  } else {
    human.colon.norm.counts.sum.zf[r,] <- apply(human.colon.normalized.counts[hg,], 2, sum.of.logs)
  }
}

zf.norm.counts.sum.human <- matrix(data = NA, ncol = ncol(zf.normalized.counts), nrow = nrow(human.colon.normalized.counts), dimnames = list(rownames(human.colon.normalized.counts), colnames(zf.normalized.counts)))

for (r in rownames(zf.norm.counts.sum.human)) {
  zg <- gene.conversion.table.scaled[gene.conversion.table.scaled$human_genes == r, "zebrafish_genes"]
  if (length(zg) == 1) {
    zf.norm.counts.sum.human[r,] <- zf.normalized.counts[zg,]
  } else {
    zf.norm.counts.sum.human[r,] <- apply(zf.normalized.counts[zg,], 2, sum.of.logs)
  }
}

zf.norm.counts.sum.human.scaled <- getZData(zf.norm.counts.sum.human)
human.colon.norm.counts.sum.zf.scaled <- getZData(human.colon.norm.counts.sum.zf)
zf.scaled.counts <- getZData(zf.normalized.counts)
human.scaled.counts <- getZData(human.colon.normalized.counts)

# Average scaled data per cluster
agg.to.df <- function(x) {
  rownames(x) <- x$Group.1
  y <- as.data.frame(t(x[,2:ncol(x)]))
  return(y)
}

av.exp.zf.scaled <- aggregate(t(zf.scaled.counts), by = list(clusters[colnames(zf.scaled.counts)]), FUN = mean)
av.exp.human.scaled <- aggregate(t(human.scaled.counts), by = list(clusters.human[colnames(human.scaled.counts)]), FUN = mean)
av.exp.zf.to.human.scaled <- aggregate(t(zf.norm.counts.sum.human.scaled), by = list(clusters[colnames(zf.norm.counts.sum.human.scaled)]), FUN = mean)
av.exp.human.to.zf.scaled  <- aggregate(t(human.colon.norm.counts.sum.zf.scaled), by = list(clusters.human[colnames(human.colon.norm.counts.sum.zf.scaled)]), FUN = mean)

av.exp.zf.scaled <- agg.to.df(av.exp.zf.scaled)
av.exp.human.scaled <- agg.to.df(av.exp.human.scaled)
av.exp.zf.to.human.scaled <- agg.to.df(av.exp.zf.to.human.scaled)
av.exp.human.to.zf.scaled <- agg.to.df(av.exp.human.to.zf.scaled)

# Calculate correlations in each direction and also try mean of the two directions
genes.cor.zf <- intersect(rownames(av.exp.zf.scaled), rownames(av.exp.human.to.zf.scaled))
genes.cor.human <- intersect(rownames(av.exp.human.scaled), rownames(av.exp.zf.to.human.scaled))

cor.zf <- as.data.frame(cor(x = av.exp.zf.scaled[genes.cor.zf,], y = av.exp.human.to.zf.scaled[genes.cor.zf,]))
cor.human <- as.data.frame(cor(x = av.exp.zf.to.human.scaled[genes.cor.human,], y = av.exp.human.scaled[genes.cor.human,]))
cor.mean <- matrix(NA, nrow = nrow(cor.zf), ncol = ncol(cor.zf), dimnames = list(rownames(cor.zf), colnames(cor.zf)))
for (i in 1:nrow(cor.zf)) cor.mean[i,] <- colMeans(rbind(cor.zf[i,], cor.human[i,]))

# Heatmaps
heatmap.colors <- colorRampPalette(c("blue", "white", "red"))(99)
# Update the colors to make the range -0.075 to 0.075 appear as white
heatmap.colors <- c(heatmap.colors[1:48], rep("#FFFFFF", 18), heatmap.colors[50:99])
# Update the colors to make -0.1 to 0.1 represented as white
pheatmap(cor.mean, color = heatmap.colors)

##Reorder the columns of the heatmap 
column_order <- c("Stem", "TA 1", "TA 2", "Cycling TA", "Secretory TA", "Immature Enterocytes 1", "Immature Enterocytes 2", "Enterocyte Progenitors", "Enterocytes", "Immature Goblet", "Goblet", "M cells", "Best4+ Enterocytes", "Enteroendocrine", "Tuft")
column_order <- c(1, 2, 3, 4, 11, 7, 5, 6, 8, 12, 13, 9, 10, 15, 14)

row_order <- c("int_progenitors-1", "int_progenitors-2", "enterocyte-1", "enterocyte-2", "enterocyte-3", "goblet-cells", "posterior-LREs", "best4+_enterocytes", "EECs", "tuft-like")
row_order <- c(1, 6, 2, 3, 4, 7, 5, 8, 9, 10)

##Specify row and column order
cor_mean_reordered <- cor.mean[row_order, column_order]

##Draw reordered heatmap
pheatmap(cor_mean_reordered, cluster_rows = FALSE, cluster_cols = FALSE, color = heatmap.colors,
         breaks = c(seq(min(cor_mean_reordered), -0.075, length.out = 58), 0, seq(0.075, max(cor_mean_reordered), length.out = 58)), legend_breaks = seq(-0.4, 0.4, by = 0.2), legend_labels = c("-0.4", "-0.2", "0", "0.2", "0.4"))

pheatmap(cor_mean_reordered, cluster_rows = FALSE, cluster_cols = FALSE, color = heatmap.colors,
         breaks = c(seq(-0.4, -0.075, length.out = 58), 0, seq(0.075, 0.4, length.out = 58)), legend_breaks = seq(-0.4, 0.4, by = 0.2), legend_labels = c("-0.4", "-0.2", "0", "0.2", "0.4"))

##calculate correlation on the zebrafish matrix using the full list of specific genes
##For zebrafish
#Step 1: Get normalized counts
norm.counts.zf <- obj.intestine@assays$RNA@data[gene.list.compare, ]

##Step 2: Scale the counts
norm.counts.zf.scaled <- getZData(norm.counts.zf)

##Step 3: Aggregate the normalized counts based on clusters
norm.counts.zf.scaled.mean <- aggregate(t(norm.counts.zf.scaled), by = list(clusters[colnames(norm.counts.zf.scaled)]), FUN = mean)

##Convert to dataframe
norm.counts.zf.scaled.mean <- agg.to.df(norm.counts.zf.scaled.mean)

##Step 5: Calculate correlation
cor.zf.exp <- as.data.frame(cor(norm.counts.zf.scaled.mean))
cor.zf.exp_reordered <- cor.zf.exp[row_order, row_order]

##Heatmap
pheatmap(cor.zf.exp_reordered, cluster_rows = FALSE, cluster_cols = FALSE, color = heatmap.colors,
         breaks = c(seq(min(cor_mean_reordered), -0.1, length.out = 58), 0, seq(0.1, max(cor_mean_reordered), length.out = 58)), legend_breaks = seq(-1, 1, by = 0.5), legend_labels = c("-1", "-0.5", "0", "0.5", "1"))

heatmap.colors <- colorRampPalette(c("blue", "white", "red"))(99)
heatmap.colors <- c(heatmap.colors[1:48], rep("#FFFFFF", 2), heatmap.colors[50:99])
pheatmap(cor.zf.exp_reordered, cluster_rows = FALSE, cluster_cols = FALSE, color = heatmap.colors)

pheatmap(cor.zf.exp_reordered, cluster_rows = FALSE, cluster_cols = FALSE, color = heatmap.colors, 
         breaks = c(seq(-1, -0.4, length.out = 49), 0, seq(0.4, 1, length.out = 49)))

##For human colon data, generate correlation data
human_genes <- gene.conversion.table.scaled[which(gene.conversion.table.scaled$zebrafish_genes %in% intersect(rownames(obj.intestine@assays$RNA@data), gene.list.compare)), "human_genes"]
human_genes <- rownames(human.colon.normalized.counts)

##Step 1: Get normalized counts
norm.counts.human <- obj@assays$RNA@data[human_genes, ]

##Step 2: Scale the counts
norm.counts.human.scaled <- getZData(norm.counts.human)

##Step 3: Aggregate the normalized counts based on clusters
norm.counts.human.scaled.mean <- aggregate(t(norm.counts.human.scaled), by = list(clusters.human[colnames(norm.counts.human.scaled)]), FUN = mean)

##Step 4: Convert to dataframe
norm.counts.human.scaled.mean <- agg.to.df(norm.counts.human.scaled.mean)

##Step 5: Calculate correlation
cor.human.exp <- as.data.frame(cor(norm.counts.human.scaled.mean))
cor.human.exp.reordered <- cor.human.exp[column_order, column_order]

##Heatmap
pheatmap(cor.human.exp, color = heatmap.colors)
pheatmap(cor.human.exp.reordered, cluster_rows = FALSE, cluster_cols = FALSE, color = heatmap.colors,
         breaks = c(seq(min(cor.human.exp.reordered), -0.01, length.out = 55), 0, seq(0.01, max(cor.human.exp.reordered), length.out = 50)))


##Now to have more confidence, now calculate correlation between the avg. expression of human genes (converted to zebrafish) and human genes (converted to zebrafish)
##The matrix was already calculated previously
head(av.exp.human.to.zf.scaled)  ##zebrafish genes as rownames and human clusters as colnames - 1768 genes
head(av.exp.zf.to.human.scaled) ##Human genes as rownames and zebrafish clusters as colnames - 2525 genes

##Calculate correlation
cor.human.to.zf.converted <- as.data.frame(cor(av.exp.human.to.zf.scaled))
cor.zf.to.human.converted <- as.data.frame(cor(av.exp.zf.to.human.scaled))

##Plot heatmap
pheatmap(cor.human.to.zf.converted, color = heatmap.colors)
pheatmap(cor.zf.to.human.converted, color = heatmap.colors)

##Plot scatter plot for zebrafish intestine best4+ cells versus human colon best4+ cells
##Subset the human expression matrix for just best4+ cells




############################ ZEBRAFISH INTESTINE VS SMALL INTESTINE ############################################
##Now read in the small intestine object
obj.si <- readRDS("~/Box/zfext/results/2022-11_human_intestine_data/Burclaff_2022_human_small_intestine_downloaded_seurat.rds")

##Plot UMAP - need to assign clusters
DimPlot(obj.si, label = T)
FeaturePlot(obj.si, c("MUC2", "AGR2", "BEST4", "OTOP2"))

##Get small intestine clusters 
cells.jejunum <- rownames(obj.si@meta.data)[which(obj.si@meta.data$tissue == "jejunum")]
cells.ileum <- rownames(obj.si@meta.data)[which(obj.si@meta.data$tissue == "ileum")]
cells.duodenum <- rownames(obj.si@meta.data)[which(obj.si@meta.data$tissue == "duodenum")]

##Get all cells for the small intestine
cells.small.intestine <- unlist(unique(list(cells.jejunum, cells.ileum, cells.duodenum)))

##Make a small intestine object - use only small intestine specific cells
obj.si <- subset(obj.si, cells = cells.small.intestine)

Idents(obj.si) <- obj.si@meta.data$type
clusters.human.SI <- Idents(obj.si)

##Gene conversion table - from Ensembl (Biomart) that Gennady generated
gene.conversion.table.full <- read.delim("~/Box/zfext/results/2022-11_human_intestine_data/ZF_human_homologs/zebrafish-to-human_gene_conversion_biomaRt.tsv", sep = "\t")

##Directly use the table with the ensembl ids and the human gene names
gene.ensembl.ids.human.table <- read.table("~/Box/zfext/results/2022-11_human_intestine_data/ZF_human_homologs/gene_orthologs/hs_dr_orthologs_ENSEMBL-SYMBOL.tsv")

##Get all the ensembl ids
ensembl_ids <- unique(rownames(obj@assays$RNA@data))
#Create named vector of gene names
gene_names <- setNames(gene.conversion.table.sub$ensembl_ids, gene.conversion.table.sub$human_genes)
gene_names <- unique(gene_names)

##Load the HGNC table that Jeff generated
gene.conversion.table <- read.delim("~/Box/zfext/results/2022-11_human_intestine_data/HGNC_human_ZF_genenamesonly.tsv.gz", sep = "\t")

# Find genes that we want to compare
gene.conversion.table.sub <- gene.ensembl.ids.human.table[, c(8, 9)]
colnames(gene.conversion.table.sub) <- c("ensembl_ids", "human_genes")
gene.conversion.table.sub <- unique(gene.conversion.table.sub)

# Step 1: Create named vector of gene names
rna_data <- as.matrix(obj@assays$RNA@data)
rna_data <- rna_data[which(rownames(rna_data) %in% gene.conversion.table.sub$ensembl_ids), ]
# Get the matching gene names from the conversion table
matching_genes <- gene.conversion.table.sub[which(rownames(rna_data) %in% gene.conversion.table.sub$ensembl_ids), "human_genes"]

rownames(rna_data) <- matching_genes
rna_data <- rna_data[rownames(rna_data) != "", ]
rna_data <- unique(rna_data)

# Convert the modified matrix back to dgCMatrix
rna_data.2 <- as(rna_data, "dgCMatrix")

obj@assays$RNA@data <- rna_data

##Scale expression per cluster
getZData <- function(matrix, genes=NULL) {
  # Z-scored data
  if (is.null(genes)) genes <- rownames(matrix)
  data.mean <- apply(matrix[genes,], 1, mean)
  data.sd <- apply(matrix[genes,], 1, sd)
  z.data <- as.matrix(sweep(sweep(matrix[genes,], 1, data.mean, "-"), 1, data.sd, "/"))
  z.data <- z.data[complete.cases(z.data),]
  return(z.data)
}

##Load gene conversion table
gene.conversion.table <- read.table("~/Box/zfext/results/2022-11_human_intestine_data/HGNC_human_ZF_genenamesonly.tsv.gz")
colnames(gene.conversion.table) <- c("zebrafish_genes", "human_genes")

# Find genes that we want to compare
gene.conversion.table.scaled <- gene.conversion.table
gene.conversion.table.scaled <- gene.conversion.table.scaled[gene.conversion.table.scaled$zebrafish_genes %in% intersect(rownames(obj.intestine@assays$RNA@data), gene.list.compare),]
gene.conversion.table.scaled <- gene.conversion.table.scaled[gene.conversion.table.scaled$human_genes != "",]
gene.conversion.table.scaled <- gene.conversion.table.scaled[gene.conversion.table.scaled$human_genes %in% rownames(obj@assays$RNA@data),]

##For zebrafish calculate average gene expression per cluster
clusters <- Idents(obj.intestine)
clusters.human <- Idents(obj)


## Get zebrafish and human normalized data
zf.normalized.counts <- as.matrix(obj.intestine@assays$RNA@data[sort(unique(gene.conversion.table.scaled$zebrafish_genes)),])
human.SI.normalized.counts <- as.matrix(obj@assays$RNA@data[sort(unique(gene.conversion.table.scaled$human_genes)),])

## Sum expression when there are multiple orthologs
sum.of.logs <- function(x) {
  log(sum(exp(x) - 1)+1)
}

human.SI.norm.counts.sum.zf <- matrix(data = NA, nrow = nrow(zf.normalized.counts), ncol = ncol(human.SI.normalized.counts), dimnames = list(rownames(zf.normalized.counts), colnames(human.SI.normalized.counts)))

for (r in rownames(human.SI.norm.counts.sum.zf)) {
  hg <- gene.conversion.table.scaled[gene.conversion.table.scaled$zebrafish_genes == r, "human_genes"]
  if (length(hg) == 1) {
    human.SI.norm.counts.sum.zf[r,] <- human.SI.normalized.counts[hg,]
  } else {
    human.SI.norm.counts.sum.zf[r,] <- apply(human.SI.normalized.counts[hg,], 2, sum.of.logs)
  }
}

zf.norm.counts.sum.human <- matrix(data = NA, ncol = ncol(zf.normalized.counts), nrow = nrow(human.SI.normalized.counts), dimnames = list(rownames(human.SI.normalized.counts), colnames(zf.normalized.counts)))

for (r in rownames(zf.norm.counts.sum.human)) {
  zg <- gene.conversion.table.scaled[gene.conversion.table.scaled$human_genes == r, "zebrafish_genes"]
  if (length(zg) == 1) {
    zf.norm.counts.sum.human[r,] <- zf.normalized.counts[zg,]
  } else {
    zf.norm.counts.sum.human[r,] <- apply(zf.normalized.counts[zg,], 2, sum.of.logs)
  }
}

zf.norm.counts.sum.human.scaled <- getZData(zf.norm.counts.sum.human)
human.SI.norm.counts.sum.zf.scaled <- getZData(human.SI.norm.counts.sum.zf)
zf.scaled.counts <- getZData(zf.normalized.counts)
human.scaled.counts <- getZData(human.SI.normalized.counts)

# Average scaled data per cluster
agg.to.df <- function(x) {
  rownames(x) <- x$Group.1
  y <- as.data.frame(t(x[,2:ncol(x)]))
  return(y)
}

av.exp.zf.scaled <- aggregate(t(zf.scaled.counts), by = list(clusters[colnames(zf.scaled.counts)]), FUN = mean)
av.exp.human.scaled <- aggregate(t(human.scaled.counts), by = list(clusters.human[colnames(human.scaled.counts)]), FUN = mean)
av.exp.zf.to.human.scaled <- aggregate(t(zf.norm.counts.sum.human.scaled), by = list(clusters[colnames(zf.norm.counts.sum.human.scaled)]), FUN = mean)
av.exp.human.to.zf.scaled  <- aggregate(t(human.SI.norm.counts.sum.zf.scaled), by = list(clusters.human[colnames(human.SI.norm.counts.sum.zf.scaled)]), FUN = mean)

av.exp.zf.scaled <- agg.to.df(av.exp.zf.scaled)
av.exp.human.scaled <- agg.to.df(av.exp.human.scaled)
av.exp.zf.to.human.scaled <- agg.to.df(av.exp.zf.to.human.scaled)
av.exp.human.to.zf.scaled <- agg.to.df(av.exp.human.to.zf.scaled)

##Subset the average expression table for just the best4+ cells
##Duplicate the rownames and save them as a column in the dataframe
av.exp.human.to.zf.scaled$genes <- rownames(av.exp.human.to.zf.scaled)
av.exp.zf.scaled$genes <- rownames(av.exp.zf.scaled)

# Calculate correlations in each direction and also try mean of the two directions
genes.cor.zf <- intersect(rownames(av.exp.zf.scaled), rownames(av.exp.human.to.zf.scaled))
genes.cor.human <- intersect(rownames(av.exp.human.scaled), rownames(av.exp.zf.to.human.scaled))

cor.zf <- as.data.frame(cor(x = av.exp.zf.scaled[genes.cor.zf,], y = av.exp.human.to.zf.scaled[genes.cor.zf,]))
cor.human <- as.data.frame(cor(x = av.exp.zf.to.human.scaled[genes.cor.human,], y = av.exp.human.scaled[genes.cor.human,]))
cor.mean <- matrix(NA, nrow = nrow(cor.zf), ncol = ncol(cor.zf), dimnames = list(rownames(cor.zf), colnames(cor.zf)))
for (i in 1:nrow(cor.zf)) cor.mean[i,] <- colMeans(rbind(cor.zf[i,], cor.human[i,]))

# Heatmaps
heatmap.colors <- colorRampPalette(c("blue", "white", "red"))(99)
heatmap.colors <- c(heatmap.colors[1:49], rep("#FFFFFF", 20), heatmap.colors[51:99])
pheatmap(cor.human, color = heatmap.colors)

##calculate correlation on the zebrafish mtrix
cor.zf.exp <- as.data.frame(cor(av.exp.zf.scaled))
pheatmap(cor.zf.exp, color = heatmap.colors)


########################## RE-PLOTTING THE SCATTER PLOTS WITH THIS NEW APPROACH ########################
##Normalized counts for the human colon, human SI and zebrafish intestine 
av.exp.zf.intestine <- aggregate(t(zf.normalized.counts), by = list(clusters[colnames(zf.normalized.counts)]), FUN = mean)

##Get normalized expression for human colon genes and convert to zebrafish ids
av.exp.human.colon <- aggregate(t(human.colon.norm.counts.sum.zf), by = list(clusters.human.colon[colnames(human.colon.norm.counts.sum.zf)]), FUN = mean)
av.exp.zf.to.human.colon <- aggregate(t(zf.normalized.counts), by = list(clusters[colnames(zf.normalized.counts)]), FUN = mean)
av.exp.human.to.zf.colon <- aggregate(t(human.colon.norm.counts.sum.zf), by = list(clusters.human.colon[colnames(human.colon.norm.counts.sum.zf)]), FUN = mean)

##Get normalized expression for human SI genes and convert to zebrafish ids
av.exp.human.SI <- aggregate(t(human.SI.norm.counts.sum.zf), by = list(clusters.human.SI[colnames(human.SI.norm.counts.sum.zf)]), FUN = mean)
av.exp.human.to.zf.SI <- aggregate(t(human.SI.norm.counts.sum.zf), by = list(clusters.human.SI[colnames(human.SI.norm.counts.sum.zf)]), FUN = mean)

##Convert the aggregated expression matrices to dataframe
av.exp.zf.intestine <- agg.to.df(av.exp.zf.intestine)
av.exp.human.colon <- agg.to.df(av.exp.human.colon)
av.exp.human.SI <- agg.to.df(av.exp.human.SI)
av.exp.human.to.zf.colon <- agg.to.df(av.exp.human.to.zf.colon)
av.exp.human.to.zf.SI <- agg.to.df(av.exp.human.to.zf.SI)
av.exp.human.to.zf.SI.2 <- na.omit(av.exp.human.to.zf.SI)

##Duplicate the rownames to create a separate column for the gene names in each of the above datasets
av.exp.zf.intestine$genes <- rownames(av.exp.zf.intestine)
av.exp.human.to.zf.colon$genes <- rownames(av.exp.human.to.zf.colon)
av.exp.human.to.zf.SI.2$genes <- rownames(av.exp.human.to.zf.SI.2)

##Subset the above dataframe with the column with gene names and expression of best4+ cells
exp.scaled.best4.si <- av.exp.human.to.zf.SI.2[, c("BEST4+", "genes")]
exp.scaled.best4.zf <- av.exp.zf.intestine[, c("best4+_enterocytes", "genes")]
exp.scaled.best4.colon <- av.exp.human.to.zf.colon[, c("Best4+ Enterocytes", "genes")]

##Make the column names similar between the colon and SI dataframes for rbind
colnames(exp.scaled.best4.colon) <- c("BEST4+", "genes")

##Now find genes that are expressed in SI best4 cells but not colon best4 cells
df.human <- rbind(exp.scaled.best4.colon, exp.scaled.best4.si)

##Now find genes that are common between the human dataframe (SI+colon) and the zebrafish dataframe
genes.common <- intersect(rownames(df.human), rownames(exp.scaled.best4.zf))

##Now make a common dataframe with comparative scaled expression of human and ZF best4+ cells
dx <- cbind(df.human[genes.common, ], exp.scaled.best4.zf[genes.common, ])
dx <- dx[, -c(2, 4)]
colnames(dx) <- c("log_exp_human", "log_exp_zf")

##Now plot scatter plot
plot(dx$log_exp_human, dx$log_exp_zf, main = "Expression in humans versus ZF ",
     xlab="human_BEST4_expression", 
     ylab="ZF_BEST4_expression", 
     cex=0.95, cex.lab= 0.5, title(cex.main=2, font.main=7))
abline(-1, 1, col='purple', lty="dashed")
abline(1, 1, col="purple", lty="dashed")
above <- dx$log_exp_zf - dx$log_exp_human > 1
points(dx$log_exp_human[above], dx$log_exp_zf[above], col="red", cex=0.5)
below <- dx$log_exp_zf - dx$log_exp_human < -1
points(dx$log_exp_human[below], dx$log_exp_zf[below], col="blue", cex=0.5)


ggplot(dx, aes(x = log_exp_human, y = log_exp_zf, label = genes)) +
  geom_point(size=2, shape=16) +
  geom_text()





