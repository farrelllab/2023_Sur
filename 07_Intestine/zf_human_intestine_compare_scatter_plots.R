##This code is written to plot the scatter plots shown in Figures 7B and Figures S9E, F

##Load libraries
library(Seurat)
library(dplyr)

##Comparing our zebrafish best4+ enterocyte expression with human colonic or large intestinal best4 enterocytes
##We used the paper - Smilie et al., 2019 Intra- and Inter-cellular Rewiring of the Human Colon during Ulcerative Colitis
##(https://www.cell.com/cell/fulltext/S0092-8674(19)30732-9?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867419307329%3Fshowall%3Dtrue)

##Load counts from the publicly available colon scRNAseq data
message(paste0(Sys.time(), ": load counts matrix"))
counts <- ReadMtx("~/Desktop/human_intestine_data/gene_sorted-Epi.matrix.mtx", cells = "Epi.barcodes2.tsv", features = "Epi.genes.tsv", feature.column = 1)
counts.keep <- counts[, cells.keep]

##Load metadata
message(paste0(Sys.time(), ": load meta.data"))
meta <- read.table("~/Desktop/human_intestine_data/all_meta.txt", header = TRUE, sep = "\t")

##Subset meta.data to only the epithelial cells
colnames(meta) <- meta[1:2, ]
meta <- meta[-1, ]
colnames(meta) <- c("name", "cluster", "nGene", "nUMI", "Subject", "Health", "Location", "group", "ID")
meta <- meta[, -9]
meta.keep <- meta[!duplicated(meta$name), ]
rownames(meta.keep) <-  meta.keep$name
meta.keep <- meta.keep[ , -1]
colnames(meta.keep) <- c("cluster", "identity.sub", "nGene", "nUMI", "Subject", "Health", "Location")

##Save the curated metadata
saveRDS(meta.keep, "~/Desktop/human_intestine_data/all_meta_reordered.rds")
meta.keep.only <- readRDS("~/Desktop/human_intestine_data/all_meta_reordered.rds")

##Create seurat object with the counts matrix
message(paste0(Sys.time(), ": create seurat object"))
obj <- CreateSeuratObject(counts = counts.keep, min.cells = 3, min.features = 1)

##Only keep cells present in the counts matrix
cells.keep <- intersect(meta$name, colnames(obj@assays$RNA@counts))
obj <- obj@assays$RNA@counts[, cells.keep]

##Subset metadata to only these cells
meta.keep.only <- meta.keep[which(rownames(meta.keep) %in% cells.keep), ]

##Add metadata columns
obj <- AddMetaData(obj, metadata = meta.keep.only)

##Save the seurat object containing  the human colon information
saveRDS(obj, "~/Desktop/human_intestine_data/human_intestine_seurat.rds")
obj <- readRDS("~/Desktop/human_intestine_data/human_intestine_seurat.rds")

##Do the usual normalization, PCA, neighbors, UMAP/TSNE
obj <- NormalizeData(object = obj)
obj <- FindVariableFeatures(object = obj)
obj <- ScaleData(object = obj)
obj <- RunPCA(obj, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)

# Dimensionality Reduction: UMAP
message(paste0(Sys.time(), ": UMAP"))
#nn.use <- floor(sqrt(ncol(obj@assays$RNA@counts))/10)*10
nn.use <- 30
obj <- RunUMAP(obj, dims=1:30, n.neighbors=nn.use)

# Clustering
message(paste0(Sys.time(), ": Clustering"))
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, nn.eps = 0, k.param = 50)

obj <- FindNeighbors(object = obj)
obj <- FindClusters(object = obj, resolution = c(2,3))
obj <- RunTSNE(object = obj, dims = 1:30, nn.use = 50)

DimPlot(obj, group.by = "RNA_snn_res.2")

DimPlot(object = obj, label = T)
DimPlot(obj, cells.highlight = WhichCells(obj, idents = c("4")))

##Load TSNE coordinates
tsne.coord <- read.table("~/Desktop/human_intestine_data/Epi.tsne.txt", header = FALSE, sep = "\t")
rownames(tsne.coord) <- tsne.coord$V1
tsne.coord <- tsne.coord[, -1]
colnames(tsne.coord) <- c("tSNE.1", "tSNE_2")

obj@reductions$umap <- NA
obj@reductions$umap@cell.embeddings <- as.matrix(tsne.coord)

##Ok subset tSNE dataframe to just the cells to keep
cells.keep <- intersect(rownames(tsne.coord), rownames(obj@meta.data))
tsne.coord <- tsne.coord[which(rownames(tsne.coord) %in% cells.keep), ]

##Saving
saveRDS(obj, "~/Desktop/human_intestine_data/human_intestine_seurat.rds")
obj <- readRDS("~/Desktop/human_intestine_data/human_intestine_seurat.rds")

##Plot some genes
FeaturePlot(obj, c("BEST4", "OTOP2", "CFTR", "ATOH1"))

##Cluster using the "cluster" column in metadata
DimPlot(obj, group.by = "cluster")
Idents(obj) <- obj@meta.data$cluster

cells.best4 <- WhichCells(obj, idents = c("Best4+"))

##Calculate markers for Best4+ cells that have an average log FC of more than 2 or 3
markers.human <- FindMarkers(obj, ident.1 = "Best4+", ident.2 = setdiff(levels(Idents(obj)), "Best4+"), logfc.threshold = 0.25, min.pct = 0.25)
saveRDS(markers.human, "~/Desktop/human_intestine_data/human_intestine_markers.rds")

##Calculate average expression of all cells within a cluster
cluster.averages <- AverageExpression(obj)
head(cluster.averages[["RNA"]])

orig.levels <- levels(obj)
Idents(obj) <- gsub(pattern = " ", replacement = "_", x = Idents(obj))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(obj) <- orig.levels
cluster.averages <- AverageExpression(obj, return.seurat = TRUE)
cluster.averages

##Plot a scatter plot between the various cell types
CellScatter(cluster.averages, cell1 = "Best4+", cell2 = "Goblet")

##Ok, now put a threshold for Best4+ enterocytes for genes that are expressed above a certain threshold
message(paste0(Sys.time(), ": Convert cluster averages to dataframe"))
df <- as.data.frame(cluster.averages)

## Get all genes expressed in the human Best4+ enterocytes that have a log-fold expression more than 3. 
genes.best4 <- rownames(df)[which(df$RNA.Best4. > 3)]

##Filter the list of genes to remove all ribosomal and mitochondrial genes
message(paste0(Sys.time(), ": Convert cluster averages to dataframe"))
filter.genes <- function(genes) {
  mt.genes <- grep("^mt-", ignore.case = T, genes, value = T)
  many.genes <- grep("\\(1 of many\\)", ignore.case = T, genes, value = T)
  cox.genes <- grep("^cox", ignore.case = T, genes, value = T)
  rps.genes <- grep("^rps", ignore.case = T, genes, value = T)
  rpl.genes <- grep("^rpl", ignore.case = T, genes, value = T)
  return(setdiff(genes, c(mt.genes, many.genes, cox.genes, rps.genes, rpl.genes)))
}

gene.list <- filter.genes(genes.best4)

##Save gene list and dataframe
write(gene.list, "~/Desktop/human_intestine_data/genes_best4_enriched.txt")
gene.list <- scan("~/Desktop/human_intestine_data/genes_best4_enriched.txt", what = "character")
saveRDS(df, "~/Desktop/human_intestine_data/human_intestine_average_exp_df.rds")

##Crop df to enriched genes and save 
dat <- df[rownames(df) %in% gene.list, ]
saveRDS(dat, "~/Desktop/human_intestine_data/human_intestine_average_exp_best4_enriched.rds")

##Ok now convert gene names to zebrafish homologs
gene.conversion.table <- read.delim("~/Desktop/human_intestine_data/ZF_human_homologs/gene_orthologs/hs_dr_orthologs_ENSEMBL-SYMBOL.tsv")
zf.homologs.best4.human <- unlist(lapply(gene.list, function(x) gene.conversion.table$dr.SYMBOL[gene.conversion.table$hs.SYMBOL == x]))
write(zf.homologs.best4.human, "~/Desktop/human_intestine_data/genes_best4_enriched_human_converted.txt")

##Now subset the zebrafish intestinal cells and get average expression of zebrafish best4+ cells
obj.int <- subset(obj1, idents = c("9", "10", "11", "12", "13", "14", "15", "16", "17", "18"))

##Do the usual normalization, PCA, neighbors, UMAP/TSNE
##obj <- NormalizeData(object = obj)
##obj <- FindVariableFeatures(object = obj)
##obj <- ScaleData(object = obj)
##obj <- RunPCA(obj, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)

# Dimensionality Reduction: UMAP
##message(paste0(Sys.time(), ": UMAP"))
#nn.use <- floor(sqrt(ncol(obj@assays$RNA@counts))/10)*10
##nn.use <- 30
##obj <- RunUMAP(obj, dims=1:30, n.neighbors=nn.use)

# Clustering
#message(paste0(Sys.time(), ": Clustering"))
#obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, nn.eps = 0, k.param = 50)

#obj <- FindNeighbors(object = obj)
#obj <- FindClusters(object = obj, resolution = c(0.75, 1))
#obj <- RunTSNE(object = obj, dims = 1:30, nn.use = 50)

#DimPlot(obj, group.by = "RNA_snn_res.1", label = T)
#Idents(obj) <- "RNA_snn_res.1"

cluster.names <- c("intestine_prog", "EC-1", "EC-2", "EC-3", "post_LREs", "progenitors/late", "goblet", "best4", "EECs", "tuft-like")
names(cluster.names) <- levels(obj.int)
obj.int <- RenameIdents(obj.int, cluster.names)
DimPlot(obj.int, label = T)

##Calculate average expression of all cells within a cluster
cluster.averages <- AverageExpression(obj.int)
head(cluster.averages[["RNA"]])

orig.levels <- levels(obj.int)
Idents(obj.int) <- gsub(pattern = " ", replacement = "_", x = Idents(obj.int))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(obj.int) <- orig.levels
cluster.averages <- AverageExpression(obj.int, return.seurat = TRUE)
cluster.averages

##Plot a scatter plot between the various cell types
CellScatter(cluster.averages, cell1 = "Best4+", cell2 = "Goblet")

##Ok, now put a threshold for Best4+ enterocytes for genes that are expressed above a certain threshold
message(paste0(Sys.time(), ": Convert cluster averages to dataframe"))
df <- as.data.frame(cluster.averages)

##Calculate markers for the different idents in obj for human colon intestinal cells
posMarkers.wilcox <- sapply(levels(Idents(obj)), function(i) FindMarkers(obj,i,assay.type='RNA',test.use="wilcox",only.pos=T, min.cells.group = 1),simplify=F)
posMarkers.roc <- sapply(levels(Idents(obj)), function(i) FindMarkers(obj,i,test.use="roc",only.pos=T),simplify=F)
markers.human <- list(wilcox = posMarkers.wilcox, roc = posMarkers.roc)
saveRDS(markers.human, file=paste0("~/Desktop/human_intestine_data/human_intestine_markers.rds"))
markers.human <- readRDS(file=paste0("~/Desktop/human_intestine_data/human_intestine_markers.rds"))

##Calculate markers for different Idents for zebrafish intestinal cells
posMarkers.wilcox <- sapply(levels(Idents(obj.int)), function(i) FindMarkers(obj.int,i,assay.type='RNA',test.use="wilcox",only.pos=T, min.cells.group = 1),simplify=F)
posMarkers.roc <- sapply(levels(Idents(obj.int)), function(i) FindMarkers(obj.int,i,test.use="roc",only.pos=T),simplify=F)
markers.zf <- list(wilcox = posMarkers.wilcox, roc = posMarkers.roc)
saveRDS(markers.zf, file=paste0("~/Desktop/human_intestine_data/zebrafish_intestine_markers.rds"))
markers.zf <- readRDS(file=paste0("~/Desktop/human_intestine_data/zebrafish_intestine_markers.rds"))

##Find genes in human dataset that has a fold change of more than 0.25
genes.human <- rownames(markers.human$wilcox$`Best4+`)[which(markers.human$wilcox$`Best4+`$avg_log2FC > 0.25)]
dat <- markers.human$wilcox$`Best4+`[genes.human, ]

genes.zf <- rownames(markers.zf$wilcox$best4)[which(markers.zf$wilcox$best4$avg_log2FC > 0.25)]
df <- markers.zf$wilcox$best4[genes.zf, ]

saveRDS(df, "~/Desktop/human_intestine_data/zf_intestine_avgLogFC.rds")
saveRDS(dat, "~/Desktop/human_intestine_data/human_intestine_avgLogFC.rds")
df <- readRDS("~/Desktop/human_intestine_data/zf_intestine_average_exp_df.rds")
dat <- readRDS("~/Desktop/human_intestine_data/human_intestine_avgLogFC.rds")

##Convert all genes to lowercase
dat$genes.lc <- tolower(rownames(dat))

##Save the human gene expression data and convert them to ZF homologs
write.csv(dat, "~/Desktop/human_intestine_data/human_best4_genes_avgLogFC.csv")
write.csv(df, "~/Desktop/human_intestine_data/zf_best4_genes_avgLogFC.csv")

##Load the two spreadhseets
dat <- read.csv("~/Desktop/human_intestine_data/human_best4_genes_avgLogFC.csv")
df <- read.csv("~/Desktop/human_intestine_data/zf_best4_genes_avgLogFC.csv")

##Draw a scatter plot as in Figure 7B
##First merge the two dataframes
##Trim dat and df to two columns
dat.human.colon <- dat[, c("zf.homologs", "avg_log2FC")]
dat.zf <- df[, c("zf.homologs", "avg_log2FC")]
dat.human.si <- dm[, c("zf.homologs", "avg_log2FC")]
dx.1 <- merge(dat.human.colon, dat.zf, by = "zf.homologs", all.y = TRUE, all.x = TRUE)
dx <- merge(dx.1, dat.human.si, by = "zf.homologs", all.y = TRUE, all.x = TRUE)

##Add rownames
rownames(dx) <- dx$zf.homologs
colnames(dx) <- c("genes", "log_exp_colon", "log_exp_zf", "log_exp_SI")


##Convert all NAs to 0
##Replace all NA value to 0
dx[is.na(dx)] <- 0

##Now plot scatter plot
plot(dx$log_exp_SI, dx$log_exp_zf)

##Now plot with genes with different colors - Figure S9F
plot(dx$log_exp_SI, dx$log_exp_zf, main = "Expression in humans versus ZF ",
     xlab="human_BEST4_expression", 
     ylab="ZF_BEST4_expression", 
     cex=0.95, cex.lab= 0.5, title(cex.main=2, font.main=7))
abline(-1, 1, col='purple', lty="dashed")
abline(1, 1, col="purple", lty="dashed")
above <- dx$log_exp_zf - dx$log_exp_SI > 1
points(dx$log_exp_SI[above], dx$log_exp_zf[above], col="red", cex=0.5)
below <- dx$log_exp_zf - dx$log_exp_SI < -1
points(dx$log_exp_SI[below], dx$log_exp_zf[below], col="blue", cex=0.5)


ggplot(dx, aes(x = log_exp_SI, y = log_exp_zf, label = genes)) +
  geom_point(size=2, shape=16) +
  geom_text()


##Now plot the individual scatter plots - ZF best4+ cells vs human colon - Figure S9E

##First merge the two dataframes
##Trim dat and df to two columns
dat.human <- dat[, c("zf.homologs", "avg_log2FC")]
dat.zf <- df[, c("zf.homologs", "avg_log2FC")]
dx <- merge(dat.human, dat.zf, by = "zf.homologs", all.y = TRUE, all.x = TRUE)

##Add rownames
rownames(dx) <- dx$zf.homologs
colnames(dx) <- c("genes", "log_exp_human", "log_exp_zf")

##Convert all NAs to 0
##Replace all NA value to 0
dx[is.na(dx)] <- 0

##Now plot scatter plot
plot(dx$log_exp_human, dx$log_exp_zf)

##Now plot with genes with different colors for overexpressed and underexpressed - Figure S9E
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

##Plot using ggplot2
ggplot(dx, aes(x = log_exp_human, y = log_exp_zf, label = genes)) +
  geom_point(size=2, shape=16) +
  geom_text()

##Now plot the individual scatter plots - ZF best4+ cells vs human small intestine - Figure S9F
##First merge the two dataframes
##Trim dat and df to two columns
dat.human.colon <- dat[, c("zf.homologs", "avg_log2FC")]
dat.zf <- df[, c("zf.homologs", "avg_log2FC")]
dat.human.si <- dm[, c("zf.homologs", "avg_log2FC")]
dx.1 <- merge(dat.human.colon, dat.zf, by = "zf.homologs", all.y = TRUE, all.x = TRUE)
dx <- merge(dx.1, dat.human.si, by = "zf.homologs", all.y = TRUE, all.x = TRUE)

##Add rownames
rownames(dx) <- dx$zf.homologs
colnames(dx) <- c("genes", "log_exp_colon", "log_exp_zf", "log_exp_SI")


##Convert all NAs to 0
##Replace all NA value to 0
dx[is.na(dx)] <- 0

##Now plot scatter plot
plot(dx$log_exp_SI, dx$log_exp_zf)

##Now plot with genes with different colors - Figure S11
plot(dx$log_exp_SI, dx$log_exp_zf, main = "Expression in humans versus ZF ",
     xlab="human_BEST4_expression", 
     ylab="ZF_BEST4_expression", 
     cex=0.95, cex.lab= 0.5, title(cex.main=2, font.main=7))
abline(-1, 1, col='purple', lty="dashed")
abline(1, 1, col="purple", lty="dashed")
above <- dx$log_exp_zf - dx$log_exp_SI > 1
points(dx$log_exp_SI[above], dx$log_exp_zf[above], col="red", cex=0.5)
below <- dx$log_exp_zf - dx$log_exp_SI < -1
points(dx$log_exp_SI[below], dx$log_exp_zf[below], col="blue", cex=0.5)


ggplot(dx, aes(x = log_exp_SI, y = log_exp_zf, label = genes)) +
  geom_point(size=2, shape=16) +
  geom_text()


