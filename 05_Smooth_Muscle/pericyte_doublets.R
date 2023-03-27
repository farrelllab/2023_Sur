##This script is to confirm that the three pericyte populations are actually distinct populations and not doublets with other cell types especially endothelial cells. 

library(Seurat)
library(Matrix)

#Load mama and annotations
message(paste0(Sys.time(), ": Loading mama"))
mama <- readRDS("/data/CSD/zfext/LTA/results/02-Clustering/merged_timepoints/obj/merged_mama_ds_slim_seurat_25Jan2022.rds")

##Add relevant columns to metadata in mama
cell.annot <- read.table("/data/CSD/zfext/LTA/results/02-Clustering/merged_timepoints/meta.data/celltype_annotations_df.tsv", sep = " ", header = T, stringsAsFactors = F)
cell.annot$tissue <- unlist(lapply(strsplit(x = cell.annot$clust, split = "\\."), function(x) x[1]))
cell.annot$clust.sg <- paste0(cell.annot$clust, "_", gsub(" ", "", cell.annot$stage.group))

mama@meta.data <- cbind(
  mama@meta.data[,setdiff(colnames(mama@meta.data), colnames(cell.annot))],
  cell.annot[rownames(mama@meta.data),]
)


##Load the non-skeletal muscle (mural-cell) seurat and markers object
obj.smc <- readRDS("/data/CSD/zfext/results/04f-SubsetsV6/obj_subsets_v6/mural_seurat.rds")
markers <- readRDS("/data/CSD/zfext/results/04f-SubsetsV6/obj_subsets_v6/mural_markers.rds")

##Plot the UMAP
DimPlot(obj.smc, label = T)  ##pericyte-1 = C20, ##pericyte-2 = C4

##Get top 5 markers for each of the three pericyte populations
markers.20 <- rownames(markers$wilcox$`20`[1:5, ])
markers.4 <- rownames(markers$wilcox$`4`[1:5, ])


##Get other cells from the global atlas that expresses these genes

##for pericyte-1, some of the top markers expressed are lmx1bb, ptger4a, and adma - find cells other than pericytes that express these genes to make artificial doublets
cells.lmx1bb <- WhichCells(mama, expression = lmx1bb > 2)
cells.adma <- WhichCells(mama, expression = adma > 2)
cells.ptger4a <- WhichCells(mama, expression = ptger4a > 2)

cells.peri.20 <- setdiff(unlist(unique(list(cells.lmx1bb, cells.adma, cells.ptger4a))), WhichCells(obj.smc, idents = c("20")))

##for pericyte-2, some of the top markers expressed are prox2, ndufa4l2a, and epas1a - find cells other than pericytes that express these genes to make artificial doublets
cells.prox2 <- WhichCells(mama, expression = prox2 > 2)
cells.ndufa4l2a <- WhichCells(mama, expression = ndufa4l2a > 2)
cells.epas1a <- WhichCells(mama, expression = epas1a > 2)

cells.peri.4 <- setdiff(unlist(unique(list(cells.prox2, cells.ndufa4l2a, cells.epas1a))), WhichCells(obj.smc, idents = c("4")))


##Only limit to cells that belong to the stages that the three pericytes belong to
##Find the stages that the cells in each pericyte category belongs to
stg.peri.1 <- table(mama@meta.data[WhichCells(obj.smc, idents = c("20")), "stage.nice"]) # 60 - 120 hpf
stg.peri.2 <- table(mama@meta.data[WhichCells(obj.smc, idents = c("4")), "stage.nice"]) # 60 - 120 hpf
stg.peri.3 <- table(mama@meta.data[WhichCells(obj.smc, idents = c("3")), "stage.nice"]) # 60 - 120 hpf
stg.peri.0 <- table(mama@meta.data[WhichCells(obj.smc, idents = c("9")), "stage.nice"]) # 60 - 120 hpf

##Subset from the above cell lists to limit to stages above 60 hrs
cells.peri.1.above.60 <- cells.peri.20[cells.peri.20 %in% rownames(mama@meta.data)[which(mama@meta.data$stage.nice > 60)]]
cells.peri.2.above.60 <- cells.peri.4[cells.peri.4 %in% rownames(mama@meta.data)[which(mama@meta.data$stage.nice > 60)]]
cells.peri.3.above.60 <- cells.peri.3[cells.peri.3 %in% rownames(mama@meta.data)[which(mama@meta.data$stage.nice > 60)]]

##To create random doubets between these cells and pericyte-0, randomize the pool of these potential doublet cells
n <- 5000
cells.peri.1.random <- sample(cells.peri.1.above.60, n, replace = T)
cells.peri.2.random <- sample(cells.peri.2.above.60, n, replace = T)
cells.peri.3.random <- sample(cells.peri.3.above.60, n, replace = T)
cells.peri.0.random <- sample(WhichCells(obj.smc, idents = c("9")), n, replace = T)

##First subset the data slot for just the cells to be used to create artificial doublets
#Subset the matrices such that the union of variable genes on the global atlas and the non-skeletal muscle atlas are used
##Get  union of variable genes
genes.to.use <- unlist(unique(list(mama@assays$RNA@var.features, obj.smc@assays$RNA@var.features)))

mat.peri.20 <- mama@assays$RNA@data[genes.to.use, which(colnames(mama@assays$RNA@data) %in% cells.peri.20)]  ##Note that the resulting matrix will only contain the rows that have at least one non-zero value in the selected columns.
mat.peri.4 <- mama@assays$RNA@data[genes.to.use, which(colnames(mama@assays$RNA@data) %in% cells.peri.4)]
mat.peri.3 <- mama@assays$RNA@data[genes.to.use, which(colnames(mama@assays$RNA@data) %in% cells.peri.3)]

##Also subset the data slot for the pericyte-0 cells 
mat.peri.0 <- mama@assays$RNA@data[genes.to.use, which(colnames(mama@assays$RNA@data) %in% cells.peri.0.random)]

##Antilog expression values
mat.peri.20 <- exp(mat.peri.20)
mat.peri.4 <- exp(mat.peri.4)
mat.peri.3 <- exp(mat.peri.3)
mat.peri.0 <- exp(mat.peri.0)


##Merge the potential doublet matrix with the pericyte-0 matrix
mat.peri.20.doublet <- cbind(mat.peri.20, mat.peri.0)
mat.peri.4.doublet <- cbind(mat.peri.4, mat.peri.0)
mat.peri.3.doublet <- cbind(mat.peri.3, mat.peri.0)


##Take mean of pairs of columns between the two subsets - pericyte-0 versus the others
# Create a matrix to store the averaged expression values
avg_mat.peri.1 <- matrix(0, nrow=5000, ncol=n)

# Loop over each column in each of the expression matrices
# iterate over each column of the matrix of cells expressing pericyte-1 genes and pericyte 0 cells
for (i in 1:n){
    # calculate the mean expression for each gene
    avg_col_peri.20 <- rowMeans(mat.peri.20.doublet[, c(cells.peri.1.random[i], cells.peri.0.random[i])])
    # assign the mean expression values to the corresponding columns of avg_mat
    avg_mat.peri.1[, i] <- avg_col_peri.20
}


# Loop over each column in each of the expression matrices
# iterate over each column of the 
avg_mat.peri.2 <- matrix(0, nrow=5000, ncol=n)
for (i in 1:n){
  # calculate the mean expression for each gene
  avg_col_peri.4 <- rowMeans(mat.peri.4.doublet[, c(cells.peri.2.random[i], cells.peri.0.random[i])])
  # assign the mean expression values to the corresponding columns of avg_mat
  avg_mat.peri.2[, i] <- avg_col_peri.4
}


##Now log each of the matrices
mat.peri.20.log.exp <- log(avg_mat.peri.1)
mat.peri.4.log.exp <- log(avg_mat.peri.2)
mat.peri.3.log.exp <- log(avg_mat.peri.3)

colnames(mat.peri.20.log.exp) <- cells.peri.1.random
colnames(mat.peri.4.log.exp) <- cells.peri.2.random
colnames(mat.peri.3.log.exp) <- cells.peri.3.random

saveRDS(mat.peri.20.log.exp, "/data/CSD/zfext/results/pericyte_doublets/pericyte_1_doublet_mat.rds")
saveRDS(mat.peri.4.log.exp, "/data/CSD/zfext/results/pericyte_doublets/pericyte_2_doublet_mat.rds")
saveRDS(mat.peri.3.log.exp, "/data/CSD/zfext/results/pericyte_doublets/pericyte_3_doublet_mat.rds")

##To find if pericyte-1, -2, and -3 are similar in gene expression to the artificial doublets with pericyte-0, a KNN neighbor matrix will be created
##Get the pericyte-1 gene expression matrix
mat.peri.1 <- mama@assays$RNA@data[genes.to.use, which(colnames(mama@assays$RNA@data) %in% WhichCells(obj.smc, idents = c("20")))]
mat.peri.2 <- mama@assays$RNA@data[genes.to.use, which(colnames(mama@assays$RNA@data) %in% WhichCells(obj.smc, idents = c("4")))]
mat.peri.3 <- mama@assays$RNA@data[genes.to.use, which(colnames(mama@assays$RNA@data) %in% WhichCells(obj.smc, idents = c("3")))]

# Merge the two matrices for each population and keep the common genes
# Combine the two matrices and keep only common genes
merged_expr_peri.1 <- cbind(mat.peri.20.log.exp, mat.peri.1)
merged_expr_peri.2 <- cbind(mat.peri.4.log.exp, mat.peri.2)
merged_expr_peri.3 <- cbind(mat.peri.3.log.exp, mat.peri.3)

##Save the merged matrices
saveRDS(merged_expr_peri.1, "/data/CSD/zfext/results/pericyte_doublets/pericyte_1_merged_exp_mat.rds")
saveRDS(merged_expr_peri.2, "/data/CSD/zfext/results/pericyte_doublets/pericyte_2_merged_exp_mat.rds")
saveRDS(merged_expr_peri.3, "/data/CSD/zfext/results/pericyte_doublets/pericyte_3_merged_exp_mat.rds")

# Extract the distance matrix
dist_matrix.peri.1 <- dist(t(merged_expr_peri.1))
dist_matrix.peri.2 <- dist(t(merged_expr_peri.2))
dist_matrix.peri.3 <- dist(t(merged_expr_peri.3))

saveRDS(knn_dist_peri.1, "/data/CSD/zfext/results/pericyte_doublets/pericyte_1_knn_dists.rds")
saveRDS(knn_dist_peri.2, "/data/CSD/zfext/results/pericyte_doublets/pericyte_2_knn_dists.rds")
saveRDS(knn_dist_peri.3, "/data/CSD/zfext/results/pericyte_doublets/pericyte_3_knn_dists.rds")

saveRDS(dist_matrix.peri.1, "/data/CSD/zfext/results/pericyte_doublets/pericyte_1_dist_mat.rds")
saveRDS(dist_matrix.peri.2, "/data/CSD/zfext/results/pericyte_doublets/pericyte_2_dist_mat.rds")
saveRDS(dist_matrix.peri.3, "/data/CSD/zfext/results/pericyte_doublets/pericyte_3_dist_mat.rds")

##Convert to matrix
matrix.peri.1 <- as.matrix(dist_matrix.peri.1)
matrix.peri.2 <- as.matrix(dist_matrix.peri.2)
matrix.peri.3 <- as.matrix(dist_matrix.peri.3)

##Find distances between each pericyte groups and between pericyte to doublets
dist.p1.p1 <- matrix.peri.1[colnames(mat.peri.1), colnames(mat.peri.1)]
dist.p1.doublet <- matrix.peri.1[colnames(mat.peri.1), colnames(mat.peri.20.log.exp)]

hist(dist.p1.p1)
hist(dist.p1.doublet)

##Plot the distance matrices between pericyte 1 to pericyte 1 cells and pericyte 1 to pericyte 1 related doublet cells
dist.p1.p1 <- matrix.peri.1[WhichCells(obj.smc, idents = c("20")), WhichCells(obj.smc, idents = c("20"))]
dist.p1.doublet <- matrix.peri.1[WhichCells(obj.smc, idents = c("20")), colnames(peri.1.doublet.mat)]

dist.p2.p2 <- matrix.peri.2[WhichCells(obj.smc, idents = c("4")), WhichCells(obj.smc, idents = c("4"))]
dist.p2.doublet <- matrix.peri.2[WhichCells(obj.smc, idents = c("4")), colnames(peri.2.doublet.mat)]

dist.p3.p3 <- matrix.peri.3[WhichCells(obj.smc, idents = c("3")), WhichCells(obj.smc, idents = c("3"))]
dist.p3.doublet <- matrix.peri.3[WhichCells(obj.smc, idents = c("3")), colnames(peri.3.doublet.mat)]

##Now plot the distribution of distances for these different combinations
#For pericyte-1, plot the distribution of distances between the actual pericyte-1 cells
library(ggplot2)
library(cowplot)

##Convert pericyte-1 to pericyte-1 and doublet distances to a dataframe
dist.p1.p1.df <- data.frame(
  distance = c(dist.p1.p1, dist.p1.doublet),
  group = c(rep("p1_to_p1", length(dist.p1.p1)), rep("p1_to_doublet", length(dist.p1.doublet)))
)
                            
# create histogram for distances between cells in pericyte 1 cluster and pericyte-1 cluster to artificial doublets
# Plot histogram of distances for p1 to p1 and p1 to doublets
ggplot(dist.p1.p1.df, aes(x = distance, fill = group)) + 
  geom_histogram(alpha = 0.5, bins = 50) +
  scale_fill_manual(values = c("blue", "red")) +
  xlab("Distance") + ylab("Count") +
  ggtitle("Distance between cells of pericyte-1 and artificial doublets")

ggsave("pericyte-1_doublets_distance.pdf", path = "~/Box/zfext/annotations_celltype_curated_newMama/mural/pericyte_doublets/")

##Convert pericyte-2 to pericyte-2 and doublet distances to a dataframe
dist.p2.p2.df <- data.frame(
  distance = c(dist.p2.p2, dist.p2.doublet),
  group = c(rep("p2_to_p2", length(dist.p2.p2)), rep("p2_to_doublet", length(dist.p2.doublet)))
)

# create histogram for distances between cells in pericyte 1 cluster and pericyte-1 cluster to artificial doublets
# Plot histogram of distances for p1 to p1 and p1 to doublets
ggplot(dist.p2.p2.df, aes(x = distance, fill = group)) + 
  geom_histogram(alpha = 0.5, bins = 50) +
  scale_fill_manual(values = c("blue", "red")) +
  xlab("Distance") + ylab("Count") +
  ggtitle("Distance between cells of pericyte-2 and artificial doublets")

ggsave("pericyte-2_doublets_distance.pdf", path = "~/Box/zfext/annotations_celltype_curated_newMama/mural/pericyte_doublets/")

