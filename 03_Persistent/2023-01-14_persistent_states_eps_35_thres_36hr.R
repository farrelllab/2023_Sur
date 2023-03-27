library(Seurat)
library(URD)

save.path <- "~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data_thres_36h/"
plot.path <- "~/Box/zfext/global_analysis_curated/transient_vs_long-term/plots_neighbors_data_thres_36h/"
base.path <- "~/Box/zfext/global_analysis_curated/transient_vs_long-term/"

#Check how many  NAs are present in the hematopoietic dataset at eps = 34 and 35
sample <- "hematopoietic"
eps <- 34
nearest.hema.34 <- readRDS(paste0(save.path, sample, "_mean_stage_diff_eps_", eps, ".rds"))
nn <- readRDS(paste0(save.path, sample, "_neighbor_matrix_", eps, "_binary.rds"))

eps <- 35
nearest.hema.35 <- readRDS(paste0(save.path, sample, "_mean_stage_diff_eps_", eps, ".rds"))
table(is.na(nearest.hema.34)) ##3368 NAs
table(is.na(nearest.hema.35)) ##2811 NAs

#Check how many  NAs are present in the glial dataset at eps = 34 and 35
sample <- "glial"
eps <- 34
nearest.glia.34 <- readRDS(paste0(save.path, sample, "_mean_stage_diff_eps_", eps, ".rds"))

eps <- 35
nearest.glia.35 <- readRDS(paste0(save.path, sample, "_mean_stage_diff_eps_", eps, ".rds"))

eps <- 39
nearest.glia.39 <- readRDS(paste0(save.path, sample, "_mean_stage_diff_eps_", eps, ".rds"))
table(is.na(nearest.glia.34)) ##34% don't have neighbors
table(is.na(nearest.glia.35)) ##29% don't have neighbors
table(is.na(nearest.glia.39)) ##19% don't have neighbors

#Check how many  NAs are present in the glial dataset at eps = 34 and 35
sample <- "neural"
eps <- 34
nearest.neural.34 <- readRDS(paste0(save.path, sample, "_mean_stage_diff_eps_", eps, ".rds"))

eps <- 35
nearest.neural.35 <- readRDS(paste0(save.path, sample, "_mean_stage_diff_eps_", eps, ".rds"))

eps <- 39
nearest.neural.39 <- readRDS(paste0(save.path, sample, "_mean_stage_diff_eps_", eps, ".rds"))
table(is.na(nearest.neural.34)) ##34% don't have neighbors
table(is.na(nearest.neural.35)) ##29% don't have neighbors


#Check how many  NAs are present in the eye dataset at eps = 34 and 35
sample <- "eye"
eps <- 34
nearest.eye.34 <- readRDS(paste0(save.path, sample, "_mean_stage_diff_eps_", eps, ".rds"))

eps <- 35
nearest.eye.35 <- readRDS(paste0(save.path, sample, "_mean_stage_diff_eps_", eps, ".rds"))

eps <- 39
nearest.eye.39 <- readRDS(paste0(save.path, sample, "_mean_stage_diff_eps_", eps, ".rds"))

table(is.na(nearest.eye.34)) ##33% don't have neighbors
table(is.na(nearest.eye.35)) ## 29% don't have neighbors
table(is.na(nearest.eye.39)) ## 19% dont't have neighbors


##Ok at eps 34, 33% cells in glia don't have neighbors, at eps 35, 29% cells don't have neighbors. Close to 40, most cells have neighbors. 
##So depends on the dataset how many cells don't have neighbors.
##Ok so, at eps = 35, atleast 70% of most datasets have neighbors, hence we pick eps = 35

sample <- "axial"
eps <- "35"
nearest.axial <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.axial <- as.data.frame(nearest.axial)
colnames(nearest.axial) <- "stage.diff"

sample <- "epidermis"
nearest.basal <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.basal <- as.data.frame(nearest.basal)
colnames(nearest.basal) <- "stage.diff"

sample <- "blastomeres"
nearest.blastula <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.blastula <- as.data.frame(nearest.blastula)
colnames(nearest.blastula) <- "stage.diff"

sample <- "endoderm"
nearest.endo <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.endo <- as.data.frame(nearest.endo)
colnames(nearest.endo) <- "stage.diff"

sample <- "eye"
nearest.eye <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.eye <- as.data.frame(nearest.eye)
colnames(nearest.eye) <- "stage.diff"

sample <- "fin"
nearest.fin <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.fin <- as.data.frame(nearest.fin)
colnames(nearest.fin) <- "stage.diff"

sample <- "glial"
nearest.glia <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.glia <- as.data.frame(nearest.glia)
colnames(nearest.glia) <- "stage.diff"

sample <- "hematopoietic"
nearest.hema <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.hema <- as.data.frame(nearest.hema)
colnames(nearest.hema) <- "stage.diff"

sample <- "ionocytes"
nearest.iono <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.iono <- as.data.frame(nearest.iono)
colnames(nearest.iono) <- "stage.diff"

sample <- "mesenchyme"
nearest.mese <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.mese <- as.data.frame(nearest.mese)
colnames(nearest.mese) <- "stage.diff"

sample <- "mural"
nearest.mural <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.mural <- as.data.frame(nearest.mural)
colnames(nearest.mural) <- "stage.diff"

sample <- "muscle"
nearest.muscle <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.muscle <- as.data.frame(nearest.muscle)
colnames(nearest.muscle) <- "stage.diff"

sample <- "neural"
nearest.neural <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.neural <- as.data.frame(nearest.neural)
colnames(nearest.neural) <- "stage.diff"

sample <- "otic"
nearest.otic <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.otic <- as.data.frame(nearest.otic)
colnames(nearest.otic) <- "stage.diff"

sample <- "periderm"
nearest.periderm <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.periderm <- as.data.frame(nearest.periderm)
colnames(nearest.periderm) <- "stage.diff"

sample <- "pgc"
nearest.pgc <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.pgc <- as.data.frame(nearest.pgc)
colnames(nearest.pgc) <- "stage.diff"

sample <- "pigment"
nearest.pigment <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.pigment <- as.data.frame(nearest.pigment)
colnames(nearest.pigment) <- "stage.diff"

sample <- "pronephros"
nearest.pron <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.pron <- as.data.frame(nearest.pron)
colnames(nearest.pron) <- "stage.diff"

sample <- "taste"
nearest.taste <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.taste <- as.data.frame(nearest.taste)
colnames(nearest.taste) <- "stage.diff"

#Combine all the dataframes - cell ids and their stage differences. 
nearest.total <- rbind(nearest.axial, nearest.basal, nearest.blastula, nearest.endo, nearest.eye, nearest.fin, nearest.glia, nearest.hema,
                       nearest.iono, nearest.mese, nearest.mural, nearest.muscle, nearest.neural, nearest.otic, nearest.periderm, nearest.pgc, 
                       nearest.pigment, nearest.pron, nearest.taste)

cells.all <- WhichCells(mama)
cells.rest <- setdiff(rownames(mama@meta.data), rownames(nearest.total))
cells.rest <- rownames(mama@meta.data)[mama@meta.data$tissue.subsets == "cephalic"]

##Load cephalic seurat object
obj <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/cephalic/cephalic_seurat.rds")
DimPlot(obj, label = T)
cells.ceph <- WhichCells(obj)

dist.mat  <- dist(t(obj@assays$RNA@data))
mat <- as.matrix(dist.mat)
for(eps in c(31, 32, 33, 34, 35, 36)){
  nearest <- sapply(rownames(mat), function(cell) mean(abs(obj@meta.data[setdiff(colnames(mat)[which(mat[cell,] <= eps)], cell), "stage.nice"] - obj@meta.data[cell, "stage.nice"])))
  nn <- sapply(rownames(mat), function(cell) setdiff(colnames(mat)[which(mat[cell, ] <= eps)], cell))
  #neighbors <- list(nearest, nn)
  #names(neighbors) <- c("mean.stage.diff", "neighbors")
  
  message(paste0(Sys.time(), ": Saving"))
  saveRDS(nearest, file = paste0(save.path, sample, "_mean_stage_diff_eps_", eps, ".rds"))
  saveRDS(nn, file = paste0(save.path, sample, "_neighbors_eps_", eps, ".rds"))
}

sample <- "cephalic"
eps <- 35
nearest.cephalic <- readRDS(file = paste0("~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data/", sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.cephalic <- as.data.frame(nearest.cephalic)
colnames(nearest.cephalic) <- "stage.diff"
nearest.all <- rbind(nearest.total, nearest.cephalic)  

cells.out <- rownames(nearest.all)[which(rownames(nearest.all) %in% rownames(mama@meta.data))]
nearest.clean <- subset(nearest.all, cells = cells.extra)
nearest.clean <- as.data.frame(nearest.all[(rownames(nearest.all) %in% cells.out), ])
rownames(nearest.clean) <- cells.out
colnames(nearest.clean) <- "stage.diff"
cells.extra <- setdiff(rownames(nearest.all), rownames(mama@meta.data))

mama@meta.data <- cbind(mama@meta.data, nearest.clean[rownames(mama@meta.data), ])
colnames(mama@meta.data)[77] <- c("stage.diff")
##Get rid of all NaN values
mama@meta.data[which(mama@meta.data$stage.diff == "NaN"), "stage.diff"] <- 0
saveRDS(mama@meta.data$stage.diff, "~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_stage.diff_data_eps_35_with_neural.rds")
stage.diff <- readRDS("~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_stage.diff_data_eps_35_with_neural.rds")
mama@meta.data$stage.diff <- stage.diff

DimPlot(mama, group.by = "stage.diff", raster = F) + NoAxes() + NoLegend() + scale_fill_continuous(type = "gradient")
library(RColorBrewer)
pdf("~/Desktop/stage.diff.pdf")
DimPlot(mama, group.by = "stage.diff", raster = F) + NoAxes() + NoLegend() + scale_fill_continuous(type = "gradient")
dev.off()

##Categorize the different cells into how long they persist for using this eps setting
##Define cells that persist for different lengths of time
cells.under.24 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff < 24)]
cells.24.36 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff >= 24 & mama@meta.data$stage.diff < 36)]
cells.36.48 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff >= 36 & mama@meta.data$stage.diff < 48)]
#cells.36.48 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff >= 36 & mama@meta.data$stage.diff < 48)]
#cells.48.60 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff >= 48 & mama@meta.data$stage.diff < 60)]
#cells.60.72 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff >= 60 & mama@meta.data$stage.diff < 72)]
#cells.72.84 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff >= 72 & mama@meta.data$stage.diff < 84)]
cells.48 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff >= 48)]

##Now assign a slot in mama and plot these cells
mama@meta.data$long_vs_short <- NA
mama@meta.data[cells.under.24, "long_vs_short"] <- "<24hr"
mama@meta.data[cells.24.36, "long_vs_short"] <- "24-36hr"
mama@meta.data[cells.36.48, "long_vs_short"] <- "36-48hr"
mama@meta.data[cells.48, "long_vs_short"] <- ">48hr"

colors.stg <- setNames(c("#479D8A", "#74B9D3", "#E7AF3D", "#E98A33"), c("<24hr", "24-36hr", "36-48hr",  ">48hr"))
dpi <- 300
png("~/Box/zfext/global_analysis_curated/transient_vs_long-term/figures/mama_transient_vs_persistent_states_data_eps_35_v4.png", height = 5*dpi, width = 5*dpi)
DimPlot(mama, group.by = "long_vs_short", cols = colors.stg, raster = F, na.value = "#f5f7f6") + NoAxes() + NoLegend()
dev.off()

saveRDS(mama@meta.data$long_vs_short, "~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_long_vs_short_stage_diff_data_eps_35.rds")
long_vs_short <- readRDS("~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_long_vs_short_stage_diff_data_eps_35.rds")
mama@meta.data$long_vs_short <- long_vs_short

##Divide into groups that includes cycling and non-cycling
cc.score <- readRDS("~/Box/zfext/02-Clustering/mama+ds_integrated/merged_mama_cellCycle_score.rds")
mama@meta.data <- cbind(mama@meta.data, cc.score)


mama@meta.data$stage.diff.cc <- NA
mama@meta.data$stage.diff.cc <- paste0(mama@meta.data$long_vs_short, "_", mama@meta.data$cc.status)
color.stg.cc <- setNames(c("#CCCCFF", "#E9967A", "#FF8000", "#FD1212", "#86B0DA", "#08519c", "#40E0D0", "#3B5A57", "#B2FF66", "#009900", "#CCCCFF", "#4728E6", "#FFCCFF", "#CC00CC", "#696969", "#404040"), c("<24hr_cycling", "<24hr_non-cycling", "24-36hr_cycling", "24-36hr_non-cycling", "36-48hr_cycling", "36-48hr_non-cycling", ">48hr_cycling", ">48hr_non-cycling"))

DimPlot(mama, group.by = "stage.diff.cc")

##Just plot the cycling states
color.stg.cycling <- setNames(c("#479D8A", "#74B9D3", "#E7AF3D", "#E98A33"), c("<24hr_cycling", "24-36hr_cycling", "36-48hr_cycling", ">48hr_cycling"))

color.stg.non.cycling <- setNames(c("#479D8A", "#74B9D3", "#E7AF3D", "#E98A33"), c("<24hr_non-cycling", "24-36hr_non-cycling", "36-48hr_non-cycling", ">48hr_non-cycling"))

dpi <- 300
png("~/Box/zfext/global_analysis_curated/transient_vs_long-term/figures/mama_transient_vs_persistent_states_cc.state_eps_40_v4.png", height = 5*dpi, width = 5*dpi)
DimPlot(mama, group.by = "stage.diff.cc", cols = color.stg.cc, raster = F, na.value = "#f5f7f6") + NoAxes()
dev.off()

dpi <- 300
png("~/Box/zfext/global_analysis_curated/transient_vs_long-term/figures/mama_transient_vs_persistent_states_only_cycling_eps_35_v4.png", height = 5*dpi, width = 5*dpi)
DimPlot(mama, group.by = "stage.diff.cc", cols = color.stg.cycling, raster = F, na.value = "#f5f7f6") + NoAxes() + NoLegend()
dev.off()

dpi <- 300
png("~/Box/zfext/global_analysis_curated/transient_vs_long-term/figures/mama_transient_vs_persistent_states_only_non-cycling_eps_35_v4.png", height = 5*dpi, width = 5*dpi)
DimPlot(mama, group.by = "stage.diff.cc", cols = color.stg.non.cycling, raster = F, na.value = "#f5f7f6") + NoAxes() + NoLegend()
dev.off()

##Find the clusters that encompass the long term states for a threshold of 36 hours
dat.annot <- readRDS(paste0(base.path, "long-term_cycling_36h_thres_eps_35_df.rds"))
dat <- readRDS(paste0(base.path, "long-term_cycling_36h_thres_eps_35_df_cell_ids.rds"))

##Get the long-term states at a time-threshold of 36 hours ## Done on the cluster
time.thresh <- 36

all.clusts <- unique(dat.annot$clust)

# Identify proportion of each cluster that is in long-term cycling or long-term non-cycling state
cluster.proportions <- t(as.data.frame(lapply(all.clusts, function(clust) {
  # Get cells from this cluster
  da <- dat.annot[dat.annot$clust ==  clust,]
  # Identify which are 'long-term'
  da$long <- da$stage.diff >= time.thresh
  da$cycling <- da$cc.status == "cycling"
  # Identify portion of cluster that is in each state
  x <- c(
    length(which(!da$long & da$cycling)),
    length(which(da$long & da$cycling)),
    length(which(!da$long & !da$cycling)),
    length(which(da$long & !da$cycling))
  )
  return (x/nrow(da))
})))
colnames(cluster.proportions) <- c("short_cycling", "long_cycling", "short_non-cycling", "long_non-cycling")
rownames(cluster.proportions) <- all.clusts
cluster.proportions <- as.data.frame(cluster.proportions)

# Idemtify clusts with long-term cycling or long-term non-cycling states
thresh.states <- 0.10
clusts.ltc <- rownames(cluster.proportions)[cluster.proportions$long_cycling >= thresh.states]
clusts.ltnc <- rownames(cluster.proportions)[cluster.proportions$`long_non-cycling` >= thresh.states]

# Only operate tissues that have a long-term state
clusts.with.lt <- unique(unlist(list(clusts.ltc, clusts.ltnc)))
tissues.with.lt <- sort(unique(URD::cutstring(c(clusts.ltc, clusts.ltnc), delimiter = "\\.", field = 1)))

# Identify members of each cluster
members <- lapply(clusts.ltc, function(c) rownames(dat)[which(dat$clust == c & dat$states.above.36 == "long-term_cycling")])
names(members) <- clusts.ltc

##Identify neighbors for these clusters that are long-term cycling
# Initialize empty list to hold neighbors results across all tissues
neighbors <- vector("list", length(clusts.ltc))
names(neighbors) <- clusts.ltc

# Loop per tissue because neighbor lists are per-tissue
tissues.with.lt <- c("axial", "endoderm", "epidermis", "eye", "fin", "glial", "hematopoietic", "ionocytes", "mesenchyme", "muscle", "neural", "otic", "pigment", "pronephros", "periderm", "taste")
for (tissue in tissues.with.lt) {
  # Load nn file
  nn <- readRDS(paste0(save.path, tissue, "_neighbor_matrix_35_binary.rds")) #Fill this part in.
  
  # Identify clusters from this tissue
  tissue.short <- substr(tissue, 1, 4)
  clusts.do <- grep(tissue.short, clusts.ltc, value = T)
  
  # Identify the neighbors of each cluster
  for (c in clusts.do) {
    neighbors[[c]] <- names(which(colSums(nn[members[[c]], ]) > 0))
  }
}

saveRDS(neighbors, paste0(base.path, "long-term_clust_with_neighbors_eps_35_all.rds"))
neighbors <- readRDS("~/Box/zfext/global_analysis_curated/transient_vs_long-term/long-term_clust_with_neighbors_eps_35_all_pct_10.rds")

clusts.neighbors <- unique(unlist(lapply(names(neighbors), function(x) dat[neighbors[[x]], "clust"])))
all.clusts <- unlist(list(clusts.neighbors, clusts.ltc))

# Create an overlap matrix
overlap.ltc <- as.data.frame(combn(clusts.ltc, 2)) # may need t()
overlap.ltc <- t(overlap.ltc)
colnames(overlap.ltc) <- c("clust1", "clust2")
overlap.df <- overlap.ltc
overlap.df$overlap <- unlist(lapply(overlap.df, function(x) {
  pair.memb <- unlist(members[x])
  pair.neighb <- unique(unlist(neighbors[x]))
  memb.in.neighb <- intersect(pair.memb, pair.neighb)
  prop <- length(pair.memb) / length(memb.in.neighb)
  return(prop)
}))


##Calculate overlap of clusters with the clusters that represent long-term cells
neighbor.clust <- vector("list", length(clusts.ltc))
names(neighbor.clust) <- clusts.ltc

##Subset dat
##Ok take one set of neighbors - hema.6
#nbr <- neighbors$hema.6

##Subset dat dataframe
#dat.sub <- dat[nbr, ]
#dat.sub <- dat.sub[which(dat.sub$clust %in% clusts.with.lt), ]

##Now get the number of cells per clusters that the neighbors belong to 
#clust.hema.6 <- as.data.frame(table(dat.sub$clust))
#colnames(clust.hema.6) <- c("neighbor.clust", "num_nn")

##Get total number of cells per cluster
#dat.subset <- dat[which(dat$tissue.subsets == tissue), ]
#total.hema.6 <- as.data.frame(table(dat.subset$clust)) 
#colnames(total.hema.6) <- c("clust", "total_num")
#total.hema.6.sub <- total.hema.6[which(total.hema.6$clust %in% clust.hema.6$neighbor.clust), ]
#clust.hema.6 <- cbind(clust.hema.6, total.hema.6.sub$total_num)
#colnames(clust.hema.6) <- c("neighbor.clust", "num_nn", "total_nn")
#clust.hema.6$nn_prop <- clust.hema.6$num_nn/clust.hema.6$total_nn
#clust.hema.6$ltc_clust <- "hema.6"

#Filter hema.6  neighbors such that all clusters where less than 50% of cells are detected as neighbors are discarded
#clust.hema.6.lt <- clust.hema.6[which(clust.hema.6$nn_prop >= 0.75), ]

##Ok now add the identity_super columnns from the annotations
#clust.hema.6.lt$clust <- clust.hema.6.lt$neighbor.clust
#clust.hema.6.lt_annot <- merge(clust.hema.6.lt, annot.confirmed[, c("clust", "identity.super", "identity.sub")], by = "clust")

##Now do this for all the long-term clusters across all tissues
dx <- data.frame()
dat.neighbors <- data.frame()
for(tissue in tissues.with.lt){
  # Identify clusters from this tissue
  tissue.short <- substr(tissue, 1, 4)
  clusts.do <- grep(tissue.short, clusts.ltc, value = T)
  dx <- data.frame()
  for(c in clusts.do){
    ##Get neighbors for the relevant cluster with long-term cycling states within tissue
    nbr <- neighbors[[c]]
    
    ##Subset main dataframe
    dat.sub <- dat[nbr, ]
    dat.sub <- dat.sub[which(dat.sub$clust %in% clusts.with.lt), ]
    
    ##Now get the number of neighbors per cluster within the tissue
    clust.subset <- as.data.frame(table(dat.sub$clust))
    colnames(clust.subset) <- c("neighbor.clust", "num_nn")
    
    ##Get total number of cells per cluster
    dat.subset <- dat[which(dat$tissue.subsets == tissue), ]
    tissue.df <- as.data.frame(table(dat.subset$clust)) 
    colnames(tissue.df) <- c("clust", "total_num")
    
    ##Limit the total number of clusters to just the long-term clusters to calculate proportion
    tissue.df.sub <- tissue.df[which(tissue.df$clust %in% clust.subset$neighbor.clust), ]
    clust.subset <- cbind(clust.subset, tissue.df.sub$total_num)
    colnames(clust.subset) <- c("neighbor.clust", "num_nn", "total_nn")
    clust.subset$nn_prop <- clust.subset$num_nn/clust.subset$total_nn
    clust.subset$ltc_clust <- c
    
    dx <- rbind(dx, clust.subset)
  }
  dat.neighbors <- rbind(dat.neighbors, dx)
}

write.table(dat.neighbors, "~/Box/zfext/global_analysis_curated/transient_vs_long-term/neighbor_clusts_long-term_states_eps_35_thres_36h_pct_10.tsv", sep = "\t")

dat.neighbors <- read.table("~/Box/zfext/global_analysis_curated/transient_vs_long-term/neighbor_clusts_long-term_states_eps_35_thres_36h_pct_10.tsv", sep = "")

##Now just limit the neighbors to long-term cycling states to find overlap between the neighbors and the long-term states
dat.neighbors.ltc <- dat.neighbors[which(dat.neighbors$neighbor.clust %in% clusts.ltc), ]
dat.neighbors.ltc.merge <- dat.neighbors.ltc[which(dat.neighbors.ltc$nn_prop >= 0.99), ]
nbr.merge <- unlist(unique(as.character(dat.neighbors.ltc.merge$neighbor.clust)))
ltc.merge <- unique(dat.neighbors.ltc.merge[dat.neighbors.ltc.merge$neighbor.clust %in% nbr.merge, "ltc_clust"])
clusts.merge <- sort(unique(unlist(list(nbr.merge, ltc.merge))))


##Clusters to merge and remove - keep one for each tissue category and remove the rest
clusts.keep.all <- list()
for(tissue in tissues.with.lt){
  tissue.short <- substr(tissue, 1, 4)
  clusts.do <- grep(tissue.short, clusts.merge, value = T)
  clust.keep <- clusts.do[1]
  clusts.keep.all <- unique(unlist(c(clusts.keep.all, clust.keep)))
  clusts.keep.all <- clusts.keep.all[!(is.na(clusts.keep.all))]
}

clusters.remove <- setdiff(clusts.merge, clusts.keep.all)

##Before doing the merge, make a copy of the original member and neighbor variables
members.orig <- members
neighbors.orig <- neighbors

##Tissues in which clusters need to be merged
clusts.keep <- list()
tissues.merge <- c("epidermis", "eye", "hematopoietic", "ionocytes", "otic", "taste")
##Now merge the cells present in these members and neighbors
for(tissue in tissues.merge){
  # Need to merge the members and neighbors within each tissue
  tissue.short <- substr(tissue, 1, 4)
  ##Find clusters within the ones that need to be merged
  clusts.to.merge <- grep(tissue.short, clusts.merge, value = T)
  
  ##Get members present in these clusters
  cells.members <- unique(unlist(lapply(clusts.to.merge, function(x) members[[x]])))
 cells.neighbors <- unique(unlist(lapply(clusts.to.merge, function(x) neighbors[[x]])))
 
 clust.keep <- clusts.to.merge[[1]]
 clusts.remove <- setdiff(clusts.to.merge, clusts.keep)
 
 members[[clust.keep]] <- cells.members
 neighbors[[clust.keep]] <- cells.neighbors
   
 for(clust.do in clusts.remove){
   members[[clust.do]] <- NA
   neighbors[[clust.do]] <- NA
 }
}
  
##Remove the member and neighbor entries that have NA values
members <- members[!is.na(members)]
neighbors <- neighbors[!is.na(neighbors)]

##Save the new neighbors and members object after merging the different clusters
saveRDS(neighbors, "~/Box/zfext/global_analysis_curated/transient_vs_long-term/long-term_clust_with_neighbors_eps_35_merged_pct_0.99.rds")
saveRDS(members, "~/Box/zfext/global_analysis_curated/transient_vs_long-term/long-term_clust_members_eps_35_merged_pct_0.99.rds")

neighbors <- readRDS("~/Box/zfext/global_analysis_curated/transient_vs_long-term/long-term_clust_with_neighbors_eps_35_merged_pct_0.99.rds")
members <- readRDS("~/Box/zfext/global_analysis_curated/transient_vs_long-term/long-term_clust_members_eps_35_merged_pct_0.99.rds")

##Subset the original dataframe such that the merged clusters are removed
dat.neighbors <- dat.neighbors[-which(dat.neighbors$ltc_clust %in% clusters.remove), ]

#Filter neighbors such that all clusters where less than 75% of cells are detected as neighbors are discarded
dat.neighbors.lt <- dat.neighbors[which(dat.neighbors$nn_prop >= 0.5), ]

##From our annotations we identified several clusters as doublets - removing those clusters from the analysis
##Get which clusters are doublets
clust.doublets <- annot.confirmed[which(annot.confirmed$identity.super %in% grep("doublets", annot.confirmed$identity.super, value = T)), "clust"]
clust.doublets <- unlist(clust.doublets$clust)

##Remove all clusters from the long-term states and neighbors that correspond to these doublets
dat.neighbors.lt.clean <- dat.neighbors.lt[!(dat.neighbors.lt$neighbor.clust %in% clust.doublets),]
dat.neighbors.lt.clean <- dat.neighbors.lt[!(dat.neighbors.lt$ltc_clust %in% clust.doublets),]

##In some cases, there are examples where early dividing cells have neighbors in their differentiated counterparts. This is a caveat of this analysis. 
##Identify these clusters and their corresponding neighbors and remove them from analysis (e.g., primitive erythroblasts)
clusts.keep <- list()
clusts.to.remove <- list()
for(cluster in unique(dat.neighbors.lt.clean$ltc_clust)){
  cells.member <- rownames(mama@meta.data)[which(mama@meta.data$clust == cluster)]
  cells.neighbors <- neighbors[[cluster]]
  stage.member <- as.data.frame(table(mama@meta.data[cells.member, "stage.nice"]))
  stage.neighbor <- as.data.frame(table(mama@meta.data[cells.neighbors, "stage.nice"]))
  colnames(stage.member) <- c("stage", "cell_num")
  colnames(stage.neighbor) <- c("stage", "cell_num")
  stage.member$cell_prop <- stage.member$cell_num/sum(stage.member$cell_num)
  stage.neighbor$cell_prop <- stage.neighbor$cell_num/sum(stage.neighbor$cell_num)
  stage.member.good <- stage.member[which(stage.member$cell_prop >= 0.03), ]
  stage.neighbor.good <- stage.neighbor[which(stage.neighbor$cell_prop >= 0.03), ]
  
  ##Now find the farthest stage that cells belong to in the member and neighbor table
  member.stage <- as.numeric(levels(stage.member.good$stage)[as.integer(stage.member.good$stage)])
  neighbor.stage <- as.numeric(levels(stage.neighbor.good$stage)[as.integer(stage.neighbor.good$stage)])
  mem.stg.min <- median(member.stage)
  nbr.stg.min <- median(neighbor.stage)
  
   if(nbr.stg.min - mem.stg.min >= 36){
    clusts.to.remove <- unlist(list(clusts.to.remove, cluster))
  }
  else{
    clusts.keep <- unlist(list(clusts.keep, setdiff(unique(dat.neighbors.lt.clean$ltc_clust), cluster)))
  }
}


##Find clusters to keep as long-term states
clusts.to.keep <- setdiff(unique(dat.neighbors.lt.clean$ltc_clust), unique(dat.neighbors.lt.clean[which(dat.neighbors.lt.clean$ltc_clust %in% clusts.to.remove), "ltc_clust"]))

##Remove clusters that are early and have neighbors in the differentiated counterparts
dat.neighbors.lt.filtered <- dat.neighbors.lt.clean[!(dat.neighbors.lt.clean$ltc_clust %in% clusts.to.remove),]
dat.neighbors.lt.filtered <- dat.neighbors.lt.filtered[!(dat.neighbors.lt.filtered$neighbor.clust %in% clust.doublets),]
#dat.neighbors.lt.filtered <- dat.neighbors.lt.filtered[!(dat.neighbors.lt.filtered$ltc_clust %in% unlist(as.list(clust.doublets))),]

##Ok now add the identity_super columnns from the annotations
##To add annotations from the annotation file, I need to use a similar column name to the annotation spreadsheet. Duplicate the ltc cluster column. 
dat.neighbors.lt.filtered$clust <- dat.neighbors.lt.filtered$ltc_clust
dat.neighbors.lt.filtered <- merge(dat.neighbors.lt.filtered, annot.confirmed[, c("clust", "identity.super", "identity.sub")], by = "clust")

##Extract long-term cycling clusters to be used for downstream analysis
clusts.ltc.filtered <- unique(dat.neighbors.lt.filtered$ltc_clust)

write.table(dat.neighbors.lt.filtered, "~/Box/zfext/global_analysis_curated/transient_vs_long-term/neighbor_clusts_long-term_states_eps_35_thres_36h_filtered_pct_10.tsv")
dat.neighbors.lt.filtered <- read.table("~/Box/zfext/global_analysis_curated/transient_vs_long-term/neighbor_clusts_long-term_states_eps_35_thres_36h_filtered_pct_10.tsv", sep = "")

##Now load the neighbor lists and overlap matrix from the files saved and transferred from the cluster
neighbors <- readRDS("~/Box/zfext/global_analysis_curated/transient_vs_long-term/long-term_clust_with_neighbors_eps_35_all.rds")

##Ok take one set of neighbors - hema.6
#nbr <- neighbors$hema.14

##Subset dat dataframe
#dat.sub <- dat[nbr, ]
#dat.sub <- dat.sub[which(dat.sub$clust == clusts.with.lt), ]

##Find the range of stages within which 80% of the neighbors fall
#low_bound <- quantile(dat.sub$stage.nice, 0.1) ##low bound of range
#high_bound <- quantile(dat.sub$stage.nice, 0.9) ##high bound of range

##Extract the elements that fall within that range
#elements_in_range <- dat.sub[which(dat.sub$stage.nice >= low_bound & dat.sub$stage.nice <= high_bound), ]


##Now iterate that over all tissues
plot.df <- as.data.frame(clusts.ltc.filtered)
for(i in plot.df$clusts.ltc.filtered){
  tissue <- unique(dat.annot[which(dat.annot$clust == i), "tissue.subsets"])
  plot.df[which(plot.df$clust == i), "tissue.subsets"] <- tissue
}

##Do that for each cluster
tissues.with.lt <- c("axial", "endoderm", "epidermis", "eye", "fin", "glial", "hematopoietic", "ionocytes", "mesenchyme", "muscle", "neural", "otic", "periderm", "taste")
for (tissue in tissues.with.lt) {
 
 # Identify clusters from this tissue
  tissue.short <- substr(tissue, 1, 4)
  clusts.do <- grep(tissue.short, clusts.ltc.filtered, value = T)
  
  ##Get neighbors for these clusters
  for (c in clusts.do){
    neighbr <- neighbors[[c]]
    dat.sub <- dat[neighbr, ]
    dat.sub <- dat.sub[which(dat.sub$clust == clusts.with.lt), ]
    
  ##Find the range of stages within which 80% of the neighbors fall
    low_bound <- quantile(dat.sub$stage.nice, 0.1) ##low bound of range
    high_bound <- quantile(dat.sub$stage.nice, 0.9) ##high bound of range
    
    ##Extract the elements that fall within that range
    elements_in_range <- dat.sub[which(dat.sub$stage.nice >= low_bound & dat.sub$stage.nice <= high_bound), ]
    
    plot.df[which(plot.df$clusts.ltc.filtered == c), "stage.start"] <- low_bound
    plot.df[which(plot.df$clusts.ltc.filtered == c), "stage.stop"] <- high_bound
  }
}


##Now get rid of the clusts where the differences are less than 36h
plot.df$diff_stg <- plot.df$stage.stop - plot.df$stage.start
plot.df.good <- plot.df[-which(plot.df$diff_stg <= 36), ]
plot.df.good <- plot.df.good[!is.na(plot.df.good$stage.start), ]
plot.df.ltc <- plot.df.good[which(plot.df.good$clusts.ltc.filtered %in% clusts.to.keep), ]

##Add the annotations
colnames(plot.df.ltc) <- c("clust", "tissue.subsets", "stage.start", "stage.stop", "diff_stg")
# Load cluster annotations 
is.empty <- function(x) {
  is.na(x) | is.null(x) | (sapply(trimws(x), nchar) == 0)
}
annot.confirmed <- readxl::read_xlsx("~/Box/zfext/02-Clustering/mama_tissue_clustering/annot/mama_annotations_full_extended_jeffedit_newtissues_curated_clustering.xlsx")
annot.confirmed <- annot.confirmed[!is.empty(annot.confirmed$clust),3:ncol(annot.confirmed)]
rownames(annot.confirmed) <- annot.confirmed$clust
annot.confirmed$confirmed <- T
plot.df.annot <- merge(plot.df.ltc, annot.confirmed[, c("clust", "identity.super", "identity.sub")], by = "clust")

##Plot all long-term states
pdf(file = paste0("cell_states_time_lengths_cyc_new.pdf"), width = 28, height = 36)
ggplot(data=plot.df.annot)+
  geom_segment(aes(x=stage.start, xend=stage.stop, y = clust, yend = clust),lwd=4)
dev.off()


##Subset into individual tissues and plot
tissue <- "eye"
df <- plot.df.annot[which(plot.df.annot$tissue.subsets == tissue), ]

pdf(file = paste0("cell_states_time_", tissue, "_lengths_cyc.pdf"), width = 28, height = 36)
ggplot(data=df)+
  geom_segment(aes(x=stage.start, xend=stage.stop, y = clust, yend = clust),lwd=4)+
  facet_grid(tissue.subsets~.)+xlab("HPF")

tissue <- "neural"
df <- plot.df.annot[which(plot.df.annot$tissue.subsets == tissue), ]

pdf(file = paste0("cell_states_time_", tissue, "_lengths_cyc.pdf"), width = 28, height = 36)
ggplot(data=df)+
  geom_segment(aes(x=stage.start, xend=stage.stop, y = clust, yend = clust),lwd=4)+
  facet_grid(tissue.subsets~.)+xlab("HPF")
dev.off()

tissue <- "epidermis"
df <- plot.df.annot[which(plot.df.annot$tissue.subsets == tissue), ]

pdf(file = paste0("cell_states_time_", tissue, "_lengths_cyc.pdf"), width = 28, height = 36)
ggplot(data=df)+
  geom_segment(aes(x=stage.start, xend=stage.stop, y = clust, yend = clust),lwd=4)+
  facet_grid(tissue.subsets~.)+xlab("HPF")
dev.off()

tissue <- "hematopoietic"
df <- plot.df.annot[which(plot.df.annot$tissue.subsets == tissue), ]
ggplot(data=df)+
  geom_segment(aes(x=stage.start, xend=stage.stop, y = clust, yend = clust),lwd=4)+
  facet_grid(tissue.subsets~.)+xlab("HPF")

tissue <- "muscle"
df <- plot.df.annot[which(plot.df.annot$tissue.subsets == tissue), ]
ggplot(data=df)+
  geom_segment(aes(x=stage.start, xend=stage.stop, y = clust, yend = clust),lwd=4)+
  facet_grid(tissue.subsets~.)+xlab("HPF")

tissue <- "glial"
df <- plot.df.annot[which(plot.df.annot$tissue.subsets == tissue), ]
ggplot(data=df)+
  geom_segment(aes(x=stage.start, xend=stage.stop, y = clust, yend = clust),lwd=4)+
  facet_grid(tissue.subsets~.)+xlab("HPF")

tissue <- "mesenchyme"
df <- plot.df.annot[which(plot.df.annot$tissue.subsets == tissue), ]
ggplot(data=df)+
  geom_segment(aes(x=stage.start, xend=stage.stop, y = clust, yend = clust),lwd=4)+
  facet_grid(tissue.subsets~.)+xlab("HPF")

tissue <- "periderm"
df <- plot.df.annot[which(plot.df.annot$tissue.subsets == tissue), ]
ggplot(data=df)+
  geom_segment(aes(x=stage.start, xend=stage.stop, y = clust, yend = clust),lwd=4)+
  facet_grid(tissue.subsets~.)+xlab("HPF")

tissue <- "taste"
df <- plot.df.annot[which(plot.df.annot$tissue.subsets == tissue), ]
ggplot(data=df)+
  geom_segment(aes(x=stage.start, xend=stage.stop, y = clust, yend = clust),lwd=4)+
  facet_grid(tissue.subsets~.)+xlab("HPF")

tissue <- "ionocytes"
df <- plot.df.annot[which(plot.df.annot$tissue.subsets == tissue), ]
ggplot(data=df)+
  geom_segment(aes(x=stage.start, xend=stage.stop, y = clust, yend = clust),lwd=4)+
  facet_grid(tissue.subsets~.)+xlab("HPF")

tissue <- "otic"
df <- plot.df.annot[which(plot.df.annot$tissue.subsets == tissue), ]
ggplot(data=df)+
  geom_segment(aes(x=stage.start, xend=stage.stop, y = clust, yend = clust),lwd=4)+
  facet_grid(tissue.subsets~.)+xlab("HPF")


##Check which of these are neighbors of each other 
tissue <- "epidermis"
tissue.short <- substr(tissue, 1, 4)
clusts.do <- grep(tissue.short, clusts.ltc.filtered, value = T)

epid.neighbors <- vector("list", length(clusts.do))
# Identify the neighbors of each cluster
for (c in clusts.do) {
  epid.neighbors[[c]] <- neighbors[[c]]
}

clusts.epid <- vector("list", length(clusts.epid))
for(c in clusts.do){
  nb <- unique(epid.neighbors[[c]])
  clusts.epid[[c]] <- unique(dat[nb, "clust"])
}


##Calculate some statistics
cells.under.24 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff <= 24)]
cells.under.12 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff <= 12)]
cells.above.48 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff <= 48)]
cells.below.36 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff <= 36)]
cells.below.21 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff <= 21 & mama@meta.data$cc.phase == "non-cycling")]


cells.72.120 <- rownames(mama@meta.data)[which(mama@meta.data$stage.nice > 72)]
cells.non.cycling.72 <- cells.72.120[cells.72.120 %in% cells.non.cycling]
table(mama@meta.data[cells.72.120, "states.above.36"])

short.non_cycling.cells.72 <- short.non_cycling.cells[which(short.non_cycling.cells %in% cells.72.120)]
table(mama@meta.data[short.non_cycling.cells.72, "tissue.subsets"])


cells.10.12 <- rownames(mama@meta.data)[which(mama@meta.data$stage.nice < 12)]
cells.21 <- rownames(mama@meta.data)[which(mama@meta.data$stage.nice <= 21)]
cells.endo <- rownames(mama@meta.data)[which(mama@meta.data$tissue.subsets == "endoderm")]
cells.to.check <- cells.10.12[cells.10.12 %in% cells.non.cycling]
cells.to.check <- cells.21[cells.21 %in% cells.non.cycling]

table(mama@meta.data[cells.to.check, "states.above.36"])

table(mama@meta.data[cells.stnc, ])

##Plot area plots for the four different states for the different tissues
##Ok first make a dataframe such that for each tissues, make a column for stage, column for the four categories, n is total cells, and percentage is the proportion
mama@meta.data$tissue_stage.diff <- paste0(mama@meta.data$tissue.subsets, "@", mama@meta.data$stage.diff.cc)
mama@meta.data$tissue_stage.diff.sg <- paste0(mama@meta.data$tissue_stage.diff, ".", mama@meta.data$stage.group)
df <- as.data.frame(table(mama@meta.data$tissue_stage.diff.sg))
colnames(df) <- c("tissue_stage.diff_sg", "cell_number")
df$tissue <- unlist(lapply(strsplit(x = as.character(df$tissue_stage.diff_sg), split = "@"), function(x) x[1]))
df$stage <- unlist(lapply(strsplit(x = as.character(df$tissue_stage.diff_sg), split = "\\."), function(x) x[2]))
df$stage.diff.sg <- unlist(lapply(strsplit(x = as.character(df$tissue_stage.diff_sg), split = "@"), function(x) x[2]))
df$stage.diff <- unlist(lapply(strsplit(x = as.character(df$stage.diff.sg), split = "\\."), function(x) x[1]))
df$stage.nice <- unlist(lapply(strsplit(x = as.character(df$stage.diff.sg), split = "\\."), function(x) x[2]))
df$stage.nice <- str_split_fixed(df$stage.nice, "-", 2)
df.good <- df[, c("tissue", "stage.nice", "stage.diff", "cell_number")]
df.good$stage.nice <- as.numeric(df.good$stage.nice[,2])
colnames(df.good) <- c("tissue", "stage", "stage.diff", "n")
#df.good$stage <- as.numeric(df.good$stage[,2])
#df.good$stage <- df.good$stage[,1]
df.good$stage <- as.numeric(df.good$stage)

saveRDS(df.good, "~/Box/zfext/global_analysis_curated/transient_vs_long-term/tissue_stage_dataframe_area_plots_with_neural.rds")


##Plot area plots for all tissues
library(dplyr)
for(tissue in unique(df.good$tissue)){
  df.tissue <- df.good[which(df.good$tissue == tissue), ]
  data <- df.tissue  %>%
    group_by(stage, stage.diff) %>%
    summarise(n = sum(n)) %>%
    mutate(percentage = n / sum(n))
  
  data$stage <- as.numeric(data$stage)
  
  p <- ggplot(data, aes(x=stage, y=percentage, fill=stage.diff)) + 
    geom_area(alpha=0.6 , size=1, colour="black")
  
  ggsave(filename = paste0(tissue, "_area_plot.pdf"), plot = p, 
         device = "pdf", path = "~/Box/zfext/global_analysis_curated/transient_vs_long-term/figures/")
  
}

df.good <- readRDS("~/Box/zfext/global_analysis_curated/transient_vs_long-term/tissue_stage_dataframe_area_plots.rds")

tissues <- unique(df.good$tissue)
stages <- sort(unique(df.good$stage))
category <- unique(df.good$stage.diff)[1:8]

df.good.norm <- do.call("rbind", lapply(tissues, function(tissue) {
  df.good.sub <- df.good[df.good$tissue == tissue,]
  mat.good.sub <- reshape2::acast(data = df.good.sub, formula = "stage ~ stage.diff", fill = 0)
  mat.good.sub.norm <- sweep(mat.good.sub, 1, rowSums(mat.good.sub), "/")
  melt.good.sub.norm <- reshape2::melt(mat.good.sub.norm)
  colnames(melt.good.sub.norm) <- c("stage", "stage.diff", "proportion")
  melt.good.sub.norm$tissue <- tissue
  return(melt.good.sub.norm)
}))

saveRDS(df.good.norm, file = "~/Box/zfext/global_analysis_curated/transient_vs_long-term/tissue_stage_dataframe_area_plots_with_neural.rds")

area.colors <- setNames(c("#DAA520", "#B8860B", "#FFA500", "#F08080", "#87CEEB", "#7FFFD4", "#00CED1", "#008080"), c("<24hr_cycling", "24-36hr_cycling", "36-48hr_cycling", ">48hr_cycling", "<24hr_non-cycling", "24-36hr_non-cycling", "36-48hr_non-cycling", ">48hr_non-cycling"))
df.good.norm$cc.status <- str_split_fixed(df.good.norm$stage.diff, pattern = "_", 2)
df.good.norm$cc.status <- as.character(df.good.norm$cc.status[,2])
colnames(df.good.norm) <- c("stage", "stage.diff", "proportion", "tissue", "cc.phase")


df.good.norm.order <- df.good.norm[order(df.good.norm$stage.diff, decreasing = FALSE),]


png(file = "~/Box/zfext/global_analysis_curated/transient_vs_long-term/figures/tissue_stage_area_plots_v3.png")
ggplot(data = df.good.norm, aes(x = stage, y = proportion, fill = interaction(stage.diff, cc.phase))) + geom_area() + facet_wrap(vars(tissue), ncol = 4, strip.position = "right") + theme_bw() + scale_fill_manual(values = c("#DAA520", "#B8860B", "#FFA500", "#FF4500", "#87CEEB", "#7FFFD4", "#00CED1", "#008080"))
dev.off()


##Save the stage.diff column as a dataframe with cell_ids (eps = 35) and export  as a dataframe
##Subset dataframe
dm <- mama@meta.data[, c("stage.nice", "stage.group", "tissue.subsets", "clust", "stage.diff", "long_vs_short", "cc.status", "states.above.36", "stage.diff.cc")]
write.table(dm, "~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_perdurance_stage_diff_eps_35_df.tsv", quote=FALSE, sep='\t')


dat.neighbors.lt.filtered[dat.neighbors.lt.filtered$ltc_clust == "pigm.7", ]
