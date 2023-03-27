##This code is designed to find long-term persistent states and short term states based on arbitrary thresholds
#Load libraries
library(Seurat)
library(dbscan)

save.path <- "~/Box/zfext/global_analysis_curated/transient_vs_long-term/knn_neighbors_data_thres_36h/"
plot.path <- "~/Box/zfext/global_analysis_curated/transient_vs_long-term/plots_neighbors_data_thres_36h/"

##Load mama
mama <- readRDS("~/Box/zfext/02-Clustering/mama+ds_integrated/merged_mama_ds_slim_seurat_25Jan2022.rds")

#Load the reductions slot for mama
mama.reductions <- readRDS("~/Desktop/merged_mama_ds_reductions_2022_01-18.rds")
mama@reductions <- mama.reductions

tissue.subsets <- readRDS("~/Box/zfext/02-Clustering/mama+ds_integrated/tissue_annot_2.rds")
mama@meta.data$tissue.annot <- tissue.subsets

# Load cell annotations to introduce for downsample purposes
cell.annot.path <- "~/Box/zfext/annotations_celltype_curated_newMama/celltype_annotations_df.tsv"
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

##Load previously calculated categories of stage.differences and long-term/short-term states
stage.diff <- readRDS("~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_stage.diff_eps_35.rds")
long_vs_short <- readRDS("~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_long_vs_short_stage_diff_5NN.rds")
mama@meta.data$long_vs_short <- long_vs_short
mama@meta.data$stage.diff <- stage.diff


######### Cell states above and below certain thresholds #############
##Now use a threshold of 36 hours - cells that have a stage difference more than 36 hours will be categorized as long-term states and ones which have stage difference less than 36 hours will be defined short term
##Define a slot in mama for long-term and short-term states
mama@meta.data$states.above.36 <- NA
mama@meta.data[which(mama@meta.data$stage.diff >= 36 & mama@meta.data$cc.status == "cycling"), "states.above.36"] <- "long-term_cycling"
mama@meta.data[which(mama@meta.data$stage.diff >= 36 & mama@meta.data$cc.status == "non-cycling"), "states.above.36"] <- "long-term_non-cycling"
mama@meta.data[which(mama@meta.data$stage.diff <= 36 & mama@meta.data$cc.status == "cycling"), "states.above.36"] <- "short-term_cycling"
mama@meta.data[which(mama@meta.data$stage.diff <= 36 & mama@meta.data$cc.status == "non-cycling"), "states.above.36"] <- "short-term_non-cycling"

colors.threshold <- setNames(c("#E9967A", "#800000", "#08519c", "#40E0D0"), c("long-term_cycling", "long-term_non-cycling", "short-term_cycling", "short-term_non-cycling"))
dpi <- 300
png("~/Box/zfext/global_analysis_curated/transient_vs_long-term/figures/mama_cell_states_thresh_36hr_eps_35.png", height = 5*dpi, width = 5*dpi)
DimPlot(mama, group.by = "states.above.36", cols = colors.threshold, raster = F) + NoAxes() + NoLegend()
dev.off()

##Save the persistent and short-term states assignment for 36 hours as a threshold
saveRDS(mama@meta.data$states.above.36, "~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_cell_states_threshold_36_eps_35.rds")

##Plot a scatterpie plot for these cells that have stage difference above 36 hours
##Load mama metadata slot
states.above.36 <- readRDS("~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_cell_states_threshold_36_eps_35.rds")
mama@meta.data$states.ab <- states.above.36
DimPlot(mama, group.by = "states.above.36", raster = F) + NoAxes()

long.cycling.cells <- rownames(mama@meta.data)[mama@meta.data$states.above.36 == "long-term_cycling"]
long.non_cycling.cells <- rownames(mama@meta.data)[mama@meta.data$states.above.36 == "long-term_non-cycling"]
short.cycling.cells <- rownames(mama@meta.data)[mama@meta.data$states.above.36 == "short-term_cycling"]
short.non_cycling.cells <- rownames(mama@meta.data)[mama@meta.data$states.above.36 == "short-term_non-cycling"]

df <- as.data.frame(table(mama@meta.data$clust))
colnames(df) <- c("clust", "total.cells")
lt_diff <- as.data.frame(table(mama@meta.data[long.non_cycling.cells, "clust"]))
lt_cyc <- as.data.frame(table(mama@meta.data[long.cycling.cells, "clust"]))
st_diff <- as.data.frame(table(mama@meta.data[short.non_cycling.cells, "clust"]))
st_cyc <- as.data.frame(table(mama@meta.data[short.cycling.cells, "clust"]))

colnames(lt_diff) <- c("clust", "lt_diff")
colnames(lt_cyc) <- c("clust", "lt_cyc")
colnames(st_diff) <- c("clust", "st_diff")
colnames(st_cyc) <- c("clust", "st_cyc")

df.full <- merge(df, lt_diff, all = T)
df.full <- merge(df.full, lt_cyc, all = T)
df.full <- merge(df.full, st_diff, all = T)
df.full <- merge(df.full, st_cyc, all = T)
colnames(df.full) <- c("clust", "total.cells", "lt_diff", "lt_cyc", "st_diff", "st_cyc")

##So we now have a dataframe that has 3 columns with the tissue_sg names, total cells and the number of long-term cells
##Convert all NA values to 0
df.full[is.na(df.full)] <- 0
df.full$lt_diff_prop <- as.numeric(df.full$lt_diff/df.full$total.cells)
df.full$lt_diff_prop <- sprintf(df.full$lt_diff_prop, fmt = '%#.3f')        # Apply sprintf function

df.full$lt_cyc_prop <- as.numeric(df.full$lt_cyc/df.full$total.cells)
df.full$lt_cyc_prop <- sprintf(df.full$lt_cyc_prop, fmt = '%#.3f')

df.full$st_cyc_prop <- as.numeric(df.full$st_cyc/df.full$total.cells)
df.full$st_cyc_prop <- sprintf(df.full$st_cyc_prop, fmt = '%#.3f')

df.full$st_diff_prop <- as.numeric(df.full$st_diff/df.full$total.cells)
df.full$st_diff_prop <- sprintf(df.full$st_diff_prop, fmt = '%#.3f')

##Which clusters have long-term states?
clusts.ltc <- unlist(df.full[which(df.full$lt_cyc_prop >= 0.15), "clust"])
clusts.ltnc <- df.full[which(df.full$lt_diff_prop >= 0.15), "clust"]

##Split the first column into 2 columns
#library(stringr)
#df.full[c("tissue", "num")] <- str_split_fixed(df.full$clust, ".", 1)
#df.full$tissue <- unlist(lapply(strsplit(x = as.character(df.full$clust), split = "\\."), function(x) x[1]))
#df.full$tissue.subsets <- cbind()
##Rearrange columns and remove original column names
#df.full <- df.full[c("tissue", "stage.group", "total.cells", "lt_diff", "lt_cyc", "st_diff", "st_cyc")]
#df.prop <- df.full[c("tissue", "stage.group", "total.cells", "lt_diff_prop", "lt_cyc_prop", "st_cyc_prop", "st_diff_prop")]

##Ok, first step, we have to get the cells that are long-term within these clusters. 
##Get the actual cell ids
cc.score <- readRDS("~/Box/zfext/02-Clustering/mama+ds_integrated/merged_mama_cellCycle_score.rds")
mama@meta.data <- cbind(mama@meta.data, cc.score)
dat <- mama@meta.data[,c("stage.nice", "stage.group", "tissue.stage", "tissue.subsets", "clust", "clust.sg", "cc.phase", "long_vs_short",  "cc.status", "stage.diff.cc", "states.above.36", "stage.diff")]
p <- DimPlot(mama, group.by = "stage.diff.cc", raster = F) + NoAxes()
color.stg.cc <- unique(ggplot_build(p)$data[[1]]$colour)

# Load cluster annotations 
is.empty <- function(x) {
  is.na(x) | is.null(x) | (sapply(trimws(x), nchar) == 0)
}
annot.confirmed <- readxl::read_xlsx("~/Box/zfext/02-Clustering/mama_tissue_clustering/annot/mama_annotations_full_extended_jeffedit_newtissues_curated_clustering.xlsx")
annot.confirmed <- annot.confirmed[!is.empty(annot.confirmed$clust),3:ncol(annot.confirmed)]
rownames(annot.confirmed) <- annot.confirmed$clust
annot.confirmed$confirmed <- T

##Merge the annotations to the dataframe
colnames(annot.confirmed$UI) <- "clust"
dat.annot <- merge(dat, annot.confirmed[, c("clust", "identity.super")], by = "clust")
dat <-  dat.annot
dat$tissue <- unlist(lapply(strsplit(x = dat$tissue.stage, split = "_"), function(x) x[1]))



##Get the long-term states at a time-threshold of 36 hours
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
thresh.states <- 0.15
clusts.ltc <- rownames(cluster.proportions)[cluster.proportions$long_cycling >= thresh.states]
clusts.ltnc <- rownames(cluster.proportions)[cluster.proportions$`long_non-cycling` >= thresh.states]
clusts.stnc <- rownames(cluster.proportions)[cluster.proportions$`short_non-cycling` >= 0.85]

# Only operate tissues that have a long-term state
clusts.with.lt <- unlist(unique(list(clusts.ltc, clusts.ltnc)))
tissues.with.lt <- sort(unique(URD::cutstring(c(clusts.ltc, clusts.ltnc), delimiter = "\\.", field = 1)))

# Identify members of each cluster
members <- lapply(clusts.ltc, function(c) rownames(cell.annot)[cell.annot$clust == c])
names(members) <- clusts.ltc

##Identify neighbors for these clusters that are long-term cycling
# Initialize empty list to hold neighbors results across all tissues
neighbors <- vector("list", length(clusts.ltc))
names(neighbors) <- clusts.ltc

# Loop per tissue because neighbor lists are per-tissue
tissues.with.lt <- c("endoderm", "epidermis", "eye", "glia", "hematopoietic", "ionocytes", "mesenchyme", "muscle", "neural", "otic", "periderm", "taste")
for (tissue in tissues.with.lt) {
  # Load nn file
  nn <- readRDS(paste0(save.path, tissue, "_neighbor_matrix_35_binary.rds")) #Fill this part in.
  
  # Identify clusters from this tissue
  tissue.short <- substr(tissue, 1, 4)
  clusts.do <- grep(tissue.short, clusts.ltc, value = T)
  
  # Identify the neighbors of each cluster
  for (c in clusts.do) {
    neighbors[[c]] <- names(rowSums(nn[members[c], ] == 0))
  }
}


# Create an overlap matrix
overlap.ltc <- as.data.frame(combn(clusts.ltc, 2)) # may need t()
overlap.ltc <- t(overlap.ltc)
colnames(overlap.ltc) <- c("clust1", "clust2")
overlap.ltc$overlap <- apply(overlap.ltc, 1, function(x) {
  pair.memb <- unlist(members[x[1:2]])
  pair.neighb <- unique(unlist(neighbors[x[1:2]]))
  memb.in.neighb <- intersect(pair.memb, pair.neighb)
  prop <- length(pair.memb) / length(memb.in.neighb)
  return(prop)
})


##Ok get neighbors for one cluster - hema.6
neighbors.6 <- names(rowSums(nn[members$hema.6, ] == 0))
neighbors.14 <- names(rowSums(nn[members$hema.14, ] == 1)) 


dat <- mama@meta.data[,c("stage.nice", "stage.group", "cc.status", "g2m.score", "s.score")]
write.table(dat, "~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_cell_cycle_scores.tsv")



