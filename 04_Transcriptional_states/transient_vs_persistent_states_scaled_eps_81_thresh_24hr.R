##This code is to use the distances calculated on the scaled gene expression space and deduce short-term and long-term 
##states at the threshold of 30h. 

##Load libraries
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(gridExtra)

##Specify paths
base.path <- "/data/CSD/zfext/LTA/results/13-EPS/"
neighbor.path.data <- "/data/CSD/zfext/LTA/results/13-EPS/knn_neighbors_data/"
neighbor.path.scaled <- "/data/CSD/zfext/LTA/results/13-EPS/knn_neighbors_scaled/"
#cells.path <- "/data/CSD/zfext/LTA/results/13-EPS/long-term_cells_eps_40/"
#save.path.data <- "/data/CSD/zfext/LTA/results/13-EPS/knn_neighbors_data/"
save.path.scaled.small <- "/data/CSD/zfext/LTA/results/13-EPS/states_scaled_eps_81/"
dist.path.scaled <- "/data/CSD/zfext/results/04f-SubsetsV6/dist_subsets_scaled/"
dist.path.data <- "/data/CSD/zfext/results/04f-SubsetsV6/dist_subsets_data/"


##Load mama
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


##Load the various stage.difference objects calculated on the cluster and plot on a UMAP how long individual cells persist for
##Ok, now generate a script such that we can sort every single cell based on their stage difference
##Load the "nearest" objects for each tissue and save the values as a slot in mama
sample <- "axial"
eps <- "81"
nearest.axial <- readRDS(file = paste0(neighbor.path.scaled, sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.axial <- as.data.frame(nearest.axial)
colnames(nearest.axial) <- "stage.diff"

sample <- "epidermis"
nearest.epidermis <- readRDS(file = paste0(neighbor.path.scaled, sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.epidermis <- as.data.frame(nearest.epidermis)
colnames(nearest.epidermis) <- "stage.diff"

sample <- "blastomeres"
nearest.blastomeres <- readRDS(file = paste0(neighbor.path.scaled, sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.blastomeres <- as.data.frame(nearest.blastomeres)
colnames(nearest.blastomeres) <- "stage.diff"

sample <- "endoderm"
nearest.endo <- readRDS(file = paste0(neighbor.path.scaled, sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.endo <- as.data.frame(nearest.endo)
colnames(nearest.endo) <- "stage.diff"

sample <- "eye"
nearest.eye <- readRDS(file = paste0(neighbor.path.scaled, sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.eye <- as.data.frame(nearest.eye)
colnames(nearest.eye) <- "stage.diff"

sample <- "fin"
nearest.fin <- readRDS(file = paste0(neighbor.path.scaled, sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.fin <- as.data.frame(nearest.fin)
colnames(nearest.fin) <- "stage.diff"

sample <- "glial"
nearest.glia <- readRDS(file = paste0(neighbor.path.scaled, sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.glia <- as.data.frame(nearest.glia)
colnames(nearest.glia) <- "stage.diff"

sample <- "hematopoietic"
nearest.hema <- readRDS(file = paste0(neighbor.path.scaled, sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.hema <- as.data.frame(nearest.hema)
colnames(nearest.hema) <- "stage.diff"

sample <- "ionocytes"
nearest.iono <- readRDS(file = paste0(neighbor.path.scaled, sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.iono <- as.data.frame(nearest.iono)
colnames(nearest.iono) <- "stage.diff"

sample <- "mesenchyme"
nearest.mese <- readRDS(file = paste0(neighbor.path.scaled, sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.mese <- as.data.frame(nearest.mese)
colnames(nearest.mese) <- "stage.diff"

sample <- "mural"
nearest.mural <- readRDS(file = paste0(neighbor.path.scaled, sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.mural <- as.data.frame(nearest.mural)
colnames(nearest.mural) <- "stage.diff"

sample <- "muscle"
nearest.muscle <- readRDS(file = paste0(neighbor.path.scaled, sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.muscle <- as.data.frame(nearest.muscle)
colnames(nearest.muscle) <- "stage.diff"

sample <- "neural"
nearest.neural <- readRDS(file = paste0(neighbor.path.scaled, sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.neural <- as.data.frame(nearest.neural)
colnames(nearest.neural) <- "stage.diff"

sample <- "otic"
nearest.otic <- readRDS(file = paste0(neighbor.path.scaled, sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.otic <- as.data.frame(nearest.otic)
colnames(nearest.otic) <- "stage.diff"

sample <- "periderm"
nearest.periderm <- readRDS(file = paste0(neighbor.path.scaled, sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.periderm <- as.data.frame(nearest.periderm)
colnames(nearest.periderm) <- "stage.diff"

sample <- "pgc"
nearest.pgc <- readRDS(file = paste0(neighbor.path.scaled, sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.pgc <- as.data.frame(nearest.pgc)
colnames(nearest.pgc) <- "stage.diff"

sample <- "pigment"
nearest.pigment <- readRDS(file = paste0(neighbor.path.scaled, sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.pigment <- as.data.frame(nearest.pigment)
colnames(nearest.pigment) <- "stage.diff"

sample <- "pronephros"
nearest.pron <- readRDS(file = paste0(neighbor.path.scaled, sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.pron <- as.data.frame(nearest.pron)
colnames(nearest.pron) <- "stage.diff"

sample <- "taste"
nearest.taste <- readRDS(file = paste0(neighbor.path.scaled, sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.taste <- as.data.frame(nearest.taste)
colnames(nearest.taste) <- "stage.diff"

sample <- "cephalic"
nearest.cephalic <- readRDS(file = paste0(neighbor.path.scaled, sample, "_mean_stage_diff_eps_", eps, ".rds"))
nearest.cephalic <- as.data.frame(nearest.cephalic)
colnames(nearest.cephalic) <- "stage.diff"

##Combine all tissue specific dataframes
nearest.total <- rbind(nearest.axial, nearest.epidermis, nearest.blastomeres, nearest.endo, nearest.eye, nearest.fin, nearest.glia, nearest.hema,
                       nearest.iono, nearest.mese, nearest.mural, nearest.muscle, nearest.neural, nearest.otic, nearest.periderm, nearest.pgc, 
                       nearest.pigment, nearest.pron, nearest.taste, nearest.cephalic)

#nearest.total <- nearest.total[WhichCells(mama), ]


mama@meta.data <- cbind(mama@meta.data, nearest.total[rownames(mama@meta.data), ])
colnames(mama@meta.data)[67] <- c("stage.diff")

##Get rid of all NaN values
mama@meta.data[which(mama@meta.data$stage.diff == "NaN"), "stage.diff"] <- 0

##Save the stage difference column 
saveRDS(mama@meta.data$stage.diff, paste0(save.path.scaled.small, "mama_stage.diff_data_eps_81.rds"))
stage.diff <- readRDS(paste0(save.path.scaled.small, "mama_stage.diff_data_eps_81.rds"))
mama@meta.data$stage.diff <- stage.diff

DimPlot(mama, group.by = "stage.diff", raster = F) + NoAxes() + NoLegend() + scale_fill_continuous(type = "gradient")
library(RColorBrewer)
dpi <- 300
pdf(paste0(save.path.scaled.small, "figures/mama_stage.diff_eps_81.pdf"), width = 5*dpi, height = 5*dpi)
DimPlot(mama, group.by = "stage.diff", raster = F) + NoAxes() + NoLegend() + scale_fill_continuous(type = "gradient")
dev.off()

##Categorize the different cells into how long they persist for using this eps setting (eps = 81)
##Define cells that persist for different lengths of time. 
cells.under.24 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff < 24)]
cells.24.36 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff >= 24 & mama@meta.data$stage.diff < 36)]
cells.36.48 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff >= 36 & mama@meta.data$stage.diff < 48)]
cells.48 <- rownames(mama@meta.data)[which(mama@meta.data$stage.diff >= 48)]

##Now assign a slot in mama and plot these cells
mama@meta.data$long_vs_short <- NA
mama@meta.data[cells.under.24, "long_vs_short"] <- "<24hr"
mama@meta.data[cells.24.36, "long_vs_short"] <- "24-36hr"
mama@meta.data[cells.36.48, "long_vs_short"] <- "36-48hr"
mama@meta.data[cells.48, "long_vs_short"] <- ">48hr"

##Generate the UMAP plot showing cells belonging to these different categories 
colors.stg <- setNames(c("#D35FB7", "#40B0A6", "#FFC20A", "#005AB5"), c("<24hr", "24-36hr", "36-48hr",  ">48hr"))
dpi <- 300
png(paste0(save.path.scaled.small, "/figures/mama_transient_vs_persistent_states_scaled_eps_81_v2.png"), height = 5*dpi, width = 5*dpi)
DimPlot(mama, group.by = "long_vs_short", cols = colors.stg, raster = F, na.value = "#f5f7f6") + NoAxes() + NoLegend()
dev.off()

saveRDS(mama@meta.data$long_vs_short, paste0(save.path.scaled.small, "mama_long_vs_short_stage_diff_scaled_eps_81.rds"))
long_vs_short <- readRDS(paste0(save.path.scaled.small, "mama_long_vs_short_stage_diff_scaled_eps_81.rds"))
mama@meta.data$long_vs_short <- long_vs_short

##Calculate cell cycle scores using Seurat

##Define G1/S phase genes
s.phase.genes <- c("mcm5", "pcna", "tyms", "mcm7", "mcm4", "rrm1", "ung1", "gins2", "mcm6", "cdca7", "dtl", "prim1", "uhrf1", "cenpu", "gmnn", "hells", 
                   "ccne2", "cdc6", "rfc2", "polr1b", "nasp", "rad51ap1", "wdr76", "slbp", "ubr7", "pold3", "msh2", "atad2", "rad51", "rrm2", "cdc45", 
                   "exo1", "tipin", "dscc1", "blm", "casbap2", "usp1", "clspn", "pola1", "chaf1b", "mrpl36", "e2f8")

##Define G2/M phase genes
g2m.phase.genes <- c("cdk1", "ube2c", "birc5", "top2a", "tpx2", "cks2", "nuf2", "mki67", "tacc3", "cenpf", "smc4", "ckap4", "kif11", "cdca3", "hmgb2", 
                     "ndc81", "cks1b", "tmpo", "pimreg", "ccnb2", "ckap2l", "ckap2", "aurkb", "bub1", "anp32e", "tubb4b", "gtse1", "kif20b", "hjurp", 
                     "jpt1", "cdc20", "ttk", "cdc25c", "kif2c", "rangap1", "ncapd2", "dlgap5", "cdca8", "cdca2", "ect2", "kif23", "hmmr", "aurka", "psrc1",
                     "anln", "lbr", "ckap5", "cenpe", "ctcf", "nek2", "g2e3", "gas2l3", "cbx5", "cenpa")

##Using the genes above, calculate the cell cycle scores
cc.score <- CellCycleScoring(mama, s.features = s.phase.genes, g2m.features = g2m.phase.genes)                                 

##Save the cc.scores
saveRDS(cc.score, "~/Box/zfext/02-Clustering/mama+ds_integrated/merged_mama_cellCycle_score.rds")          

##Now divide into groups that includes the above categories as well as their cycling status (i.e., cycling and non-cycling)
cc.score <- readRDS(paste0(save.path.scaled.small, "merged_mama_cellCycle_score.rds"))
mama@meta.data <- cbind(mama@meta.data, cc.score)

##Find which cells belong to G1/S phase and which belong to G2/M phase
cells.s.phase <- rownames(mama@meta.data)[mama@meta.data$s.score > 0]
cells.g2m.phase <- rownames(mama@meta.data)[mama@meta.data$g2m.score > 0]

##Define a slot in the global metadata assigning cells to the different cell.cycle phases
mama@meta.data$cc.phase <- NA
mama@meta.data[cells.s.phase, "cc.phase"] <- "G1/S-phase"
mama@meta.data[cells.g2m.phase, "cc.phase"] <- "G2M-phase"


##Use the above information to classify cells into cycling and non-cycling groups
##Classify cells into either cycling or non-cycling based on the above cell cycle scores
cells.cycling <- rownames(mama@meta.data)[which(mama@meta.data$cc.phase == "G1/S-phase")]
cells.also.cycling <- rownames(mama@meta.data)[which(mama@meta.data$cc.phase == "G2M-phase")]

##Total cycling  cells
cells.cycling.total <- unlist(unique(list(cells.cycling, cells.also.cycling)))

##Non cycling  cells
#cells.non.cycling <- rownames(mama@meta.data)[which(is.na(mama@meta.data$cc.phase))]
cells.non.cycling <- setdiff(WhichCells(mama), cells.cycling.total)

##Add these cycling and non-cycling groups to the global dataset metadata
mama@meta.data$cc.status <- NA
mama@meta.data[cells.non.cycling, "cc.status"] <- "non-cycling"
mama@meta.data[cells.cycling.total, "cc.status"] <- "cycling"

##Add non-cycling cells to the "cc.phase" column
mama@meta.data[cells.non.cycling, "cc.phase"] <- "non-cycling"

color.stg.cc <- setNames(c("#4728E6", "#F2D9D0", "#3B5A57"), c("G1/S-phase", "G2M-phase", "non-cycling"))

##Plot the UMAP showing cells in the different cell cycle phases - Figure S2A, B
png(paste0(save.path.scaled.small, "/figures/mama_cell_cycle_score_v2.png"), width = 5*dpi, height = 5*dpi)
DimPlot(mama, group.by = "cc.phase", raster = F, cols = color.stg.cc) + NoAxes() + NoLegend()
dev.off()

##Plot the cells categorized into their length of persistence and cycling  status - Figures S2A-B
mama@meta.data$stage.diff.cc <- NA
mama@meta.data$stage.diff.cc <- paste0(mama@meta.data$long_vs_short, "_", mama@meta.data$cc.status)
color.stg.cc <- setNames(c("#CCCCFF", "#4728E6", "#FF8100", "#FD1212", "#86B0DA", "#08519c", "#40E0D0", "#3B5A57"), c("<24hr_cycling", "<24hr_non-cycling", "24-36hr_cycling", "24-36hr_non-cycling", "36-48hr_cycling", "36-48hr_non-cycling", ">48hr_cycling", ">48hr_non-cycling"))

DimPlot(mama, group.by = "stage.diff.cc", cols = color.stg.cc)

##Just plot the cycling states
color.stg.cycling <- setNames(c("#D35FB7", "#40B0A6", "#FFC20A", "#005AB5"), c("<24hr_cycling", "24-36hr_cycling", "36-48hr_cycling", ">48hr_cycling"))

color.stg.non.cycling <- setNames(c("#D35FB7", "#40B0A6", "#FFC20A", "#005AB5"), c("<24hr_non-cycling", "24-36hr_non-cycling", "36-48hr_non-cycling", ">48hr_non-cycling"))


dpi <- 300
png(paste0(save.path.scaled.small, "/figures/mama_stage_diff_cc_state_eps_81.png"), height = 5*dpi, width = 5*dpi)
DimPlot(mama, group.by = "stage.diff.cc", cols = color.stg.cc, raster = F, na.value = "#f5f7f6") + NoAxes()
dev.off()

##Figure S2A
dpi <- 300
png(paste0(save.path.scaled.small, "/figures/mama_stage_diff_cycling_eps_81.png"), height = 5*dpi, width = 5*dpi)
DimPlot(mama, group.by = "stage.diff.cc", cols = color.stg.cycling, raster = F, na.value = "#f5f7f6") + NoAxes() + NoLegend()
dev.off()

##Figure S2B
dpi <- 300
png(paste0(save.path.scaled.small, "/figures/mama_stage_diff_non-cycling_eps_81.png"), height = 5*dpi, width = 5*dpi)
DimPlot(mama, group.by = "stage.diff.cc", cols = color.stg.non.cycling, raster = F, na.value = "#f5f7f6") + NoAxes() + NoLegend()
dev.off()

##################################################################################################################################################
##For downstream analyses, we classified cells into “short-term” and “long-term” states based on a threshold of 36 hours—transcriptional states present ≥36 hours were considered “long-term” (Figure 2B). This threshold was chosen to balance focusing on states whose duration was rare while avoiding approaching the upper limit possible for this analysis on a time course of 120 hours. Additionally, each cell was classified as “cycling” or “non-cycling” based on its expression of transcripts associated with different cell cycle phases. 

##For similar analysis using other thresholds, see the optimization of threshold script

##################################################################################################################################################

##Now use a threshold of 24hours - cells that have a stage difference more than 24hours will be categorized as long-term states and ones which have stage difference less than 24hours will be defined short term
##Define a slot in mama for long-term and short-term states
mama@meta.data$states.above.24 <- NA
mama@meta.data[which(mama@meta.data$stage.diff >= 24 & mama@meta.data$cc.status == "cycling"), "states.above.24"] <- "long-term_cycling"
mama@meta.data[which(mama@meta.data$stage.diff >= 24 & mama@meta.data$cc.status == "non-cycling"), "states.above.24"] <- "long-term_non-cycling"
mama@meta.data[which(mama@meta.data$stage.diff <= 24 & mama@meta.data$cc.status == "cycling"), "states.above.24"] <- "short-term_cycling"
mama@meta.data[which(mama@meta.data$stage.diff <= 24 & mama@meta.data$cc.status == "non-cycling"), "states.above.24"] <- "short-term_non-cycling"

##Generate the UMAP plot for Figure 2B
colors.threshold <- setNames(c("#E9967A", "#810000", "#08519c", "#40E0D0"), c("long-term_cycling", "long-term_non-cycling", "short-term_cycling", "short-term_non-cycling"))
dpi <- 300
png(paste0(save.path.scaled.small, "/figures/mama_states_above_24_eps_81.png"), height = 5*dpi, width = 5*dpi)
DimPlot(mama, group.by = "states.above.24", cols = colors.threshold, raster = F) + NoAxes() + NoLegend()
dev.off()

##Save the persistent and short-term states assignment for 24 hours as a threshold
saveRDS(mama@meta.data$states.above.24, paste0(save.path.scaled.small, "/figures/mama_cell_states_threshold_24_eps_81.rds"))
states.above.24 <- readRDS(paste0(save.path.scaled.small, "/figures/mama_cell_states_threshold_30_eps_81.rds"))
mama@meta.data$states.above.24 <- states.above.24

##Subset a dataframe from mama metadata
##Load in the precalculated cell cycle score on mama
cc.score <- readRDS("~/Box/zfext/02-Clustering/mama+ds_integrated/merged_mama_cellCycle_score.rds")
mama@meta.data <- cbind(mama@meta.data, cc.score)

##Make a column in mama meta.data with the maximum score between S and G2M scores
library(dplyr)
mama@meta.data$cc.score <- pmax(mama@meta.data$s.score, mama@meta.data$g2m.score)

##Subset dataframe
dat <- mama@meta.data[,c("stage.nice", "stage.group", "tissue.subsets", "clust", "clust.sg", "cc.phase", "long_vs_short",  "cc.status", "stage.diff.cc", "states.above.24", "stage.diff", "cc.score")]
p <- DimPlot(mama, group.by = "stage.diff.cc", raster = F) + NoAxes()
color.stg.cc <- unique(ggplot_build(p)$data[[1]]$colour)

# Load cluster annotations 
is.empty <- function(x) {
  is.na(x) | is.null(x) | (sapply(trimws(x), nchar) == 0)
}
annot.confirmed <- read.csv(paste0(base.path, "mama_annotations_full_extended_jeffedit_newtissues_curated_clustering.csv"))
annot.confirmed <- annot.confirmed[!is.empty(annot.confirmed$clust),3:ncol(annot.confirmed)]
rownames(annot.confirmed) <- annot.confirmed$clust
annot.confirmed$confirmed <- T

##Merge the annotations to the dataframe
colnames(annot.confirmed$UI) <- "clust"
dat.annot <- merge(dat, annot.confirmed[, c("clust", "identity.super")], by = "clust")
dat <-  dat.annot
dat$tissue <- unlist(lapply(strsplit(x = dat$tissue.stage, split = "_"), function(x) x[1]))

##Save the dataframe
saveRDS(dat, paste0(save.path.scaled.small, "long-term_cycling_scaled_24h_thres_0.10_eps_81_df_cell_ids.rds"))
saveRDS(dat.annot, paste0(save.path.scaled.small, "long-term_cycling_scaled_24h_thres_0.10_eps_81_df.rds"))


##Find the clusters that encompass the long term states for a threshold of 24 hours
dat.annot <- readRDS(paste0(save.path.scaled.small, "long-term_cycling_scaled_30h_thres_0.10_eps_81_df_cell_ids.rds"))
dat <- readRDS(paste0(save.path.scaled.small, "long-term_cycling_scaled_30h_thres_0.10_eps_81_df.rds"))

##Get the long-term states at a time-threshold of 24 hours ## Done on the cluster
time.thresh <- 24

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

# Identify clusts with long-term cycling or long-term non-cycling states
thresh.states <- 0.10
clusts.ltc <- rownames(cluster.proportions)[cluster.proportions$long_cycling >= thresh.states]
clusts.ltnc <- rownames(cluster.proportions)[cluster.proportions$`long_non-cycling` >= thresh.states]

# Only operate tissues that have a long-term state
clusts.with.lt <- unique(unlist(list(clusts.ltc, clusts.ltnc)))
tissues.with.lt <- sort(unique(URD::cutstring(c(clusts.ltc, clusts.ltnc), delimiter = "\\.", field = 1)))

# Identify members of each cluster
members <- lapply(clusts.ltc, function(c) rownames(dat)[which(dat$clust == c & dat$states.above.24 == "long-term_cycling")])
names(members) <- clusts.ltc

##Identify neighbors for these clusters that are long-term cycling
# Initialize empty list to hold neighbors results across all tissues
neighbors <- vector("list", length(clusts.ltc))
names(neighbors) <- clusts.ltc

# Loop per tissue because neighbor lists are per-tissue
tissues.with.lt <- c("epidermis", "eye", "fin", "glial", "hematopoietic", "mesenchyme", "ionocytes", "muscle", "neural", "otic", "periderm", "pigment", "taste")
for (tissue in tissues.with.lt) {
  # Load nn file
  nn.scaled <- readRDS(paste0(neighbor.path.scaled, tissue, "_neighbor_matrix_", eps, "_binary.rds"))
  
  # Identify clusters from this tissue
  tissue.short <- substr(tissue, 1, 4)
  clusts.do <- grep(tissue.short, clusts.ltc, value = T)
  
  # Identify the neighbors of each cluster
  for (c in clusts.do) {
    neighbors[[c]] <- names(which(colSums(nn.scaled[members[[c]], ]) > 0))
  }
}

saveRDS(neighbors, paste0(save.path.scaled.small, "long-term_cyc_clust_with_neighbors_scaled_eps_81_all_pct_10_thres_24hr.rds"))
saveRDS(members, paste0(save.path.scaled.small, "long-term_cyc_clust_with_members_scaled_eps_81_all_pct_10.rds"))

##Load the full members and nighbors list
neighbors <- readRDS(paste0(save.path.scaled.small, "long-term_cyc_clust_with_neighbors_scaled_eps_81_all_pct_10.rds"))
members <- readRDS(paste0(save.path.scaled.small, "long-term_cyc_clust_with_members_scaled_eps_81_all_pct_10.rds"))

##Get all the clusters that the neighbors belong to
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
names(neighbor.clust) <- sort(clusts.ltc)

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

write.table(dat.neighbors, paste0(save.path.scaled.small, "neighbor_clusts_ltc_scaled_eps_81_thres_30h_pct_10.tsv"), sep = "\t")

dat.neighbors <- read.table(paste0(save.path.scaled.small, "neighbor_clusts_ltc_scaled_eps_81_thres_30h_pct_10.tsv"), sep = "")

##Now just limit the neighbors to long-term cycling states to find overlap between the neighbors and the long-term states
dat.neighbors.ltc <- dat.neighbors[which(dat.neighbors$neighbor.clust %in% clusts.ltc), ]

##Before doing the merge, make a copy of the original member and neighbor variables
members.orig <- members
neighbors.orig <- neighbors

##Merge long-term cycling clusters that  are neighbors of each other and have >= 99% overlap
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

##Tissues in which clusters need to be merged
clusts.keep <- list()
tissues.merge <- c("eye", "hematopoietic")
##Now merge the cells present in these members and neighbors
for(tissue in tissues.merge){
  # Need to merge the members and neighbors within each tissue
  tissue.short <- substr(tissue, 1, 4)
  ##Find clusters within the ones that need to be merged
  clusts.to.merge <- grep(tissue.short, clusts.merge, value = T)
  
  ##Get members present in these clusters
  cells.members <- unique(unlist(lapply(clusts.to.merge, function(x) members[[x]])))
  ##Get neighbors present in these clusters
  cells.neighbors <- unique(unlist(lapply(clusts.to.merge, function(x) neighbors[[x]])))
  
  ##Keep only one of the clusters out of the ones to merge
  clust.keep <- unlist(clusts.to.merge[[1]])
  clusts.remove <- setdiff(clusts.to.merge, clust.keep)
  
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
saveRDS(neighbors, paste0(save.path.scaled.small, "long-term_clust_with_neighbors_scaled_eps_81_merged_pct_0.99_merged_thres_24h.rds"))
saveRDS(members, paste0(save.path.scaled.small, "long-term_clust_members_scaled_eps_81_merged_pct_0.99_merged_thres_24h.rds"))

neighbors <- readRDS(paste0(save.path.scaled.small, "long-term_clust_with_neighbors_scaled_eps_81_merged_pct_0.99_merged.rds"))
members <- readRDS(paste0(save.path.scaled.small, "long-term_clust_members_scaled_eps_81_merged_pct_0.99_merged.rds"))

##Subset the original dataframe such that the merged clusters are removed
dat.neighbors <- dat.neighbors[which(dat.neighbors$ltc_clust %in% names(members)), ]

#Filter neighbors such that all clusters where less than 50% of cells are detected as neighbors are discarded
#dat.neighbors.lt <- dat.neighbors[which(dat.neighbors$nn_prop >= 0.5), ]

##From our annotations we identified several clusters as doublets - removing those clusters from the analysis
##Get which clusters are doublets
clust.doublets <- annot.confirmed[which(annot.confirmed$identity.super %in% grep("doublets", annot.confirmed$identity.super, value = T)), "clust"]
clust.doublets <- unlist(clust.doublets$clust)

##Remove all clusters from the long-term states and neighbors that correspond to these doublets
dat.neighbors.lt.clean <- dat.neighbors[!(dat.neighbors$neighbor.clust %in% clust.doublets),]
dat.neighbors.lt.clean <- dat.neighbors[!(dat.neighbors$ltc_clust %in% clust.doublets),]

##In some cases, there are examples where early dividing cells have neighbors in their differentiated counterparts. This is a caveat of this analysis. 
##Identify these clusters and their corresponding neighbors and remove them from analysis (e.g., primitive erythroblasts)
clusts.keep <- list()
clusts.to.remove <- list()
for(cluster in unique(dat.neighbors.lt.clean$ltc_clust)){
  cells.member <- members[[cluster]]
  cells.neighbors <- neighbors[[cluster]]
  stage.member <- as.data.frame(table(mama@meta.data[cells.member, "stage.nice"]))
  stage.neighbor <- as.data.frame(table(mama@meta.data[cells.neighbors, "stage.nice"]))
  member.cc <- as.data.frame(table(mama@meta.data[cells.member, "cc.status"]))
  neighbor.cc <- as.data.frame(table(mama@meta.data[cells.neighbors, "cc.status"]))
  colnames(stage.member) <- c("stage", "cell_num")
  colnames(stage.neighbor) <- c("stage", "cell_num")
  colnames(member.cc) <- c("cc.status", "cell_num")
  colnames(neighbor.cc) <- c("cc.status", "cell_num")
  
  stage.member$cell_prop <- stage.member$cell_num/sum(stage.member$cell_num)
  stage.neighbor$cell_prop <- stage.neighbor$cell_num/sum(stage.neighbor$cell_num)
  member.cc$cell_prop <- member.cc$cell_num/sum(member.cc$cell_num)
  neighbor.cc$cell_prop <- neighbor.cc$cell_num/sum(neighbor.cc$cell_num)
  stage.member.good <- stage.member[which(stage.member$cell_prop >= 0.03), ]
  stage.neighbor.good <- stage.neighbor[which(stage.neighbor$cell_prop >= 0.03), ]
  
  
  ##Now find the farthest stage that cells belong to in the member and neighbor table
  member.stage <- as.numeric(levels(stage.member.good$stage)[as.integer(stage.member.good$stage)])
  neighbor.stage <- as.numeric(levels(stage.neighbor.good$stage)[as.integer(stage.neighbor.good$stage)])
  mem.stg.min <- median(member.stage)
  nbr.stg.min <- median(neighbor.stage)
  
  nbr.proportion.cycling <- neighbor.cc[neighbor.cc$cc.status == "cycling", "cell_prop"]
  
  if(nbr.stg.min - mem.stg.min >= 30){
    if(nbr.proportion.cycling <= 0.75){
      clusts.to.remove <- unlist(list(clusts.to.remove, cluster))
    }
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

write.table(dat.neighbors.lt.filtered, paste0(save.path.scaled.small, "neighbor_clusts_long-term_states_scaled_eps_81_thres_24h_filtered_pct_10.tsv"))
dat.neighbors.lt.filtered <- read.table(paste0(save.path.scaled.small, "neighbor_clusts_long-term_states_scaled_eps_81_thres_24h_filtered_pct_10.tsv"), sep = "")

##Now load the neighbor lists and overlap matrix from the files saved and transferred from the cluster
##neighbors <- readRDS("~/Box/zfext/global_analysis_curated/transient_vs_long-term/long-term_clust_with_neighbors_eps_35_all.rds")

##Ok take one set of neighbors - hema.6
#nbr <- neighbors$hema.14

##Subset dat dataframe
#dat.sub <- dat[nbr, ]
#dat.sub <- dat.sub[which(dat.sub$clust == clusts.with.lt), ]

##Find the range of stages within which 81% of the neighbors fall
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

for(tissue in tissues.with.lt) {
  
  # Identify clusters from this tissue
  tissue.short <- substr(tissue, 1, 4)
  clusts.do <- grep(tissue.short, clusts.ltc.filtered, value = T)
  
  ##Get neighbors for these clusters
  for (c in clusts.do){
    neighbr <- neighbors[[c]]
    dat.sub <- dat[neighbr, ]
    dat.sub <- dat.sub[which(dat.sub$clust %in% clusts.with.lt), ]
    
    ##Find the range of stages within which 80% of the neighbors fall
    low_bound <- quantile(dat.sub$stage.nice, 0.1) ##low bound of range
    high_bound <- quantile(dat.sub$stage.nice, 0.9) ##high bound of range
    
    ##Extract the elements that fall within that range
    elements_in_range <- dat.sub[which(dat.sub$stage.nice >= low_bound & dat.sub$stage.nice <= high_bound), ]
    
    plot.df[which(plot.df$clusts.ltc.filtered == c), "stage.start"] <- low_bound
    plot.df[which(plot.df$clusts.ltc.filtered == c), "stage.stop"] <- high_bound
  }
}


##Add the annotations
colnames(plot.df) <- c("clust", "tissue.subsets", "stage.start", "stage.stop")
# Load cluster annotations 
is.empty <- function(x) {
  is.na(x) | is.null(x) | (sapply(trimws(x), nchar) == 0)
}
annot.confirmed <- readxl::read_xlsx("~/Box/zfext/02-Clustering/mama_tissue_clustering/annot/mama_annotations_full_extended_jeffedit_newtissues_curated_clustering.xlsx")
annot.confirmed <- annot.confirmed[!is.empty(annot.confirmed$clust),3:ncol(annot.confirmed)]
rownames(annot.confirmed) <- annot.confirmed$clust
annot.confirmed$confirmed <- T
plot.df.annot <- merge(plot.df, annot.confirmed[, c("clust", "identity.super", "identity.sub")], by = "clust")

##Plot all long-term states - segment plot shown in Figure 2F
pdf(file = paste0(save.path.scaled.small, "/figures/cell_states_time_lengths_cyc_scaled_eps_81_thres_30h.pdf"), width = 28, height = 36)
ggplot(data=plot.df.annot)+
  geom_segment(aes(x=stage.start, xend=stage.stop, y = clust, yend = clust),lwd=4)
dev.off()


##To calculate what percentage of cells belongs to each category, we also plotted area plots shown in Figure S2D
##Plot area plots for all tissues
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

saveRDS(df.good, paste0(save.path.scaled.small, "tissue_stage_scaled_eps_81_thres_0.1_area_plots.rds"))

##Read the above dataframe
df.good <- readRDS(paste0(save.path.scaled.small, "tissue_stage_scaled_eps_81_thres_0.1_area_plots.rds"))

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
         device = "pdf", path = paste0(save.path.scaled.small, "/figures/"))
  
}

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

saveRDS(df.good.norm, file = paste0(save.path.scaled.small, "tissue_stage_dataframe_area_plots_eps_81.rds"))

area.colors <- setNames(c("#DAA520", "#B8860B", "#FFA500", "#F08181", "#87CEEB", "#7FFFD4", "#00CED1", "#008181"), c("<24hr_cycling", "24-36hr_cycling", "36-48hr_cycling", ">48hr_cycling", "<24hr_non-cycling", "24-36hr_non-cycling", "36-48hr_non-cycling", ">48hr_non-cycling"))
df.good.norm$cc.status <- str_split_fixed(df.good.norm$stage.diff, pattern = "_", 2)
df.good.norm$cc.status <- as.character(df.good.norm$cc.status[,2])
colnames(df.good.norm) <- c("stage", "stage.diff", "proportion", "tissue", "cc.phase")


df.good.norm.order <- df.good.norm[order(df.good.norm$stage.diff, decreasing = FALSE),]


png(file = paste0(save.path.scaled.small, "/figures/tissue_stage_area_plots_v3.png"))
ggplot(data = df.good.norm, aes(x = stage, y = proportion, fill = interaction(stage.diff, cc.phase))) + geom_area() + facet_wrap(vars(tissue), ncol = 4, strip.position = "right") + theme_bw() + scale_fill_manual(values = c("#DAA520", "#B8860B", "#FFA500", "#FF4500", "#87CEEB", "#7FFFD4", "#00CED1", "#008181"))
dev.off()

##Now for some known cell types, using the scaled data and eps of 81, plot the stage_diff vs cell_cycle score plots

##Figure 2C, 2D, 2E, 2G-I,  S2E-I
############################## MUSCLE ##################################################
##Plot scatter plots and the neighbor distribution plots for the well known subsets - fast-muscle, slow-muscle, and satellite cells
##Satellite cells are in the muscle dataset, so load muscle dataset
tissue <- "muscle"
eps <- "81"
#nn <- readRDS(paste0(neighbor.path.scaled, tissue, "_neighbor_matrix_", eps, "_binary.rds")) #Fill this part in.
obj <- readRDS(paste0("/data/CSD/zfext/results/04f-SubsetsV6/obj_subsets_v6/", tissue, "_seurat.rds"))

##Set idents
Idents(mama) <- mama@meta.data$clust

##Using the function above, plot the same for satellite cells
##Satellite cells 
cells.satellite <- WhichCells(mama, idents = c("musc.4", "musc.19", "musc.9", "musc.15", "musc.7"))
cells.fast.muscle <- WhichCells(mama, idents = c("musc.3", "musc.4", "musc.8", "musc.10", "musc.27", "musc.5", "musc.6"))
cells.slow.muscle <- WhichCells(mama, idents = c("musc.3", "musc.24", "musc.17", "musc.11"))
cells.sclerotome <- WhichCells(mama, idents = c("musc.4", "musc.19", "musc.21", "musc.13", "musc.16"))

##New populations
cells.cephalic <- WhichCells(mama, idents = c("musc.26"))
cells.fibroblasts <- WhichCells(mama, idents = c("musc.9"))

##Subset mama metadata
dm <- mama@meta.data[cells.satellite, ]
dm <- mama@meta.data[cells.fast.muscle, ]
dm <- mama@meta.data[cells.slow.muscle, ]
dm <- mama@meta.data[cells.sclerotome, ]
dm <- mama@meta.data[cells.cephalic, ]
dm <- mama@meta.data[cells.fibroblasts, ]

##Take the max of s.score and g2m.score and save as cc.score
dm$cc.score <- pmax(dm$s.score, dm$g2m.score)

##Plot scatter plot
range(dm$stage.nice)

sg_colors_4.120 <- c("#f4be1d", "#bd8700", "#8B4513", "#FF66FF", "#FF1493", "#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#6600CC", "#330066")
sg_colors_10.12 <- c("#FF66FF", "#FF1493", "#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#6600CC", "#330066")
sg_colors_14.21 <- c("#FF1493", "#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#6600CC", "#330066")
sg_colors_24_120 <- c("#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#330066", "#6600CC", "#4C0099")
sg_colors_48_120 <- c("#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#330066", "#6600CC", "#4C0099")


##Scatter Plot _stage.group"#D2B48C", "#CD853F", "#8B4513"
cell_type <- "satellite-cells"
dpi <- 300
png(paste0("/data/CSD/zfext/LTA/results/13-EPS/figures/tissue_scatter_plots/", cell_type, "_scatter_plot_stage.group_eps_81_v2.png"), height = 5*dpi, width = 5*dpi)
# Create a scatter plot with color based on cluster
ggplot(dm, aes(x = cc.score, y = stage.diff, color = stage.group)) +
  geom_point(alpha = 0.6, size = 7) +
  scale_color_manual(values = sg_colors_14.21) + 
  geom_hline(yintercept = 24, linetype = "dashed", color = "red", size = 4) +  # Add horizontal line
  geom_vline(xintercept = 0.0, linetype = "dashed", color = "blue", size = 4) +  # Add vertical line
  labs(x = "Cell Cycle Score", y = "Stage Difference") +
  ggtitle("Scatter Plot of Cell Cycle Score vs Stage Difference") + theme_minimal() +
  scale_x_continuous(limits = c(-0.5, 1)) + scale_y_continuous(limits = c(0, 80))
dev.off()

##Scatter Plot: cluster
cluster_colors <- c("#00CED1",  "#08519c", "#008000", "#8FBC8F", "#800000", "#00FA9A", "#FF00FF", "#008080", "#DB7093",  "#FF1493",  "#2b8cbe", "#fcb163", "#02818a", "#bdbdbd")
pdf(paste0("/data/CSD/zfext/LTA/results/13-EPS/figures/tissue_scatter_plots/", cell_type, "_scatter_plot_cluster.pdf"))
# Create a scatter plot with color based on cluster
set.seed(43)
ggplot(dm, aes(x = cc.score, y = stage.diff, color = clust)) +
  geom_point(alpha = 0.6, size = 7) +
  scale_color_manual(values = sample(cluster_colors)) + 
  geom_hline(yintercept = 24, linetype = "dashed", color = "red") +  # Add horizontal line
  geom_vline(xintercept = 0.0, linetype = "dashed", color = "blue") +  # Add vertical line
  labs(x = "Cell Cycle Score", y = "Stage Difference") +
  ggtitle("Scatter Plot of Cell Cycle Score vs Stage Difference") + theme_minimal() + 
  scale_x_continuous(limits = c(-0.5, 1)) + scale_y_continuous(limits = c(0, 110))
dev.off()

##Get the stages to plot
stages.to.plot <- unique(dm$stage.group)

##Get all neighbor information
for(stages in stages.to.plot){
  plot_neighbor_stage("muscle", "fibroblasts", stages)
}

####### THE FUNCTION BELOW PLOTS THE EXPRESSION OF ANY GENE ON A CC.SCORE vs STAGE.DIFF scatter plot
##Write a function for the above such that it plots the expression of a gene on this scatter plot
library(ggplot2)
library(ggnewscale)
URD.colors <- c("#B2B2B2", "#9BABC2", "#7D9FD1", "#5A90E0", "#307DF0", "#0065FF", 
                "#0078FF", "#008DFF", "#00A1FF", "#00B5FF", "#00CAFF", "#00DEFF", 
                "#00F2FF", "#27FFD7", "#8CFF71", "#F1FF0D", "#FFEE00", "#FFDB00", 
                "#FFC900", "#FFB700", "#FFA500", "#FF9200", "#FF8000", "#FF6D00", 
                "#FF5B00", "#FF4800", "#FF3000", "#FF2400", "#FF1200", "#FF0000"
)

plot_gene_scatter <- function(gene){
  ##Retrieve gene expression from expression data matrix
  gene.expression <- as.data.frame(obj@assays$RNA@data[gene, rownames(dm)])
  ##Rename column names
  colnames(gene.expression) <- "exp_fold_change"
  
  ##Subset dataframe from global dataset with cc.scores and stage.diff information already present
  dt <- dm
  dt$exp_fold_change <- gene.expression$exp_fold_change
  
  ##Set scale and plot the gene expression
  # Create the scatter plot
  ##Define the URD colors as the default paletter
  color_palette <- URD.colors
  
  ##Set alpha
  dt$alpha <- 1
  dt[dt$exp_fold_change < 0.5, "alpha"] <- 1 - (0.8/0.5) * (0.5 - dt[dt$exp_fold_change < 0.5, "exp_fold_change"])
  
  ##Plot the scatter plot with gene expression
   ggplot(dt, aes(x = cc.score, y = stage.diff, color = exp_fold_change)) +
    geom_point(alpha = dt$alpha, size = 7) +
    scale_color_gradientn(colors = URD.colors) +
    #scale_color_identity(guide = "legend", breaks = expression_values, labels = expression_values) +
    geom_hline(yintercept = 24, linetype = "dashed", color = "red", size = 3) +  # Add horizontal line
    geom_vline(xintercept = 0.0, linetype = "dashed", color = "blue", size = 3) +  # Add vertical line
    labs(x = "Cell Cycle Score", y = "Stage Difference", color = "Gene\nExpr.\n(log2)") +
    ggtitle(paste0(gene, "_expression_scaled_eps_81")) + theme_minimal() +
    scale_x_continuous(limits = c(-0.5, 1)) + scale_y_continuous(limits = c(0, 80))
}

genes.to.plot <- c("meox1", "myf5", "myhz1.1", "myhz1.2", "myhz1.3", "mymk", "myog", "pax3a", "pax3b", "pax7a", 
                   "tbx16l", "tbx6", "tcf15", "tnni2b.2", "twist1a", "ripply2", "ripply1")
dpi <- 300
gene.plot <- "tcf21"
png(paste0("/data/CSD/zfext/LTA/results/13-EPS/figures/expression_plots_eps_81/", cell_type, "_expression_plot_", gene.plot, ".png"), width = 5*dpi, height = 5*dpi)
plot_gene_scatter(gene.plot)
dev.off()

################## EYE #########################
sample <- "eye"
tissue <- "eye"
nn.data <- readRDS(paste0(save.path.data, sample, "_neighbor_matrix_35_binary.rds")) #Fill this part in.
obj <- readRDS(paste0("/data/CSD/zfext/results/04f-SubsetsV6/obj_subsets_v6/", sample, "_seurat.rds"))

##Using the function above, plot the same for satellite cells
##Get cells for eye cells
cells.optic.field <- WhichCells(mama, idents = c("eye.17", "eye.26", "eye.2", "eye.37", "eye.11", "eye.33", "eye.3", "eye.18", "eye.31", "eye.16", "eye.23", "eye.13", "eye.38"))
cells.photoreceptors <- WhichCells(mama, idents = c("eye.16", "eye.20", "eye.14", "eye.5", "eye.1", "eye.28"))
cells.rgc <- WhichCells(mama, idents = c("eye.23"))
cells.rpe <- WhichCells(mama, idents = c("eye.22", "eye.41", "eye.35"))
cells.muller <- WhichCells(mama, idents = c("eye.3", "eye.29"))
cells.lens <- WhichCells(mama, idents = c("eye.34"))

##Subset mama metadata
dm <- mama@meta.data[cells.optic.field, ]
dm <- mama@meta.data[cells.rgc, ]
dm <- mama@meta.data[cells.muller, ]
dm <- mama@meta.data[cells.lens, ]

##Take the max of s.score and g2m.score and save as cc.score
dm$cc.score <- pmax(dm$s.score, dm$g2m.score)

##Plot scatter plot
range(dm$stage.nice)

######## PLOTTING NEIGHBORS AND QUERY CELLS FOR VARIOUS EYE SUBTYPES #################################
cell_type <- "muller_glia"

##Get the stages to plot
stages.to.plot <- unique(dm$stage.group)
stages.to.plot <- unique(dm$stage.group)[2:12]


URD.colors <- c("#B2B2B2", "#9BABC2", "#7D9FD1", "#5A90E0", "#307DF0", "#0065FF", 
                "#0078FF", "#008DFF", "#00A1FF", "#00B5FF", "#00CAFF", "#00DEFF", 
                "#00F2FF", "#27FFD7", "#8CFF71", "#F1FF0D", "#FFEE00", "#FFDB00", 
                "#FFC900", "#FFB700", "#FFA500", "#FF9200", "#FF8000", "#FF6D00", 
                "#FF5B00", "#FF4800", "#FF3600", "#FF2400", "#FF1200", "#FF0000"
)

max(dm$stage.diff)

dpi <- 300
png(paste0("/data/CSD/zfext/LTA/results/13-EPS/figures/tissue_scatter_plots/", cell_type, "_scatter_plot_stage.group_eps_81.png"), height = 5*dpi, width = 5*dpi)
# Create a scatter plot with color based on cluster
ggplot(dm, aes(x = cc.score, y = stage.diff, color = stage.group)) +
  geom_point(alpha = 0.6, size = 7) +
  scale_color_manual(values = sg_colors_14.21) + 
  geom_hline(yintercept = 24, linetype = "dashed", color = "red", size = 4) +  # Add horizontal line
  geom_vline(xintercept = 0.0, linetype = "dashed", color = "blue", size = 4) +  # Add vertical line
  labs(x = "Cell Cycle Score", y = "Stage Difference") +
  ggtitle("Scatter Plot of Cell Cycle Score vs Stage Difference") + theme_minimal() +
  scale_x_continuous(limits = c(-0.5, 1)) + scale_y_continuous(limits = c(0, 80))
dev.off()

##Scatter Plot: cluster
cluster_colors <- c("#00CED1",  "#08519c", "#008000", "#8FBC8F", "#800000", "#00FA9A", "#FF00FF", "#008080", "#DB7093",  "#FF1493",  "#2b8cbe", "#fcb163", "#02818a", "#bdbdbd")
pdf(paste0("/data/CSD/zfext/LTA/results/13-EPS/figures/tissue_scatter_plots/", cell_type, "_scatter_plot_cluster.pdf"))
# Create a scatter plot with color based on cluster
set.seed(43)
ggplot(dm, aes(x = cc.score, y = stage.diff, color = clust)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = sample(cluster_colors)) + 
  geom_hline(yintercept = 24, linetype = "dashed", color = "red") +  # Add horizontal line
  geom_vline(xintercept = 0.0, linetype = "dashed", color = "blue") +  # Add vertical line
  labs(x = "Cell Cycle Score", y = "Stage Difference") +
  ggtitle("Scatter Plot of Cell Cycle Score vs Stage Difference") + theme_minimal() + 
  scale_x_continuous(limits = c(-0.5, 1)) + scale_y_continuous(limits = c(0, 80))
dev.off()

dpi <- 300
gene.plot <- "igf2bp2b"
png(paste0("/data/CSD/zfext/LTA/results/13-EPS/figures/expression_plots_eps_81/", cell_type, "_expression_plot_", gene.plot, ".png"), width = 5*dpi, height = 5*dpi)
plot_gene_scatter(gene.plot)
dev.off()


################## EPIDERMIS #########################
sample <- "epidermis"
tissue <- "epidermis"
#nn.data <- readRDS(paste0(save.path.data, sample, "_neighbor_matrix_35_binary.rds")) #Fill this part in.
obj <- readRDS(paste0("/data/CSD/zfext/results/04f-SubsetsV6/obj_subsets_v6/", sample, "_seurat.rds"))

##Using the function above, plot the same for satellite cells
##Get cells for eye cells
cells.taste <- WhichCells(mama, idents = c("epid.21"))
cells.epidermis.posterior <- WhichCells(mama, idents = c("epid.27", "epid.12", "epid.26", "epid.7", "epid.14", "epid.18", "epid.13", "epid.3", "epid.1"))

##Subset mama metadata
dm <- mama@meta.data[cells.taste, ]
dm <- mama@meta.data[cells.epidermis.posterior, ]

##Take the max of s.score and g2m.score and save as cc.score
dm$cc.score <- pmax(dm$s.score, dm$g2m.score)

##Plot scatter plot
range(dm$stage.nice)

######## PLOTTING NEIGHBORS AND QUERY CELLS FOR VARIOUS EYE SUBTYPES #################################
cell_type <- "taste"

##Get the stages to plot
stages.to.plot <- unique(dm$stage.group)
stages.to.plot <- unique(dm$stage.group)[2:12]

sg_colors_10.12 <- c("#f4be1d", "#bd8700", "#FF66FF", "#FF1493",  "#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#330066")
sg_colors_14.21 <- c("#bd8700", "#FF66FF", "#FF1493",  "#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#330066")
sg_colors_24_120 <- c("#FF66FF", "#FF1493",  "#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#330066")
sg_colors_36_120 <- c("#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#330066")
sg_colors_48_120 <- c("#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#330066", "#6600CC", "#4C0099")



URD.colors <- c("#B2B2B2", "#9BABC2", "#7D9FD1", "#5A90E0", "#307DF0", "#0065FF", 
                "#0078FF", "#008DFF", "#00A1FF", "#00B5FF", "#00CAFF", "#00DEFF", 
                "#00F2FF", "#27FFD7", "#8CFF71", "#F1FF0D", "#FFEE00", "#FFDB00", 
                "#FFC900", "#FFB700", "#FFA500", "#FF9200", "#FF8000", "#FF6D00", 
                "#FF5B00", "#FF4800", "#FF3600", "#FF2400", "#FF1200", "#FF0000"
)

max(dm$stage.diff)

dpi <- 300
pdf(paste0("/data/CSD/zfext/LTA/results/13-EPS/figures/tissue_scatter_plots/", cell_type, "_scatter_plot_stage.group_eps_81_v2.pdf"), height = 26, width = 34)
# Create a scatter plot with color based on cluster
ggplot(dm, aes(x = cc.score, y = stage.diff, color = stage.group)) +
  geom_point(alpha = 0.6, size = 7) +
  scale_color_manual(values = sg_colors_4.120) + 
  geom_hline(yintercept = 24, linetype = "dashed", color = "red", size = 4) +  # Add horizontal line
  geom_vline(xintercept = 0.0, linetype = "dashed", color = "blue", size = 4) +  # Add vertical line
  labs(x = "Cell Cycle Score", y = "Stage Difference") +
  ggtitle("Scatter Plot of Cell Cycle Score vs Stage Difference") + theme_minimal() +
  scale_x_continuous(limits = c(-0.5, 1)) + scale_y_continuous(limits = c(0, 65))
dev.off()

##Scatter Plot: cluster
cluster_colors <- c("#00CED1",  "#08519c", "#008000", "#8FBC8F", "#800000", "#00FA9A", "#FF00FF", "#008080", "#DB7093",  "#FF1493",  "#2b8cbe", "#fcb163", "#02818a", "#bdbdbd")
pdf(paste0("/data/CSD/zfext/LTA/results/13-EPS/figures/tissue_scatter_plots/", cell_type, "_scatter_plot_cluster.pdf"))
# Create a scatter plot with color based on cluster
set.seed(43)
ggplot(dm, aes(x = cc.score, y = stage.diff, color = clust)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = sample(cluster_colors)) + 
  geom_hline(yintercept = 24, linetype = "dashed", color = "red") +  # Add horizontal line
  geom_vline(xintercept = 0.0, linetype = "dashed", color = "blue") +  # Add vertical line
  labs(x = "Cell Cycle Score", y = "Stage Difference") +
  ggtitle("Scatter Plot of Cell Cycle Score vs Stage Difference") + theme_minimal() + 
  scale_x_continuous(limits = c(-0.5, 1)) + scale_y_continuous(limits = c(0, 65))
dev.off()

plot_gene_scatter <- function(gene){
  ##Retrieve gene expression from expression data matrix
  gene.expression <- as.data.frame(obj@assays$RNA@data[gene, rownames(dm)])
  ##Rename column names
  colnames(gene.expression) <- "exp_fold_change"
  
  ##Subset dataframe from global dataset with cc.scores and stage.diff information already present
  dt <- dm
  dt$exp_fold_change <- gene.expression$exp_fold_change
  
  ##Set scale and plot the gene expression
  # Create the scatter plot
  ##Define the URD colors as the default paletter
  color_palette <- URD.colors
  
  ##Set alpha
  dt$alpha <- 1
  dt[dt$exp_fold_change < 0.5, "alpha"] <- 1 - (0.8/0.5) * (0.5 - dt[dt$exp_fold_change < 0.5, "exp_fold_change"])
  
  ##Plot the scatter plot with gene expression
  ggplot(dt, aes(x = cc.score, y = stage.diff, color = exp_fold_change)) +
    geom_point(alpha = dt$alpha, size = 4) +
    scale_color_gradientn(colors = URD.colors) +
    #scale_color_identity(guide = "legend", breaks = expression_values, labels = expression_values) +
    geom_hline(yintercept = 24, linetype = "dashed", color = "red", size = 3) +  # Add horizontal line
    geom_vline(xintercept = 0.0, linetype = "dashed", color = "blue", size = 3) +  # Add vertical line
    labs(x = "Cell Cycle Score", y = "Stage Difference", color = "Gene\nExpr.\n(log2)") +
    ggtitle(paste0(gene, "_expression_scaled_eps_81")) + theme_minimal() +
    scale_x_continuous(limits = c(-0.5, 1)) + scale_y_continuous(limits = c(0, 80))
}


dpi <- 300
gene.plot <- "frem3"
png(paste0("/data/CSD/zfext/LTA/results/13-EPS/figures/expression_plots_eps_81/", cell_type, "_expression_plot_", gene.plot, ".png"), width = 5*dpi, height = 5*dpi)
plot_gene_scatter(gene.plot)
dev.off()

plot_gene_scatter("sox8b")

################## GLIAL #########################
sample <- "glial"
tissue <- "glial"
#nn.data <- readRDS(paste0(save.path.data, sample, "_neighbor_matrix_35_binary.rds")) #Fill this part in.
obj <- readRDS(paste0("/data/CSD/zfext/results/04f-SubsetsV6/obj_subsets_v6/", sample, "_seurat.rds"))

##Using the function above, plot the same for satellite cells
##Get cells for eye cells
##Using the function above, plot the same for liver, intestine and pancreas
cells.radial.glia <- WhichCells(mama, idents = c("glia.7", "glia.18", "glia.27"))
cells.floor.plate <- WhichCells(mama, idents = c("glia.15"))

##Subset mama metadata
dm <- mama@meta.data[cells.radial.glia, ]
dm <- mama@meta.data[cells.floor.plate, ]

##Take the max of s.score and g2m.score and save as cc.score
dm$cc.score <- pmax(dm$s.score, dm$g2m.score)

##Plot scatter plot
range(dm$stage.nice)

######## PLOTTING NEIGHBORS AND QUERY CELLS FOR VARIOUS EYE SUBTYPES #################################
cell_type <- "radial_glia"

##Get the stages to plot
stages.to.plot <- unique(dm$stage.group)
stages.to.plot <- unique(dm$stage.group)[2:12]

sg_colors_10.12 <- c("#f4be1d", "#bd8700", "#FF66FF", "#FF1493",  "#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#330066")
sg_colors_14.21 <- c("#bd8700", "#FF66FF", "#FF1493",  "#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#330066")
sg_colors_24_120 <- c("#FF66FF", "#FF1493",  "#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#330066")
sg_colors_36_120 <- c("#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#330066")
sg_colors_48_120 <- c("#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#330066", "#6600CC", "#4C0099")



URD.colors <- c("#B2B2B2", "#9BABC2", "#7D9FD1", "#5A90E0", "#307DF0", "#0065FF", 
                "#0078FF", "#008DFF", "#00A1FF", "#00B5FF", "#00CAFF", "#00DEFF", 
                "#00F2FF", "#27FFD7", "#8CFF71", "#F1FF0D", "#FFEE00", "#FFDB00", 
                "#FFC900", "#FFB700", "#FFA500", "#FF9200", "#FF8000", "#FF6D00", 
                "#FF5B00", "#FF4800", "#FF3600", "#FF2400", "#FF1200", "#FF0000"
)

max(dm$stage.diff)

dpi <- 300
png(paste0("/data/CSD/zfext/LTA/results/13-EPS/figures/tissue_scatter_plots/", cell_type, "_scatter_plot_stage.group_eps_81.png"), height = 5*dpi, width = 5*dpi)
# Create a scatter plot with color based on cluster
ggplot(dm, aes(x = cc.score, y = stage.diff, color = stage.group)) +
  geom_point(alpha = 0.6, size = 7) +
  scale_color_manual(values = sg_colors_14.21) + 
  geom_hline(yintercept = 24, linetype = "dashed", color = "red", size = 4) +  # Add horizontal line
  geom_vline(xintercept = 0.0, linetype = "dashed", color = "blue", size = 4) +  # Add vertical line
  labs(x = "Cell Cycle Score", y = "Stage Difference") +
  ggtitle("Scatter Plot of Cell Cycle Score vs Stage Difference") + theme_minimal() +
  scale_x_continuous(limits = c(-0.5, 1)) + scale_y_continuous(limits = c(0, 90))
dev.off()

##Scatter Plot: cluster
cluster_colors <- c("#00CED1",  "#08519c", "#008000", "#8FBC8F", "#800000", "#00FA9A", "#FF00FF", "#008080", "#DB7093",  "#FF1493",  "#2b8cbe", "#fcb163", "#02818a", "#bdbdbd")
pdf(paste0("/data/CSD/zfext/LTA/results/13-EPS/figures/tissue_scatter_plots/", cell_type, "_scatter_plot_cluster.pdf"))
# Create a scatter plot with color based on cluster
set.seed(43)
ggplot(dm, aes(x = cc.score, y = stage.diff, color = clust)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = sample(cluster_colors)) + 
  geom_hline(yintercept = 24, linetype = "dashed", color = "red") +  # Add horizontal line
  geom_vline(xintercept = 0.0, linetype = "dashed", color = "blue") +  # Add vertical line
  labs(x = "Cell Cycle Score", y = "Stage Difference") +
  ggtitle("Scatter Plot of Cell Cycle Score vs Stage Difference") + theme_minimal() + 
  scale_x_continuous(limits = c(-0.5, 1)) + scale_y_continuous(limits = c(0, 65))
dev.off()

plot_gene_scatter <- function(gene){
  ##Retrieve gene expression from expression data matrix
  gene.expression <- as.data.frame(obj@assays$RNA@data[gene, rownames(dm)])
  ##Rename column names
  colnames(gene.expression) <- "exp_fold_change"
  
  ##Subset dataframe from global dataset with cc.scores and stage.diff information already present
  dt <- dm
  dt$exp_fold_change <- gene.expression$exp_fold_change
  
  ##Set scale and plot the gene expression
  # Create the scatter plot
  ##Define the URD colors as the default paletter
  color_palette <- URD.colors
  
  ##Set alpha
  dt$alpha <- 1
  dt[dt$exp_fold_change < 0.5, "alpha"] <- 1 - (0.8/0.5) * (0.5 - dt[dt$exp_fold_change < 0.5, "exp_fold_change"])
  
  ##Plot the scatter plot with gene expression
  ggplot(dt, aes(x = cc.score, y = stage.diff, color = exp_fold_change)) +
    geom_point(alpha = dt$alpha, size = 7) +
    scale_color_gradientn(colors = URD.colors) +
    #scale_color_identity(guide = "legend", breaks = expression_values, labels = expression_values) +
    geom_hline(yintercept = 24, linetype = "dashed", color = "red", size = 3) +  # Add horizontal line
    geom_vline(xintercept = 0.0, linetype = "dashed", color = "blue", size = 3) +  # Add vertical line
    labs(x = "Cell Cycle Score", y = "Stage Difference", color = "Gene\nExpr.\n(log2)") +
    ggtitle(paste0(gene, "_expression_scaled_eps_81")) + theme_minimal() +
    scale_x_continuous(limits = c(-0.5, 1)) + scale_y_continuous(limits = c(0, 90))
}


dpi <- 300
gene.plot <- "slc6a15"
png(paste0("/data/CSD/zfext/LTA/results/13-EPS/figures/expression_plots_eps_81/", cell_type, "_expression_plot_", gene.plot, ".png"), width = 5*dpi, height = 5*dpi)
plot_gene_scatter(gene.plot)
dev.off()

plot_gene_scatter("sox8b")

############ ENDODERM ###########
##Load endoderm dataset
sample <- "endoderm"
tissue <- "endoderm"
#nn.data <- readRDS(paste0(save.path.data, sample, "_neighbor_matrix_35_binary.rds")) #Fill this part in.
obj <- readRDS(paste0("/data/CSD/zfext/results/04f-SubsetsV6/obj_subsets_v6/", sample, "_seurat.rds"))

##Using the function above, plot the same for satellite cells
##Get cells for eye cells
cells.liver <- WhichCells(mama, idents = c("endo.1", "endo.2", "endo.3", "endo.4", "endo.5", "endo.6", "endo.7", "endo.8", "endo.27"))
cells.intestine <- WhichCells(mama, idents = c("endo.9", "endo.10", "endo.11", "endo.12", "endo.13", "endo.14", "endo.15", "endo.16", "endo.17", "endo.18"))

##Subset mama metadata
dm <- mama@meta.data[cells.liver, ]
dm <- mama@meta.data[cells.intestine, ]

##Take the max of s.score and g2m.score and save as cc.score
dm$cc.score <- pmax(dm$s.score, dm$g2m.score)

##Plot scatter plot
range(dm$stage.nice)

sg_colors_10.12 <- c("#f4be1d", "#bd8700", "#FF66FF", "#FF1493",  "#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#330066")
sg_colors_14.21 <- c("#bd8700", "#FF66FF", "#FF1493",  "#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#330066")
sg_colors_24_120 <- c("#FF66FF", "#FF1493",  "#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#330066")
sg_colors_36_120 <- c("#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#330066")
sg_colors_48_120 <- c("#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#330066", "#6600CC", "#4C0099")



URD.colors <- c("#B2B2B2", "#9BABC2", "#7D9FD1", "#5A90E0", "#307DF0", "#0065FF", 
                "#0078FF", "#008DFF", "#00A1FF", "#00B5FF", "#00CAFF", "#00DEFF", 
                "#00F2FF", "#27FFD7", "#8CFF71", "#F1FF0D", "#FFEE00", "#FFDB00", 
                "#FFC900", "#FFB700", "#FFA500", "#FF9200", "#FF8000", "#FF6D00", 
                "#FF5B00", "#FF4800", "#FF3600", "#FF2400", "#FF1200", "#FF0000"
)

max(dm$stage.diff)

cell_type <- "intestine"

dpi <- 300
png(paste0("/data/CSD/zfext/LTA/results/13-EPS/figures/tissue_scatter_plots/", cell_type, "_scatter_plot_stage.group_eps_81_v2.png"), height = 5*dpi, width = 5*dpi)
# Create a scatter plot with color based on cluster
ggplot(dm, aes(x = cc.score, y = stage.diff, color = stage.group)) +
  geom_point(alpha = 0.6, size = 7) +
  scale_color_manual(values = sg_colors_14.21) + 
  geom_hline(yintercept = 24, linetype = "dashed", color = "red", size = 4) +  # Add horizontal line
  geom_vline(xintercept = 0.0, linetype = "dashed", color = "blue", size = 4) +  # Add vertical line
  labs(x = "Cell Cycle Score", y = "Stage Difference") +
  ggtitle("Scatter Plot of Cell Cycle Score vs Stage Difference") + theme_minimal() +
  scale_x_continuous(limits = c(-0.5, 1)) + scale_y_continuous(limits = c(0, 65))
dev.off()

##Scatter Plot: cluster
cluster_colors <- c("#00CED1",  "#08519c", "#008000", "#8FBC8F", "#800000", "#00FA9A", "#FF00FF", "#008080", "#DB7093",  "#FF1493",  "#2b8cbe", "#fcb163", "#02818a", "#bdbdbd")
pdf(paste0("/data/CSD/zfext/LTA/results/13-EPS/figures/tissue_scatter_plots/", cell_type, "_scatter_plot_cluster.pdf"))
# Create a scatter plot with color based on cluster
set.seed(43)
ggplot(dm, aes(x = cc.score, y = stage.diff, color = clust)) +
  geom_point(alpha = 0.6, size = 7) +
  scale_color_manual(values = sample(cluster_colors)) + 
  geom_hline(yintercept = 24, linetype = "dashed", color = "red") +  # Add horizontal line
  geom_vline(xintercept = 0.0, linetype = "dashed", color = "blue") +  # Add vertical line
  labs(x = "Cell Cycle Score", y = "Stage Difference") +
  ggtitle("Scatter Plot of Cell Cycle Score vs Stage Difference") + theme_minimal() + 
  scale_x_continuous(limits = c(-0.5, 1)) + scale_y_continuous(limits = c(0, 65))
dev.off()

plot_gene_scatter <- function(gene){
  ##Retrieve gene expression from expression data matrix
  gene.expression <- as.data.frame(obj@assays$RNA@data[gene, rownames(dm)])
  ##Rename column names
  colnames(gene.expression) <- "exp_fold_change"
  
  ##Subset dataframe from global dataset with cc.scores and stage.diff information already present
  dt <- dm
  dt$exp_fold_change <- gene.expression$exp_fold_change
  
  ##Set scale and plot the gene expression
  # Create the scatter plot
  ##Define the URD colors as the default paletter
  color_palette <- URD.colors
  
  ##Set alpha
  dt$alpha <- 1
  dt[dt$exp_fold_change < 0.5, "alpha"] <- 1 - (0.8/0.5) * (0.5 - dt[dt$exp_fold_change < 0.5, "exp_fold_change"])
  
  ##Plot the scatter plot with gene expression
  ggplot(dt, aes(x = cc.score, y = stage.diff, color = exp_fold_change)) +
    geom_point(alpha = dt$alpha, size = 7) +
    scale_color_gradientn(colors = URD.colors) +
    #scale_color_identity(guide = "legend", breaks = expression_values, labels = expression_values) +
    geom_hline(yintercept = 24, linetype = "dashed", color = "red", size = 3) +  # Add horizontal line
    geom_vline(xintercept = 0.0, linetype = "dashed", color = "blue", size = 3) +  # Add vertical line
    labs(x = "Cell Cycle Score", y = "Stage Difference", color = "Gene\nExpr.\n(log2)") +
    ggtitle(paste0(gene, "_expression_scaled_eps_81")) + theme_minimal() +
    scale_x_continuous(limits = c(-0.5, 1)) + scale_y_continuous(limits = c(0, 65))
}


dpi <- 300
gene.plot <- "scg3"
png(paste0("/data/CSD/zfext/LTA/results/13-EPS/figures/expression_plots_eps_81/", cell_type, "_expression_plot_", gene.plot, ".png"), width = 5*dpi, height = 5*dpi)
plot_gene_scatter(gene.plot)
dev.off()

plot_gene_scatter("sox8b")

############ TASTE-OLFACTORY ###########
##Load endoderm dataset
sample <- "taste"
tissue <- "taste"
#nn.data <- readRDS(paste0(save.path.data, sample, "_neighbor_matrix_35_binary.rds")) #Fill this part in.
obj <- readRDS(paste0("/data/CSD/zfext/results/04f-SubsetsV6/obj_subsets_v6/", sample, "_seurat.rds"))

##Using the function above, plot the same for satellite cells
##Get cells for eye cells
cells.olfactory <- WhichCells(mama, idents = c("tast.4", "tast.12"))

##Subset mama metadata
dm <- mama@meta.data[cells.olfactory, ]

##Take the max of s.score and g2m.score and save as cc.score
dm$cc.score <- pmax(dm$s.score, dm$g2m.score)

##Plot scatter plot
range(dm$stage.nice)

sg_colors_10.12 <- c("#f4be1d", "#bd8700", "#FF66FF", "#FF1493",  "#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#330066")
sg_colors_14.21 <- c("#bd8700", "#FF66FF", "#FF1493",  "#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#330066")
sg_colors_24_120 <- c("#FF66FF", "#FF1493",  "#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#330066")
sg_colors_36_120 <- c("#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#330066")
sg_colors_48_120 <- c("#FFB266", "#FF8000",  "#66CC00", "#00994C", "#0066CC", "#0000FF", "#B266FF", "#330066", "#6600CC", "#4C0099")



URD.colors <- c("#B2B2B2", "#9BABC2", "#7D9FD1", "#5A90E0", "#307DF0", "#0065FF", 
                "#0078FF", "#008DFF", "#00A1FF", "#00B5FF", "#00CAFF", "#00DEFF", 
                "#00F2FF", "#27FFD7", "#8CFF71", "#F1FF0D", "#FFEE00", "#FFDB00", 
                "#FFC900", "#FFB700", "#FFA500", "#FF9200", "#FF8000", "#FF6D00", 
                "#FF5B00", "#FF4800", "#FF3600", "#FF2400", "#FF1200", "#FF0000"
)

max(dm$stage.diff)

cell_type <- "olfactory_epithelium"

dpi <- 300
png(paste0("/data/CSD/zfext/LTA/results/13-EPS/figures/tissue_scatter_plots/", cell_type, "_scatter_plot_stage.group_eps_81_v2.png"), height = 5*dpi, width = 5*dpi)
# Create a scatter plot with color based on cluster
ggplot(dm, aes(x = cc.score, y = stage.diff, color = stage.group)) +
  geom_point(alpha = 0.6, size = 7) +
  scale_color_manual(values = sg_colors_4.120) + 
  geom_hline(yintercept = 24, linetype = "dashed", color = "red", size = 4) +  # Add horizontal line
  geom_vline(xintercept = 0.0, linetype = "dashed", color = "blue", size = 4) +  # Add vertical line
  labs(x = "Cell Cycle Score", y = "Stage Difference") +
  ggtitle("Scatter Plot of Cell Cycle Score vs Stage Difference") + theme_minimal() +
  scale_x_continuous(limits = c(-0.5, 1)) + scale_y_continuous(limits = c(0, 70))
dev.off()

##Scatter Plot: cluster
cluster_colors <- c("#00CED1",  "#08519c", "#008000", "#8FBC8F", "#800000", "#00FA9A", "#FF00FF", "#008080", "#DB7093",  "#FF1493",  "#2b8cbe", "#fcb163", "#02818a", "#bdbdbd")
pdf(paste0("/data/CSD/zfext/LTA/results/13-EPS/figures/tissue_scatter_plots/", cell_type, "_scatter_plot_cluster.pdf"))
# Create a scatter plot with color based on cluster
set.seed(43)
ggplot(dm, aes(x = cc.score, y = stage.diff, color = clust)) +
  geom_point(alpha = 0.6, size = 7) +
  scale_color_manual(values = sample(cluster_colors)) + 
  geom_hline(yintercept = 24, linetype = "dashed", color = "red") +  # Add horizontal line
  geom_vline(xintercept = 0.0, linetype = "dashed", color = "blue") +  # Add vertical line
  labs(x = "Cell Cycle Score", y = "Stage Difference") +
  ggtitle("Scatter Plot of Cell Cycle Score vs Stage Difference") + theme_minimal() + 
  scale_x_continuous(limits = c(-0.5, 1)) + scale_y_continuous(limits = c(0, 70))
dev.off()

plot_gene_scatter <- function(gene){
  ##Retrieve gene expression from expression data matrix
  gene.expression <- as.data.frame(obj@assays$RNA@data[gene, rownames(dm)])
  ##Rename column names
  colnames(gene.expression) <- "exp_fold_change"
  
  ##Subset dataframe from global dataset with cc.scores and stage.diff information already present
  dt <- dm
  dt$exp_fold_change <- gene.expression$exp_fold_change
  
  ##Set scale and plot the gene expression
  # Create the scatter plot
  ##Define the URD colors as the default paletter
  color_palette <- URD.colors
  
  ##Set alpha
  dt$alpha <- 1
  dt[dt$exp_fold_change < 0.5, "alpha"] <- 1 - (0.8/0.5) * (0.5 - dt[dt$exp_fold_change < 0.5, "exp_fold_change"])
  
  ##Plot the scatter plot with gene expression
  ggplot(dt, aes(x = cc.score, y = stage.diff, color = exp_fold_change)) +
    geom_point(alpha = dt$alpha, size = 7) +
    scale_color_gradientn(colors = URD.colors) +
    #scale_color_identity(guide = "legend", breaks = expression_values, labels = expression_values) +
    geom_hline(yintercept = 24, linetype = "dashed", color = "red", size = 3) +  # Add horizontal line
    geom_vline(xintercept = 0.0, linetype = "dashed", color = "blue", size = 3) +  # Add vertical line
    labs(x = "Cell Cycle Score", y = "Stage Difference", color = "Gene\nExpr.\n(log2)") +
    ggtitle(paste0(gene, "_expression_scaled_eps_81")) + theme_minimal() +
    scale_x_continuous(limits = c(-0.5, 1)) + scale_y_continuous(limits = c(0, 65))
}


dpi <- 300
gene.plot <- "matn3b"
png(paste0("/data/CSD/zfext/LTA/results/13-EPS/figures/expression_plots_eps_81/", cell_type, "_expression_plot_", gene.plot, ".png"), width = 5*dpi, height = 5*dpi)
plot_gene_scatter(gene.plot)
dev.off()

plot_gene_scatter("sox8b")

##Get number of cells with a stage difference above 24
cells.long.term <- rownames(mama@meta.data)[mama@meta.data$stage.diff > 24]
cells.long.term.cycling <- rownames(mama@meta.data)
