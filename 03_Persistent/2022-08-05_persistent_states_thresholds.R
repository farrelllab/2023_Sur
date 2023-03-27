##This code is designed to find long-term persistent states and short term states based on arbitrary thresholds
#Load libraries
library(Seurat)
library(dbscan)

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
stage.diff <- readRDS("~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_stage.diff_eps_5NN.rds")
long_vs_short <- readRDS("~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_long_vs_short_stage_diff_5NN.rds")
mama@meta.data$long_vs_short <- long_vs_short
mama@meta.data$stage.diff <- stage.diff


######### Cell states above and below certain thresholds #############

##Start with a threshold of 60 hours - cells that have a stage difference more than 60 hours will be categorized as long-term states and ones which have stage difference less than 60 hours will be defined short term
##Define a slot in mama for long-term and short-term states
mama@meta.data$states.above.60 <- NA
mama@meta.data[which(mama@meta.data$stage.diff >= 60 & mama@meta.data$cc.status == "cycling"), "states.above.60"] <- "long-term_cycling"
mama@meta.data[which(mama@meta.data$stage.diff >= 60 & mama@meta.data$cc.status == "non-cycling"), "states.above.60"] <- "long-term_non-cycling"
mama@meta.data[which(mama@meta.data$stage.diff <= 60 & mama@meta.data$cc.status == "cycling"), "states.above.60"] <- "short-term_cycling"
mama@meta.data[which(mama@meta.data$stage.diff <= 60 & mama@meta.data$cc.status == "non-cycling"), "states.above.60"] <- "short-term_non-cycling"

colors.threshold <- setNames(c("#E9967A", "#800000", "#08519c", "#40E0D0"), c("long-term_cycling", "long-term_non-cycling", "short-term_cycling", "short-term_non-cycling"))
dpi <- 300
png("~/Box/zfext/global_analysis_curated/transient_vs_long-term/figures/mama_cell_states_thresh_60_eps_5NN.png", height = 5*dpi, width = 5*dpi)
DimPlot(mama, group.by = "long_v_short", cols = colors.threshold, raster = F) + NoAxes()
dev.off()

##Save the persistent and short-term states assignment for 60 hours as a threshold
saveRDS(mama@meta.data$states.above.60, "~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_cell_states_threshold_60_eps_5NN.rds")

##Plot a scatterpie plot for these cells that have stage difference above 60 hours
##Load mama metadata slot
states.above.60 <- readRDS("~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_cell_states_threshold_60_eps_5NN.rds")
mama@meta.data$long_v_short <- states.above.60
DimPlot(mama, group.by = "long_v_short", cols = colors.threshold, raster = F) + NoAxes()

long.cycling.cells <- rownames(mama@meta.data)[mama@meta.data$long_v_short == "long-term_cycling"]
long.non_cycling.cells <- rownames(mama@meta.data)[mama@meta.data$long_v_short == "long-term_non-cycling"]
short.cycling.cells <- rownames(mama@meta.data)[mama@meta.data$long_v_short == "short-term_cycling"]
short.non_cycling.cells <- rownames(mama@meta.data)[mama@meta.data$long_v_short == "short-term_non-cycling"]

df <- as.data.frame(table(mama@meta.data$tissue.sg))
colnames(df) <- c("tissue_sg", "total.cells")
lt_diff <- as.data.frame(table(mama@meta.data[long.non_cycling.cells, "tissue.sg"]))
lt_cyc <- as.data.frame(table(mama@meta.data[long.cycling.cells, "tissue.sg"]))
st_diff <- as.data.frame(table(mama@meta.data[short.non_cycling.cells, "tissue.sg"]))
st_cyc <- as.data.frame(table(mama@meta.data[short.cycling.cells, "tissue.sg"]))

colnames(lt_diff) <- c("tissue_sg", "lt_diff")
colnames(lt_cyc) <- c("tissue_sg", "lt_cyc")
colnames(st_diff) <- c("tissue_sg", "st_diff")
colnames(st_cyc) <- c("tissue_sg", "st_cyc")

df.full <- merge(df, lt_diff, all = T)
df.full <- merge(df.full, lt_cyc, all = T)
df.full <- merge(df.full, st_diff, all = T)
df.full <- merge(df.full, st_cyc, all = T)
colnames(df.full) <- c("tissue_sg", "total.cells", "lt_diff", "lt_cyc", "st_diff", "st_cyc")

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


##Split the first column into 2 columns
library(stringr)
df.full[c("tissue", "stage.group")] <- str_split_fixed(df.full$tissue_sg, "_", 2)
##Rearrange columns and remove original column names
#df.full <- df.full[c("tissue", "stage.group", "total.cells", "lt_diff", "lt_cyc", "st_diff", "st_cyc")]
df.prop <- df.full[c("tissue", "stage.group", "total.cells", "lt_diff_prop", "lt_cyc_prop", "st_cyc_prop", "st_diff_prop")]


##Convert the 3 column dataframe into a matrix for a dotplot
library(tidyr)
library(reshape2)
library(viridis)
library(tidyverse)
library(plotly)
library(purrr)

master.df <- df.prop %>% 
  mutate(tissue = recode_factor(tissue,
                                "axia" = "1",
                                "basa" = "2",
                                "blas" = "3",
                                "mesoen" = "21",
                                "endo" = "4",
                                "eye" = "5",
                                "fin" = "6",
                                "gast" = "7",
                                "glia" = "8",
                                "hema" = "9",
                                "iono" = "10",
                                "mese" = "11",
                                "musc" = "12",
                                "neur" = "13",
                                "mura" = "14",
                                "otic" = "15",
                                "peri" = "16",
                                "pgc" = "17",
                                "pigm" = "18",
                                "pron" = "19",
                                "tast" = "20",
                                "unknown" = "22",
                                .ordered = T), 
         tissue = fct_rev(tissue),
         stage = fct_relevel(stage.group)) 

master.df <- master.df %>% mutate(stage = fct_relevel(stage), tissue = fct_rev(tissue))

master.df$stage.nice <- str_split_fixed(master.df$stage, "-", 2)
master.df <- master.df[c("tissue", "stage.nice", "total.cells", "lt_diff_prop", "lt_cyc_prop", "st_diff_prop", "st_cyc_prop")]
master.df$stage <- as.numeric(master.df$stage.nice[,1])
master.df <- master.df[c("tissue", "stage", "lt_cyc_prop", "lt_diff_prop", "st_cyc_prop", "st_diff_prop")]
master.df$tissue <- as.numeric(master.df$tissue)
master.df$stage <- as.numeric(master.df$stage)
master.df$lt_diff_prop <- as.numeric(master.df$lt_diff_prop)
master.df$lt_cyc_prop <- as.numeric(master.df$lt_cyc_prop)
master.df$st_diff_prop <- as.numeric(master.df$st_diff_prop)
master.df$st_cyc_prop <- as.numeric(master.df$st_cyc_prop)

##Plot the scatter pie plots using cell proportions
library(scatterpie)
pdf("mama_cell_states_60h_scatter_pie_eps_5NN.pdf", height = 10, width = 12)
plot.pie <- ggplot() + geom_scatterpie(mapping = aes(x = stage, y = tissue), data = master.df, cols = c("lt_cyc_prop", "lt_diff_prop", "st_cyc_prop", "st_diff_prop"), pie_scale = 0.15, sorted_by_radius = FALSE, size = 0.1) + coord_fixed() + geom_dotplot(binwidth = 1) + labs(y = "Tissue", x = "Hour-Post-Fertilization")
plot(plot.pie)
dev.off()

#Next redo the above script with a threshold of 48 hours - cells that have a stage difference more than 48 hours will be categorized as long-term states and ones which have stage difference less than 48 hours will be defined short term
##Define a slot in mama for long-term and short-term states
mama@meta.data$states.above.48 <- NA
mama@meta.data[which(mama@meta.data$stage.diff >= 48 & mama@meta.data$cc.status == "cycling"), "states.above.48"] <- "long-term_cycling"
mama@meta.data[which(mama@meta.data$stage.diff >= 48 & mama@meta.data$cc.status == "non-cycling"), "states.above.48"] <- "long-term_non-cycling"
mama@meta.data[which(mama@meta.data$stage.diff <= 48 & mama@meta.data$cc.status == "cycling"), "states.above.48"] <- "short-term_cycling"
mama@meta.data[which(mama@meta.data$stage.diff <= 48 & mama@meta.data$cc.status == "non-cycling"), "states.above.48"] <- "short-term_non-cycling"

dpi <- 300
png("~/Box/zfext/global_analysis_curated/transient_vs_long-term/figures/mama_cell_states_thresh_48hr_data_eps_40_v3_with_neural.png", height = 5*dpi, width = 5*dpi)
colors.threshold <- setNames(c("#E9967A", "#800000", "#08519c", "#40E0D0"), c("long-term_cycling", "long-term_non-cycling", "short-term_cycling", "short-term_non-cycling"))
DimPlot(mama, group.by = "states.above.48", cols = colors.threshold, raster = F, na.value = "#f5f7f6") + NoAxes() + NoLegend()
dev.off() 

##Save the persistent and short-term states assignment for 48 hours as a threshold
saveRDS(mama@meta.data$states.above.48, "~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_cell_states_threshold_48_eps_40_with_neural.rds")

##Plot a scatterpie plot for these cells that have stage difference above 48 hours
##Load mama metadata slot
states.above.48 <- readRDS("~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_cell_states_threshold_48_eps_40_except_neural.rds")
mama@meta.data$states.above.48 <- states.above.48
DimPlot(mama, group.by = "states.above.48", raster = F) + NoAxes()

long.cycling.cells <- rownames(mama@meta.data)[which(mama@meta.data$states.above.48 == "long-term_cycling")]
long.non_cycling.cells <- rownames(mama@meta.data)[which(mama@meta.data$states.above.48 == "long-term_non-cycling")]
short.cycling.cells <- rownames(mama@meta.data)[which(mama@meta.data$states.above.48 == "short-term_cycling")]
short.non_cycling.cells <- rownames(mama@meta.data)[which(mama@meta.data$states.above.48 == "short-term_non-cycling")]

df <- as.data.frame(table(mama@meta.data$tissue.sg))
colnames(df) <- c("tissue_sg", "total.cells")
lt_diff <- as.data.frame(table(mama@meta.data[long.non_cycling.cells, "tissue.sg"]))
lt_cyc <- as.data.frame(table(mama@meta.data[long.cycling.cells, "tissue.sg"]))
st_diff <- as.data.frame(table(mama@meta.data[short.non_cycling.cells, "tissue.sg"]))
st_cyc <- as.data.frame(table(mama@meta.data[short.cycling.cells, "tissue.sg"]))

colnames(lt_diff) <- c("tissue_sg", "lt_diff")
colnames(lt_cyc) <- c("tissue_sg", "lt_cyc")
colnames(st_diff) <- c("tissue_sg", "st_diff")
colnames(st_cyc) <- c("tissue_sg", "st_cyc")

df.full <- merge(df, lt_diff, all = T)
df.full <- merge(df.full, lt_cyc, all = T)
df.full <- merge(df.full, st_diff, all = T)
df.full <- merge(df.full, st_cyc, all = T)
colnames(df.full) <- c("tissue_sg", "total.cells", "lt_diff", "lt_cyc", "st_diff", "st_cyc")

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


##Split the first column into 2 columns
library(stringr)
df.full[c("tissue", "stage.group")] <- str_split_fixed(df.full$tissue_sg, "_", 2)
##Rearrange columns and remove original column names
#df.full <- df.full[c("tissue", "stage.group", "total.cells", "lt_diff", "lt_cyc", "st_diff", "st_cyc")]
df.prop <- df.full[c("tissue", "stage.group", "total.cells", "lt_diff_prop", "lt_cyc_prop", "st_cyc_prop", "st_diff_prop")]


##Convert the 3 column dataframe into a matrix for a dotplot
library(tidyr)
library(reshape2)
library(viridis)
library(tidyverse)
library(plotly)
library(purrr)

master.df <- df.prop %>% 
  mutate(tissue = recode_factor(tissue,
                                "axia" = "axial",
                                "basa" = "2",
                                "blas" = "3",
                                "mesoen" = "21",
                                "endo" = "4",
                                "eye" = "5",
                                "fin" = "6",
                                "gast" = "7",
                                "glia" = "8",
                                "hema" = "9",
                                "iono" = "10",
                                "mese" = "11",
                                "musc" = "12",
                                "neur" = "13",
                                "mura" = "14",
                                "otic" = "15",
                                "peri" = "16",
                                "pgc" = "17",
                                "pigm" = "18",
                                "pron" = "19",
                                "tast" = "20",
                                "unknown" = "22",
                                .ordered = T), 
         tissue = fct_rev(tissue),
         stage = fct_relevel(stage.group)) 
master.df <- df.prop

master.df <- master.df %>% mutate(stage = fct_relevel(stage), tissue = fct_rev(tissue))

master.df$stage.nice <- str_split_fixed(master.df$stage, "-", 2)
master.df <- master.df[c("tissue", "stage.nice", "total.cells", "lt_diff_prop", "lt_cyc_prop", "st_diff_prop", "st_cyc_prop")]
master.df$stage <- as.numeric(master.df$stage.nice[,1])
master.df <- master.df[c("tissue", "stage", "lt_cyc_prop", "lt_diff_prop", "st_cyc_prop", "st_diff_prop")]
master.df$tissue <- as.numeric(master.df$tissue)
master.df$stage <- as.numeric(master.df$stage)
master.df$lt_diff_prop <- as.numeric(master.df$lt_diff_prop)
master.df$lt_cyc_prop <- as.numeric(master.df$lt_cyc_prop)
master.df$st_diff_prop <- as.numeric(master.df$st_diff_prop)
master.df$st_cyc_prop <- as.numeric(master.df$st_cyc_prop)

library(ggplot2)
p <- ggplot(master.df, aes(x=stage, y=lt_cyc_prop, fill=tissue)) + 
  geom_area(alpha=0.2 , size=1, colour="black")

##Create a dataframe for each tissue and plot the short-term and long-term states
df.all <- as.data.frame(table(mama@meta.data[long.cycling.cells, "clust.sg"]))
colnames(df.all) <- c("clust.sg", "long_cyc")
df.all$tissue <- unlist(lapply(strsplit(x = as.character(df.all$clust.sg), split = "\\."), function(x) x[1]))
df.all$stage.group <- unlist(lapply(strsplit(x = as.character(df.all$clust.sg), split = "_"), function(x) x[2]))
df.all$clust <- unlist(lapply(strsplit(x = as.character(df.all$clust.sg), split = "_"), function(x) x[1]))
df.all$stage.nice <- str_split_fixed(df.all$stage, "-", 2)
df.all$stage <- as.numeric(df.all$stage.nice[,1])

df.axial <- df.all[df.all$tissue == "axia", ]
p <- ggplot(df.axial, aes(x=stage, y=long_cyc, fill=clust)) + 
  geom_area(alpha=0.2 , size=1, colour="black")

##Plot the scatter pie plots using cell proportions
library(scatterpie)
pdf("mama_cell_states_48h_scatter_pie_eps_5NN.pdf", height = 10, width = 12)
plot.pie <- ggplot() + geom_scatterpie(mapping = aes(x = stage, y = tissue), data = master.df, cols = c("lt_cyc_prop", "lt_diff_prop", "st_cyc_prop", "st_diff_prop"), pie_scale = 0.15, sorted_by_radius = FALSE, size = 0.1) + coord_fixed() + geom_dotplot(binwidth = 1) + labs(y = "Tissue", x = "Hour-Post-Fertilization")
plot(plot.pie)
dev.off()


##Now use a threshold of 36 hours - cells that have a stage difference more than 36 hours will be categorized as long-term states and ones which have stage difference less than 36 hours will be defined short term
##Define a slot in mama for long-term and short-term states
mama@meta.data$states.above.36 <- NA
mama@meta.data[which(mama@meta.data$stage.diff >= 36 & mama@meta.data$cc.status == "cycling"), "states.above.36"] <- "long-term_cycling"
mama@meta.data[which(mama@meta.data$stage.diff >= 36 & mama@meta.data$cc.status == "non-cycling"), "states.above.36"] <- "long-term_non-cycling"
mama@meta.data[which(mama@meta.data$stage.diff <= 36 & mama@meta.data$cc.status == "cycling"), "states.above.36"] <- "short-term_cycling"
mama@meta.data[which(mama@meta.data$stage.diff <= 36 & mama@meta.data$cc.status == "non-cycling"), "states.above.36"] <- "short-term_non-cycling"

colors.threshold <- setNames(c("#E9967A", "#800000", "#08519c", "#40E0D0"), c("long-term_cycling", "long-term_non-cycling", "short-term_cycling", "short-term_non-cycling"))
dpi <- 300
png("~/Box/zfext/global_analysis_curated/transient_vs_long-term/figures/mama_cell_states_thresh_36hr_eps_35_v3.png", height = 5*dpi, width = 5*dpi)
DimPlot(mama, group.by = "states.above.36", cols = colors.threshold, raster = F) + NoAxes() + NoLegend()
dev.off()

##Save the persistent and short-term states assignment for 36 hours as a threshold
saveRDS(mama@meta.data$states.above.36, "~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_cell_states_threshold_36_eps_40.rds")

##Plot a scatterpie plot for these cells that have stage difference above 36 hours
##Load mama metadata slot
states.above.36 <- readRDS("~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_cell_states_threshold_36_eps_40.rds")
mama@meta.data$long_v_short <- states.above.36
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


##Split the first column into 2 columns
library(stringr)
df.full[c("tissue", "num")] <- str_split_fixed(df.full$clust, ".", 1)
df.full$tissue <- unlist(lapply(strsplit(x = as.character(df.full$clust), split = "\\."), function(x) x[1]))
df.full$tissue.subsets <- cbind()
##Rearrange columns and remove original column names
#df.full <- df.full[c("tissue", "stage.group", "total.cells", "lt_diff", "lt_cyc", "st_diff", "st_cyc")]
df.prop <- df.full[c("tissue", "stage.group", "total.cells", "lt_diff_prop", "lt_cyc_prop", "st_cyc_prop", "st_diff_prop")]


##Convert the 3 column dataframe into a matrix for a dotplot
library(tidyr)
library(reshape2)
library(viridis)
library(tidyverse)
library(plotly)
library(purrr)

master.df <- df.prop %>% 
  mutate(tissue = recode_factor(tissue,
                                "axia" = "1",
                                "basa" = "2",
                                "blas" = "3",
                                "mesoen" = "21",
                                "endo" = "4",
                                "eye" = "5",
                                "fin" = "6",
                                "gast" = "7",
                                "glia" = "8",
                                "hema" = "9",
                                "iono" = "10",
                                "mese" = "11",
                                "musc" = "12",
                                "neur" = "13",
                                "mura" = "14",
                                "otic" = "15",
                                "peri" = "16",
                                "pgc" = "17",
                                "pigm" = "18",
                                "pron" = "19",
                                "tast" = "20",
                                "unknown" = "22",
                                .ordered = T), 
         tissue = fct_rev(tissue),
         stage = fct_relevel(stage.group)) 

master.df <- master.df %>% mutate(stage = fct_relevel(stage), tissue = fct_rev(tissue))

master.df$stage.nice <- str_split_fixed(master.df$stage, "-", 2)
master.df <- master.df[c("tissue", "stage.nice", "total.cells", "lt_diff_prop", "lt_cyc_prop", "st_diff_prop", "st_cyc_prop")]
master.df$stage <- as.numeric(master.df$stage.nice[,1])
master.df <- master.df[c("tissue", "stage", "lt_cyc_prop", "lt_diff_prop", "st_cyc_prop", "st_diff_prop")]
master.df$tissue <- as.numeric(master.df$tissue)
master.df$stage <- as.numeric(master.df$stage)
master.df$lt_diff_prop <- as.numeric(master.df$lt_diff_prop)
master.df$lt_cyc_prop <- as.numeric(master.df$lt_cyc_prop)
master.df$st_diff_prop <- as.numeric(master.df$st_diff_prop)
master.df$st_cyc_prop <- as.numeric(master.df$st_cyc_prop)

master.df <- master.df %>% mutate(stage = fct_relevel(stage), tissue = fct_rev(tissue))

##Plot the scatter pie plots using cell proportions
library(scatterpie)
pdf("mama_cell_states_36h_scatter_pie_eps_5NN.pdf", height = 10, width = 12)
plot.pie <- ggplot() + geom_scatterpie(mapping = aes(x = stage, y = tissue), data = master.df, cols = c("lt_cyc_prop", "lt_diff_prop", "st_cyc_prop", "st_diff_prop"), pie_scale = 0.15, sorted_by_radius = FALSE, size = 0.1) + coord_fixed() + geom_dotplot(binwidth = 1) + labs(y = "Tissue", x = "Hour-Post-Fertilization")
plot(plot.pie)
dev.off()


##Do the same with a threshold of 24 hours - cells that have a stage difference more than 24 hours will be categorized as long-term states and ones which have stage difference less than 24 hours will be defined short term
##Define a slot in mama for long-term and short-term states
mama@meta.data$states.above.24 <- NA
mama@meta.data[which(mama@meta.data$stage.diff >= 24 & mama@meta.data$cc.status == "cycling"), "states.above.24"] <- "long-term_cycling"
mama@meta.data[which(mama@meta.data$stage.diff >= 24 & mama@meta.data$cc.status == "non-cycling"), "states.above.24"] <- "long-term_non-cycling"
mama@meta.data[which(mama@meta.data$stage.diff <= 24 & mama@meta.data$cc.status == "cycling"), "states.above.24"] <- "short-term_cycling"
mama@meta.data[which(mama@meta.data$stage.diff <= 24 & mama@meta.data$cc.status == "non-cycling"), "states.above.24"] <- "short-term_non-cycling"

colors.threshold <- setNames(c("#E9967A", "#800000", "#08519c", "#40E0D0"), c("long-term_cycling", "long-term_non-cycling", "short-term_cycling", "short-term_non-cycling"))
dpi <- 300
png("~/Box/zfext/global_analysis_curated/transient_vs_long-term/figures/mama_cell_states_thresh_24hr_eps_40.png", height = 5*dpi, width = 5*dpi)
DimPlot(mama, group.by = "states.above.24", cols = colors.threshold, raster = F) + NoAxes() + NoLegend()
dev.off()

##Save the persistent and short-term states assignment for 24 hours as a threshold
saveRDS(mama@meta.data$states.above.24, "~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_cell_states_threshold_24_eps_40.rds")

##Plot a scatterpie plot for these cells that have stage difference above 24 hours
##Load mama metadata slot
states.above.24 <- readRDS("~/Box/zfext/global_analysis_curated/transient_vs_long-term/mama_cell_states_threshold_24_eps_40.rds")
mama@meta.data$long_v_short <- states.above.24
DimPlot(mama, group.by = "states.above.24", raster = F) + NoAxes() + NoLegend()

long.cycling.cells <- rownames(mama@meta.data)[mama@meta.data$states.above.48 == "long-term_cycling"]
long.non_cycling.cells <- rownames(mama@meta.data)[mama@meta.data$states.above.48 == "long-term_non-cycling"]
short.cycling.cells <- rownames(mama@meta.data)[mama@meta.data$states.above.48 == "short-term_cycling"]
short.non_cycling.cells <- rownames(mama@meta.data)[mama@meta.data$states.above.48 == "short-term_non-cycling"]

df <- as.data.frame(table(mama@meta.data$tissue.sg))
colnames(df) <- c("tissue_sg", "total.cells")
lt_diff <- as.data.frame(table(mama@meta.data[long.non_cycling.cells, "tissue.sg"]))
lt_cyc <- as.data.frame(table(mama@meta.data[long.cycling.cells, "tissue.sg"]))
st_diff <- as.data.frame(table(mama@meta.data[short.non_cycling.cells, "tissue.sg"]))
st_cyc <- as.data.frame(table(mama@meta.data[short.cycling.cells, "tissue.sg"]))

colnames(lt_diff) <- c("tissue_sg", "lt_diff")
colnames(lt_cyc) <- c("tissue_sg", "lt_cyc")
colnames(st_diff) <- c("tissue_sg", "st_diff")
colnames(st_cyc) <- c("tissue_sg", "st_cyc")

df.full <- merge(df, lt_diff, all = T)
df.full <- merge(df.full, lt_cyc, all = T)
df.full <- merge(df.full, st_diff, all = T)
df.full <- merge(df.full, st_cyc, all = T)
colnames(df.full) <- c("tissue_sg", "total.cells", "lt_diff", "lt_cyc", "st_diff", "st_cyc")

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


##Split the first column into 2 columns
library(stringr)
df.full[c("tissue", "stage.group")] <- str_split_fixed(df.full$tissue_sg, "_", 2)
##Rearrange columns and remove original column names
#df.full <- df.full[c("tissue", "stage.group", "total.cells", "lt_diff", "lt_cyc", "st_diff", "st_cyc")]
df.prop <- df.full[c("tissue", "stage.group", "total.cells", "lt_diff_prop", "lt_cyc_prop", "st_cyc_prop", "st_diff_prop")]


##Convert the 3 column dataframe into a matrix for a dotplot
library(tidyr)
library(reshape2)
library(viridis)
library(tidyverse)
library(plotly)
library(purrr)

master.df <- df.prop %>% 
  mutate(tissue = recode_factor(tissue,
                                "axia" = "1",
                                "basa" = "2",
                                "blas" = "3",
                                "mesoen" = "21",
                                "endo" = "4",
                                "eye" = "5",
                                "fin" = "6",
                                "gast" = "7",
                                "glia" = "8",
                                "hema" = "9",
                                "iono" = "10",
                                "mese" = "11",
                                "musc" = "12",
                                "neur" = "13",
                                "mura" = "14",
                                "otic" = "15",
                                "peri" = "16",
                                "pgc" = "17",
                                "pigm" = "18",
                                "pron" = "19",
                                "tast" = "20",
                                "unknown" = "22",
                                .ordered = T), 
         tissue = fct_rev(tissue),
         stage = fct_relevel(stage.group)) 

master.df <- master.df %>% mutate(stage = fct_relevel(stage), tissue = fct_rev(tissue))

master.df$stage.nice <- str_split_fixed(master.df$stage, "-", 2)
master.df <- master.df[c("tissue", "stage.nice", "total.cells", "lt_diff_prop", "lt_cyc_prop", "st_diff_prop", "st_cyc_prop")]
master.df$stage <- as.numeric(master.df$stage.nice[,1])
master.df <- master.df[c("tissue", "stage", "lt_cyc_prop", "lt_diff_prop", "st_cyc_prop", "st_diff_prop")]
master.df$tissue <- as.numeric(master.df$tissue)
master.df$stage <- as.numeric(master.df$stage)
master.df$lt_diff_prop <- as.numeric(master.df$lt_diff_prop)
master.df$lt_cyc_prop <- as.numeric(master.df$lt_cyc_prop)
master.df$st_diff_prop <- as.numeric(master.df$st_diff_prop)
master.df$st_cyc_prop <- as.numeric(master.df$st_cyc_prop)

##Plot the scatter pie plots using cell proportions
library(scatterpie)
pdf("mama_cell_states_24h_scatter_pie_eps_5NN.pdf", height = 10, width = 12)
plot.pie <- ggplot() + geom_scatterpie(mapping = aes(x = stage, y = tissue), data = master.df, cols = c("lt_cyc_prop", "lt_diff_prop", "st_cyc_prop", "st_diff_prop"), pie_scale = 0.15, sorted_by_radius = FALSE, size = 0.1) + coord_fixed() + geom_dotplot(binwidth = 1) + labs(y = "Tissue", x = "Hour-Post-Fertilization")
plot(plot.pie)
dev.off()



