library(Seurat)
library(URD)

#Grabbing smooth muscle cells from different sources:
#Smooth muscle cell markers: acta2/tagln, myl9a, lmod1b, cnn1b, desmb etc. 
#Let's see where these genes are expressed?

plotDim(mama, "foxc1b")
plotDim(mama, "iter_cluster", label.clusters = T, legend = F)
plotDimHighlight(mama, "iter_cluster", "6.1", legend = F)

##Read in new mural cell object
obj.smc <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/mural/obj_seurat/mural_seurat.rds")
markers <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/mural/obj_seurat/mural_markers.rds")

stage.colors.new <- c(
  colorRampPalette(c("#BA55D3", "#C71585", "#FF1493"))(7), # 14, 16, 18, 21`
  rep(RColorBrewer::brewer.pal(6, "Oranges"), each = 2), # 24, 26 / 28, 30 / 32, 34 / 36, 38 / 40, 42 / 44, 46
  rep(RColorBrewer::brewer.pal(6, "Greens"), each = 2), # 48, 50 / 52, 54 / 56, 58 / 60, 62 / 64, 66 / 68, 70
  rep(RColorBrewer::brewer.pal(6, "Blues"), each = 2),  # 72, 74 / 76, 78 / 80, 82 / 84, 86 / 88, 90 / 92, 94
  rep(RColorBrewer::brewer.pal(7, "Purples"), each = 2) # 96, 98 / 100, 102 / 104, 106 / 108, 110 / 112, 114 / 116, 118 / 120
)

dpi <- 300
png(file = "mural_cells_UMAP_stage_colored.png", width = 5 * dpi, height = 5 * dpi)
DimPlot(obj.smc, group.by = "stage.nice", cols = stage.colors.new, pt.size = 2)
dev.off()

DimPlot(obj.smc, label = T)
DimPlot(obj.smc, group.by = "stage.group")
DimPlot(obj.smc, group.by = "RNA_snn_res.3")

##Ok, RNA_snn_res.3 was a better representation of the clustering of these cells
Idents(obj.smc) <- obj.smc@meta.data$RNA_snn_res.3

posMarkers.wilcox <- sapply(levels(Idents(obj)), function(i) FindMarkers(obj,i,assay.type='RNA',test.use="wilcox",only.pos=T, min.cells.group = 1),simplify=F)
posMarkers.roc <- sapply(levels(Idents(obj)), function(i) FindMarkers(obj,i,test.use="roc",only.pos=T),simplify=F)
markers <- list(wilcox = posMarkers.wilcox, roc = posMarkers.roc)
saveRDS(markers, file="~/Box/zfext/annotations_celltype_curated_newMama/mural_cells/obj_seurat/mural-cells_markers.rds")
saveRDS(obj, file="~/Box/zfext/annotations_celltype_curated_newMama/mural_cells/obj_seurat/mural-cells_seurat.rds")

## I. How many pericyte populations do we detect? Can we distinguish them transcriptionally?
##Ok, there are several pieces in the atlas that express the pericyte marker ndufa4l2a. So get those cells out and use them as a separate cluster.
##Ok, break the cluster 7 into the two pieces based on the two branches of the URD trajectory
Idents(obj.smc) <- obj.smc@meta.data$RNA_snn_res.2
##cells.peri.1 loaded from the previously saved workspace includes the pericyte signature cells originating from the NCC-derived transient state
##Load a second object with idents set as RNA_snn_res.3, but save the new cluserting on obj.sm
obj.sm <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/mural_cells/obj_seurat/mural-cells_seurat.rds")
cells.9 <- WhichCells(obj.sm, idents = c("9"))
cells.peri.new <- setdiff(cells.peri.1, cells.9)
##Check where the new cluster is defined
DimPlot(obj.smc, cells.highlight = cells.ndu)
cells.ndu <- WhichCells(obj.sm, idents = c("6"), expression = ndufa4l2a > 4)
##Set the respective idents for the new pericyte populations
obj.sm <- SetIdent(obj.sm, cells = cells.ndu, value = "27")
obj.sm <- SetIdent(obj.sm, cells = cells.art.short, value = "28")
Idents(obj.smc) <- obj.smc@meta.data$RNA_snn_res.3
obj.sm <- SetIdent(obj.sm, cells = cells.art.short, value = "28")
obj.sm <- SetIdent(obj.sm, cells = WhichCells(obj.smc, idents = c("22")), value = "vSMC_artery_1")

cells.artery.other <- setdiff(cells.art.long, WhichCells(obj.sm, idents = c("22")))
DimPlot(obj.sm, cells.highlight = cells.artery.other)

obj.sm <- obj.smc
DimPlot(obj.sm, label = T, pt.size = 1.2)
FeaturePlot(obj.smc, features = c("ndufa4l2a", "il13ra2", "tesca", "abcc9"))
FeaturePlot(obj.smc, features = c("bgna", "cnn1a", "cln6b", "tdrd9"))
FeaturePlot(obj.smc, features = c("bgna", "ndufa4l2a", "kcnk18", "smpx"))

cluster.names.new <- c("vSMC-aorta", "pericyte_3", "pericyte_2", "mesenchyme_int", "cycling", "neural", "cycling", "i-SMCs", "sb-SMCs", "vSMC-transient", "pericyte-transient", "cycling", "vi-SMC", "epicardium", "vSMC-OFT", "HSCs", "cardiac_muscle", "myofibroblasts", "myocardium", "neural", "mesenchyme_int", "HSCs", "pericytes-brain", "vSMC-coro", "neural", "SMC-post", "NCC-precursors", "fibroblasts", "epidermis")
names(cluster.names.new) <- levels(obj.smc)
obj.sm <- RenameIdents(obj.smc, cluster.names.new) 

cells.smc.transient <- WhichCells(obj.smc, idents = c("9"), expression = acta2 > 3.5)
DimPlot(obj.smc, cells.highlight = cells.smc.transient)
obj.sm <- SetIdent(obj.sm, cells = cells.smc.transient, value = "vSMC-transient")

cells.20 <- WhichCells(obj.smc, idents = c("20"))
obj.sm <- SetIdent(obj.sm, cells = cells.20, value = "cycling")

DimPlot(obj.sm, label = T)
saveRDS(obj.sm, "~/Box/zfext/annotations_celltype_curated_newMama/mural_cells/obj_seurat/mural-cells_seurat_reclustered.rds")
obj.sm <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/mural_cells/obj_seurat/mural-cells_seurat_reclustered.rds")

## II. Ok, is there a fourth population of pericytes that is at the tip of cluster 21 of the newly reclustered mural-cell atlas?
##In that case, let's take those cells out and compare them to the other pericyte populations
##How distinct are these? I know that these also express acta2 along with ndufa4l2a which is a cell state that has never been observed before. 
## Generally, ndufa4l2a is known to be a pericyte-specific markers while acta2 is known to be a SMC-specific marker

##Grab the ndufa4l2a expressing cells from the cluster 21
cells.ndu <- WhichCells(obj.smc, idents = c("21", "11"), expression = ndufa4l2a > 1)   ##14 cells
cells.smc <- setdiff(WhichCells(obj.smc, idents = c("21", "11")), cells.ndu)
cells.smc.rest <- WhichCells(obj.smc, idents = c("2", "12", "15"))
cells.smc.all <- unlist(unique(list(cells.smc, cells.smc.rest)))
##Create a sub-object with all pericyte populations across the atlas
##Are there three pericyte populations?
##Let's compare the three pericyte populations and the transient state - that's becoz the early pericyte markers are also not known
cells.peri.1 <- WhichCells(obj.smc, idents = c("9"), expression = ndufa4l2a > 1)
cells.peri.2 <- WhichCells(obj.smc, idents = c("20"))
cells.peri.3 <- WhichCells(obj.smc, idents = c("4"))
cells.peri.4 <- WhichCells(obj.smc, idents = c("3"))

cells.all <- unlist(unique(list(cells.peri.1, cells.peri.2, cells.peri.3, cells.peri.4)))

obj.peri <- subset(obj.smc, cells = cells.all)

obj.peri@meta.data$pericyte.clusters <- NA
obj.peri@meta.data[cells.peri.1, "pericyte.clusters"] <- "pericyte-transient"
obj.peri@meta.data[cells.peri.2, "pericyte.clusters"] <- "pericyte_0"
obj.peri@meta.data[cells.peri.3, "pericyte.clusters"] <- "pericytes_1"
obj.peri@meta.data[cells.peri.4, "pericyte.clusters"] <- "pericytes_2"
obj.peri@meta.data[cells.ndu, "pericyte.clusters"] <- "mFB"

Idents(obj.peri) <- obj.peri@meta.data$pericyte.clusters

library(dplyr)
markers.peri <- FindAllMarkers(obj.peri, assay = "RNA", logfc.threshold = 0.25, only.pos = T, min.pct = 0.25)
top20 <- markers.peri %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
gene.list <- unlist(unique(list(top20$gene, "pdgfrb", "ndufa4l2a", "abcc9", "kcnj8", "plp1b", "bgna", "pld6", "cnn1a", "pthlha", "fbp2", "tdrd9", "cln6b")))

filter.genes <- function(genes) {
  mt.genes <- grep("^mt-", ignore.case = T, genes, value = T)
  many.genes <- grep("\\(1 of many\\)", ignore.case = T, genes, value = T)
  cox.genes <- grep("^cox", ignore.case = T, genes, value = T)
  rp.genes <- grep("^rp", ignore.case = T, genes, value = T)
  hb.genes <- grep("^hb", ignore.case = T, genes, value = T)
  si.genes <- grep("^si:", ignore.case = T, genes, value = T)
  zgc.genes <- grep("^zgc:", ignore.case = T, genes, value = T)
  loc.genes <- grep("^XLOC", ignore.case = T, genes, value = T)
  bx.genes <- grep("^BX", ignore.case = T, genes, value = T)
  cr.genes <- grep("^CT", ignore.case = T, genes, value = T)
  return(setdiff(genes, c(mt.genes, many.genes, cox.genes, rp.genes, hb.genes, si.genes, zgc.genes, loc.genes, cr.genes, bx.genes)))
}

gene.list <- top20$gene
gene.list <- filter.genes(gene.list)

library(viridis)
pdf("pericyte_comparisons_dotplot_5.pdf", width = 12, height = 28)
DotPlot(obj.peri, features = unique(gene.list), dot.scale = 10, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

markers.peri.trans <- markers.peri$gene[which(markers.peri$cluster == "pericyte-transient" & markers.peri$pct.2 < 0.08)]
markers.peri.1 <- markers.peri$gene[which(markers.peri$cluster == "pericytes_0" & markers.peri$pct.2 < 0.08)]
markers.peri.2 <- markers.peri$gene[which(markers.peri$cluster == "pericytes_1" & markers.peri$pct.2 < 0.08)]
markers.peri.3 <- markers.peri$gene[which(markers.peri$cluster == "pericytes_2" & markers.peri$pct.2 < 0.08)]

gene.list <- unlist(unique(list(markers.peri.trans, markers.peri.1, markers.peri.2, markers.peri.3)))
gene.list.ordered <- factor(gene.list, levels = unique(c(markers.peri.trans, markers.peri.1, markers.peri.2, markers.peri.3, "acta2", "tagln", "myl9a")))
pdf("pericyte_comparisons_dotplot_5.pdf", width = 12, height = 28)
DotPlot(obj.peri, features = unique(gene.list.ordered), dot.scale = 10, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()



library(BioVenn)
venn <- draw.venn(markers.peri.1, markers.peri.2, markers.peri.3)

FeaturePlot(obj.smc, c("nkx2.3", "pgfb", "nptna", "rgs4"))


##Make a small dotplot highlighting some of the key differentially expressed markers between the pericyte populatios
gene.list <- c("ndufa4l2a", "abcc9", "rasgef1ba", "pdgfra", "pdgfrb", "ctgfa", "bgna", "dcn", "nt5c1aa", "esama", "epas1a", "adma", "olfml3b", "lmx1bb", "nptna", "sytl4", "pgfb", "acta2", "tagln", "myl9a")
pdf("pericyte_comparisons_dotplot_small_v3.pdf", width = 8, height = 10)
DotPlot(obj.peri, features = gene.list, idents = c("pericyte-transient", "pericyte_0", "pericytes_1", "pericytes_2", "vSMCs"), dot.scale = 10, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

##Conclusion: Differential expression of these populations revealed that these cells in cluster 21 differentially expressed ndufa4l2a from the rest of the SMCs in the same branch, but expressed acta2, tagln and all other markers that are expressed in vSMCs. Now as these cells doesn't have as many pericyte markers that are expressed here, maybe these are a SMC-->pericyte transient state (but still very much of a vSMC fate)

## III. Immature pericytes and vSMC markers: So far, immature and mature pericytes cannot be distinguished due to the lack of markers. When researchers have observed pericytes, it has not been clear where they came from and whether they are already mature

FeaturePlot(obj.smc, c("sytl4", "epas1a", "ebf1b", "pgfb")) ## expressed in the pericyte transient state overlapping ndufa4l2a and abcc9 expression and turns off in other pericyte clusters

##Similarly what maybe some markers of immature vSMCs?
FeaturePlot(obj.smc, c("foxc1a", "mef2cb", "slc7a2", "adrb3a")) ## slc7a2 and mef2cb seems to be good markers of this vSMC-transient state. 

##Check what stages these pericyte and SMc transient states belong to
cells.smc.transient <- WhichCells(obj.smc, idents = c("12"))
cells.peri.transient <- WhichCells(obj.smc, idents = c("9"))
table(obj.smc@meta.data[cells.peri.transient, "stage.nice"])
##Both can be found at 5 dpf. 

## IV. Now that we detected four different pericyte populations, how are these pericyte populations different? 
## I mean I already detected a bunch of genes different between them and plotted a dotplot, but what does that actually mean?
## Hence, I need some type of a Gene enrichment analysis - GSEA, KEGG, GO, Wikipathways?? What is the best way? 

## Maybe try using Msigdb2 first through clusterProfiler
library(clusterProfiler)
##Do pathway enrichment using MSigDB
install.packages("msigdbr")
#Load package
library(msigdbr)
##Load all gene sets for zebrafish
all_gene_sets = msigdbr(species = "zebrafish")
head(all_gene_sets)

msigdbr_species()

##We can retrieve specific genesets as well
h_gene_sets = msigdbr(species = "zebrafish", category = "H")
head(h_gene_sets)

cgp_gene_sets = msigdbr(species = "zebrafish", category = "C2", subcategory = "CGP")
head(cgp_gene_sets)

##Helper function to show the different collections
msigdbr_collections()

##The msigdbr() function output is a data frame and can be manipulated using more standard methods
all_gene_sets %>%
  dplyr::filter(gs_cat == "H") %>% 
  head()

##Let's look at the GO/KEGG enrichment of the pericyte genes for the 2 pericyte populations
##See R-script - mural_cells_GEA.R
genes.peri.1 <- markers.peri$gene[markers.peri$cluster == "pericyte-transient"]
genes.peri.2 <- markers.peri$gene[markers.peri$cluster == "pericytes_brain"]
genes.peri.3 <- markers.peri$gene[markers.peri$cluster == "pericytes_2"]
genes.peri.4 <- markers.peri$gene[markers.peri$cluster == "pericytes_3"]
genes.peri.5 <- markers.peri$gene[markers.peri$cluster == "pericytes_4"]
genes.peri.6 <- markers.peri$gene[markers.peri$cluster == "vSMCs"]

genes.peri.1 <- filter.genes(genes.peri.1)

##The msigdbr output can be used with various popular pathway analysis packages
msigdbr_t2s = all_gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
msigdbr_df_2 <- as.data.frame(msigdbr_t2s)
x.peri.1 <- enricher(gene = genes.peri.1, TERM2GENE = msigdbr_t2s)
x.peri.2 <- enricher(gene = genes.peri.2, TERM2GENE = msigdbr_t2s)
x.peri.3 <- enricher(gene = genes.peri.3, TERM2GENE = msigdbr_t2s)
x.peri.4 <- enricher(gene = genes.peri.4, TERM2GENE = msigdbr_t2s)
x.peri.5 <- enricher(gene = genes.peri.5, TERM2GENE = msigdbr_t2s)
x.peri.6 <- enricher(gene = genes.peri.6, TERM2GENE = msigdbr_t2s)
y <- GSEA(entrez_ids_peri.1$ENTREZID, TERM2GENE = msigdbr_t2s)

msigdbr_t2g = all_gene_sets %>% dplyr::distinct(gs_name, entrez_gene, gene_symbol) %>% as.data.frame()
msigdbr_df <- as.data.frame(msigdbr_t2g)
xx <- enricher(gene = gene.list, TERM2GENE = msigdbr_t2s)
y <- GSEA(gene.list, TERM2GENE = msigdbr_t2s)

gene.list.entrez.ids <- msigdbr_df$entrez_gene[which(msigdbr_df$gene_symbol == gene.list)]

BiocManager::install("enrichplot")
library(enrichplot)
library(ggnewscale)

p1 <- cnetplot(xx, foldChange=genes.peri.1)
p2 <- cnetplot(xx, categorySize="pvalue", foldChange=gene.list)
p3 <- cnetplot(xx, foldChange=gene.list, circular = TRUE, colorEdge = TRUE) 

barplot(pt, showCategory=20)

p1 <- cnetplot(xx, node_label="category", 
               cex_label_category = 0.2) 
p2 <- cnetplot(xx, node_label="gene", 
               cex_label_gene = 0.8) 
p3 <- cnetplot(xx, node_label="all") 
p4 <- cnetplot(xx, node_label="none", 
               color_category='firebrick', 
               color_gene='steelblue') 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

msigdbr_t2g = all_gene_sets %>% dplyr::distinct(gene_symbol, entrez_gene) %>% as.data.frame()


msigdbr_t2g = all_gene_sets %>% dplyr::distinct(gs_name, entrez_gene, gene_symbol) %>% as.data.frame()
msigdbr_df <- as.data.frame(msigdbr_t2g)
#Pericyte-transient
pt <- enricher(gene = genes.peri.1, TERM2GENE = msigdbr_t2g)
pt1 <- cnetplot(pt, showCategory = 10, foldChange=genes.peri.1, cex_label_category = 0.2)
pt2 <- cnetplot(pt, showCategory = 10, categorySize="pvalue", foldChange=genes.peri.1, cex_label_category = 0.5)
pt3 <- cnetplot(pt, showCategory = 10, foldChange=genes.peri.1, circular = TRUE, colorEdge = TRUE, cex_label_category = 0.5) 
cowplot::plot_grid(pt1, pt2, pt3, ncol=3, labels=LETTERS[1:3])

dpi <- 300
png("pericyte_brain_msigDB.png", width = 5*dpi, height = 5*dpi)
cnetplot(pt, showCategory = 10, foldChange=genes.peri.1, circular = TRUE, colorEdge = TRUE, cex_label_category = 0.5, node_label = "gene") + NoLegend()
cnetplot(pb, showCategory = 10, foldChange=genes.peri.2, circular = TRUE, colorEdge = TRUE, cex_label_category = 0.5, node_label = "gene") + NoLegend()
dev.off()

#Pericyte-brain
pb <- enricher(gene = genes.peri.2, TERM2GENE = msigdbr_t2g)
pb1 <- cnetplot(pb, foldChange=genes.peri.2)
pb2 <- cnetplot(pb, categorySize="pvalue", foldChange=genes.peri.2)
pb3 <- cnetplot(pb, showCategory = 10, foldChange=genes.peri.2, circular = TRUE, colorEdge = TRUE, cex_label_category = 0.5) 
cowplot::plot_grid(pb1, pb2, pb3, ncol=3, labels=LETTERS[1:3])

msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
gsva(gset.idx.list = msigdbr_list)

kegg_organism <- "dre"
kk2 <- gseKEGG(geneList = entrez_ids, 
               organism = kegg_organism,
               nPerm = 10000,
               minGSSize = 3,
               maxGSSize = 800,
               pvalueCutoff = 0.5,
               pAdjustMethod = "none",
               keyType = "ncbi-geneid")

##Reactome PA
##Enrichment analysis is a widely used approach to identify biological themes. ReactomePA implemented enrichPathway() that uses hypergeometric model to assess whether the number of selected genes associated with a reactome pathway is larger than expected.
BiocManager::install("ReactomePA")
library(ReactomePA)
entrez_ids_peri.1 <- clusterProfiler::bitr(geneID = genes.peri.1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Dr.eg.db")
entrez_ids_peri.2 <- clusterProfiler::bitr(geneID = genes.peri.2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Dr.eg.db")
head(entrez_ids_peri.1)
x <- enrichPathway(entrez_ids_peri.1$ENTREZID, organism = "zebrafish", pvalueCutoff = 0.05, readable = T)
de <- sort(entrez_ids_peri.1$ENTREZID, decreasing = T)
y <- gsePathway(de, organism = "zebrafish", 
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                verbose = TRUE)

##Differentiating between intestinal and SB smooth muscle cells
obj.int <- subset(obj.smc, idents = c("sb-SMCs", "i-SMCs"))
markers.int <- FindAllMarkers(obj.int, assay = "RNA", logfc.threshold = 0.25, only.pos = T, min.pct = 0.25, test.use = "wilcox")
top30 <- markers.int %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

p <- DotPlot(obj.int, features = top30$gene) + coord_flip()
p

tf.list <- read.delim(file="~/Box/zfext/02-Clustering/2021-03 Iterative Clustering/gene_info/tfs/2021-06-25_zebrafish_LTA_TFs.txt", header = T, sep = "")
markers.int <- rownames(markers$wilcox$`6`)
markers.tfs.only <- markers.int[which(markers.int %in% tf.list$Symbol), ]

##Highlight the four pericyte populations
obj.smc@meta.data$pericyte.clusters <- NA
obj.smc@meta.data[cells.peri.1, "pericyte.clusters"] <- "pericyte-transient"
obj.smc@meta.data[cells.peri.2, "pericyte.clusters"] <- "pericyte-brain"
obj.smc@meta.data[cells.peri.3, "pericyte.clusters"] <- "pericyte-2"
obj.smc@meta.data[cells.peri.4, "pericyte.clusters"] <- "pericyte-3"

DimPlot(obj.smc, group.by = "pericyte.clusters", na.value = "#F5F5F5")

##Now grab the five distinct SMC populations and compare between them
##Ok now perform a differential expression between all vascular SMC subtypes. 
##
obj.sm <- obj.smc
Idents(obj.sm) <- obj.sm@meta.data$annotations
obj.smooth <- subset(obj.sm, idents = c("vSMC_aorta", "vSMC-OFT", "vSMC-transient", "vSMC-artery_1", "vSMC-artery_2"))
obj.smooth <- subset(obj.smc, idents = c("2", "12", "15", "11", "21"))
markers.smc <- FindAllMarkers(obj.smooth, assay = "RNA", logfc.threshold = 0.25, only.pos = T, min.pct = 0.25)
library(dplyr)
top30 <- markers.smc %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
gene.list <- unique(top30$gene)
gene.list <- filter.genes(gene.list)

pdf("vSMC_comparison_dotplot_3.pdf", width = 10, height = 8)
gene.list <- c("slc2a1b", "slc16a3", "cxcl14", "htra1a", "elnb", "loxa", "ifitm1", "ccl25b", "mcamb", "mylkb", "grem2b", "aldh1a2", "sfrp2", "nkx3.2", "cpn1", "acta2", "tagln", "myl9a")
DotPlot(obj.smooth, features = gene.list, dot.scale = 10) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

Idents(obj.smc) <- obj.smc@meta.data$RNA_snn_res.3

##Now compare between vascular SMCs, visceral SMCs and pericytes
cells.vascular <- WhichCells(obj.sm, idents = c("vSMC-artery_1", "vSMC-artery_2", "vSMC-OFT", "vSMC_aorta"))
cells.visceral <- WhichCells(obj.smc, idents = c("7", "10", "8"))
cells.pericytes <- WhichCells(obj.smc, idents = c("20", "4"))

obj.sm <- subset(obj.smc, cells = unlist(unique(list(cells.vascular, cells.visceral, cells.pericytes))))

obj.sm@meta.data$smc.class <- NA
obj.sm@meta.data[cells.vascular, "smc.class"] <- "vascular"
obj.sm@meta.data[cells.visceral, "smc.class"] <- "visceral"
obj.sm@meta.data[cells.pericytes, "smc.class"] <- "pericytes"

colors <- setNames(c("#71B000", "#FF6C91", "#00B7E8"), c("vascular", "visceral", "pericytes"))
DimPlot(obj.sm, group.by = "smc.class", cols = colors)

Idents(obj.sm) <- obj.sm@meta.data$smc.class
my_levels <- c("vSMC-transient", "vSMC_aorta", "vSMC-OFT", "vSMC-artery_1", "vSMC-artery_2", "vi-SMCs", "i-SMCs", "sb-SMCs", "pericyte_transient", "pericytes-brain", "pericyte_2", "pericyte_3", "NCC-precursors", "mesenchyme_int", "cycling", "unknown", "cardiac_muscle", "myocardium", "epicardium", "HSCs", "fibroblasts", "myofibroblasts?", "SMC-post")
# Relevel object@ident
Idents(obj.sm) <- factor(x = Idents(obj.sm), levels = my_levels)

Idents(obj.sm) <- obj.sm@meta.data$smc.class

genes.smc.peri <- FindAllMarkers(obj.sm, assay = "RNA", only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
library(dplyr)
top15 <- genes.smc.peri %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
gene.list <- unlist(top20$gene)

filter.genes <- function(genes) {
  mt.genes <- grep("^mt-", ignore.case = T, genes, value = T)
  many.genes <- grep("\\(1 of many\\)", ignore.case = T, genes, value = T)
  cox.genes <- grep("^cox", ignore.case = T, genes, value = T)
  rp.genes <- grep("^rp", ignore.case = T, genes, value = T)
  hb.genes <- grep("^hb", ignore.case = T, genes, value = T)
  si.genes <- grep("^si:", ignore.case = T, genes, value = T)
  zgc.genes <- grep("^zgc:", ignore.case = T, genes, value = T)
  loc.genes <- grep("^LOC", ignore.case = T, genes, value = T)
  bx.genes <- grep("^BX", ignore.case = T, genes, value = T)
  cr.genes <- grep("^CT", ignore.case = T, genes, value = T)
  return(setdiff(genes, c(mt.genes, many.genes, cox.genes, rp.genes, hb.genes, si.genes, zgc.genes, loc.genes, cr.genes, bx.genes)))
}

gene.list <- filter.genes(gene.list)
gene.list <- c("wfdc2", "fbln5", "fbp2", "ldha", "slc16a3", "loxa", "desmb", "cald1b", "csrp1b", "nkx2.3", "smtna", "smtnb", "abcc9", "ndufa4l2a", "pdgfrb", "bgnb", "rgs5b", "foxl1")

pdf("pericyte_vaSMC_viSMC_dotplot.pdf", width = 8, height = 10)
DotPlot(obj.sm, features = gene.list, dot.scale = 10, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()



# Relevel object@ident
object@ident <- factor(x = object@ident, levels = my_levels)
VlnPlot(obj.sm, genes.peri, idents = c("vSMC-transient", "vSMC_aorta", "vSMC-OFT", "vSMC-artery_1", "vSMC-artery_2", "vi-SMCs", "i-SMCs", "sb-SMCs", "pericyte_transient", "pericytes-brain", "pericyte_2", "pericyte_3"), stack = TRUE, flip=TRUE, combine = FALSE) + NoLegend()

dpi <- 300
png("SMC_pericyte_vlnplot_noAxes_v2.png", height = 5*dpi, width = 5*dpi)
VlnPlot(obj.sm, genes.peri, idents = c("vSMC-transient", "vSMC_aorta", "vSMC-OFT", "vSMC-artery_1", "vSMC-artery_2", "vi-SMCs", "i-SMCs", "sb-SMCs", "pericyte_transient", "pericytes-brain", "pericyte_2", "pericyte_3"), stack = TRUE, flip=TRUE, combine = FALSE) + NoLegend()
dev.off()

dpi <- 300
png("UMAP_with_SMC_pericyte.png", width = 5*dpi, height = 5*dpi)
DimPlot(obj.sm, group.by = "smc.class", cols = colors) + NoLegend() + NoAxes()
dev.off()

Idents(obj.smc) <- obj.sm@meta.data$smc.class
obj.peri <- subset(obj.smc, idents = c("9", "20", "4", "3"))
genes.peri <- FindAllMarkers(obj.peri, assay = "RNA", only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
library(dplyr)
top10 <- genes.peri %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
gene.list <- unlist(top10$gene)

filter.genes <- function(genes) {
  mt.genes <- grep("^mt-", ignore.case = T, genes, value = T)
  many.genes <- grep("\\(1 of many\\)", ignore.case = T, genes, value = T)
  cox.genes <- grep("^cox", ignore.case = T, genes, value = T)
  rp.genes <- grep("^rp", ignore.case = T, genes, value = T)
  hb.genes <- grep("^hb", ignore.case = T, genes, value = T)
  si.genes <- grep("^si:", ignore.case = T, genes, value = T)
  zgc.genes <- grep("^zgc:", ignore.case = T, genes, value = T)
  loc.genes <- grep("^LOC", ignore.case = T, genes, value = T)
  bx.genes <- grep("^BX", ignore.case = T, genes, value = T)
  cr.genes <- grep("^CT", ignore.case = T, genes, value = T)
  return(setdiff(genes, c(mt.genes, many.genes, cox.genes, rp.genes, hb.genes, si.genes, zgc.genes, loc.genes, cr.genes, bx.genes)))
}

gene.list <- filter.genes(gene.list)
gene.list <- c("ggt5a", "pdgfab", "pla1a", "myocd", "pdlim3a", "pdlim3b", "rprmb", "nptna", "vcla", "nkx3.3", "nkx2.3", "smtna", "smtnb", "desmb", "pdlim1", "foxc1b", "vim", "fbln5", "wfdc2", "cd248a", "loxl5b", "loxa", "cnn1b", "tagln", "acta2", "ndufa4l2a", "abcc9", "pdgfrb", "rasl12", "kcnj8", "kcne4", "notch2", "notch3", "bgnb", "tm4sf18", "pros1", "c1qtnf5", "rasgef1ba")
library(viridis)
pdf("pericyte_vs_SMC_dotplot.pdf", width = 10, height = 16)
DotPlot(obj.peri, features = gene.list, dot.scale = 10) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

obj.vi <- subset(obj.smc, idents = c("8", "10", "13"))
genes.vi <- FindAllMarkers(obj.vi, assay = "RNA", only.pos = T, min.pct = 0.75, logfc.threshold = 0.5)
library(dplyr)
top5 <- genes.vi %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
gene.list <- unlist(top5$gene)

gene.list <- c("foxf2a", "gucy1a1", "kcnk18", "actc1a", "fsta", "smpx", "npnt", "stom", "il13ra2", "tesca", "myl9b", "acta2", "fhl1a", "krt92", "itga4", "postnb", "rgs2", "fhl3b", "pdgfra", "fgfr2", "bambia", "tmem88b", "cxcl12a", "inka1a", "smtna", "smtnb", "desmb")
pdf("visceralSMC_dotplot_v2.pdf", width = 10, height = 16)
DotPlot(obj.vi, features = gene.list, dot.scale = 10, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

FeaturePlot(obj.smc, features = c("smpx", "pln", "slitrk6", "slc1a9"))

##Figure out what genes are expressed in that transdifferentiating population of visceral SMCs
DimPlot(obj.smc, label = T)
FeaturePlot(obj.smc, c("actc1a", "fsta", "gucy1a1", "npnt"))
cells.8 <- WhichCells(obj.smc, idents = c("8"))
cells.actc1a <- WhichCells(obj.smc, idents = c("8"), expression = actc1a > 1)
cells.fsta <- WhichCells(obj.smc, idents = c("8"), expression = fsta > 1)
cells.to.move <- unlist(unique(list(cells.actc1a, cells.fsta)))
cells.middle <- setdiff(cells.8, cells.to.move)

obj.vi@meta.data$cell.annot <- NA
obj.vi@meta.data[cells.middle, "cell.annot"] <- "transdifferentiating"
obj.vi@meta.data[cells.to.move, "cell.annot"] <- "viSMC-int"
obj.vi@meta.data[WhichCells(obj.smc, idents = c("10")), "cell.annot"] <- "viSMC-sb"
obj.vi@meta.data[WhichCells(obj.smc, idents = c("13")), "cell.annot"] <- "viSMC-prog"

Idents(obj.vi) <- obj.vi@meta.data$cell.annot
genes.vi <- FindAllMarkers(obj.vi, assay = "RNA", only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
library(dplyr)
top30 <- genes.vi %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
gene.list <- unique(unlist(top20$gene))

FeaturePlot(obj.smc, c("daam1a", "ppp1r12a", "ptch2"))

##Plot a differential gene expression violin plot between all visceral SMC populations including  the unknown ones
obj.vi <- subset(obj.smc, idents = c("8", "10", "13", "16", "22", "23", "2", "11", "12", "15", "21"))
genes.vi <- FindAllMarkers(obj.vi, assay = "RNA", only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
my_levels <- c("8", "10", "13", "16", "22", "23", "2", "11", "12", "15", "21")
Idents(obj.vi) <- factor(x = Idents(obj.vi), levels = my_levels)
library(dplyr)
top5 <- genes.vi %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
gene.list <- unique(unlist(list(top5$gene, c("il13ra2", "tesca", "kcnk18", "fsta", "desmb", "smtnb", "smtna", "smpx", "col6a1", "col6a2", "nrg1", "foxl1", "hoxa13a", "hoxa13b", "evx1", "corin", "rergla", "wscd2", "sfrp2", "pax1a", "ldha", "acta2", "tagln", "fbp2", "vim", "cxcl12a", "elnb", "loxa", "fbln5", "pfkpa", "cnn1a"))))
gene.list <- c("kcnk18", "fsta", "actc1a", "acta1b", "acta1a", "cald1b", "CR753876.1", "il13ra2", "tesca", "fhl1a", "csrp1b", "nkx3.3", "prdx1", "rprmb", "hsd11b2", "mdka", "smtnb", "smtna", "zgc:195173", "vwde", "foxl1", "col4a5", "igfbp7", "col6a1", "col6a2", "nrg1", "evx1", "hoxa13a", "hoxa13b", "bambia", "lrrc17", "corin", "rergla", "krt4", "wscd2", "cyt1l", "sfrp2", "wfdc2", "pax1a", "cxcl12a", "vim", "crip1", "ldha", "pgk1", "aldocb", "eno3", "fbp2", "prss35", "rbp4", "acta2", "tagln", "si:dkey-164f24.2", "elnb", "loxa", "si:dkey-57k2.6", "fbln5", "pgam1a", "eno1a", "aldob", "ak1", "pfkpa", "cnn1a")
pdf("allSMC_vascular_visceral_vlnplot_v1.pdf", width = 10, height = 16)
VlnPlot(obj.vi, unique(gene.list), idents = factor(levels(Idents(obj.vi)), levels = c("8", "10", "13", "16", "22", "23", "2", "11", "12", "15", "21")), stack = TRUE, flip=TRUE, combine = FALSE) + NoLegend()
dev.off()


pdf("allSMC_vascular_visceral_dotplot_v1.pdf", width = 18, height = 28)
DotPlot(obj.vi, features = gene.list, dot.scale = 10, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()



##Color these SMCs with separate colors
##Find markers of these individual SMC populations
##vSMC-aorta markers
FeaturePlot(obj.smc, features = c("sfrp2", "pax1a", "si:dkeyp-106c3.1", "htra3b"))

genes.smc.1 <- markers.smc$gene[markers.smc$cluster == "2"]
genes.smc.2 <- markers.smc$gene[markers.smc$cluster == "12"]
genes.smc.3 <- markers.smc$gene[markers.smc$cluster == "15"]
genes.smc.4 <- markers.smc$gene[markers.smc$cluster == "11"]
genes.smc.5 <- markers.smc$gene[markers.smc$cluster == "21"]
genes.smc.6 <- markers.smc$gene[markers.smc$cluster == "vSMC-artery"]

##Make a slot in the main object to highlight the vascular SMC clusters
obj.smc@meta.data$smc_types <- NA
obj.smc@meta.data[WhichCells(obj.smc, idents = c("12")), "smc_types"] <- "vSMC-transient"
obj.smc@meta.data[WhichCells(obj.smc, idents = c("2")), "smc_types"] <- "vSMC-aorta"
obj.smc@meta.data[WhichCells(obj.smc, idents = c("15")), "smc_types"] <- "vSMC-OFT"
obj.smc@meta.data[WhichCells(obj.smc, idents = c("11")), "smc_types"] <- "vSMC-artery_1"
obj.smc@meta.data[WhichCells(obj.smc, idents = c("21")), "smc_types"] <- "vSMC-artery_2"

DimPlot(obj.smc, group.by = "smc_types") + NoLegend() + NoAxes()

obj.smc@meta.data$smc_types <- NA
obj.smc@meta.data[WhichCells(obj.smc, idents = c("13")), "smc_types"] <- "viSMC_progenitors"
obj.smc@meta.data[WhichCells(obj.smc, idents = c("10")), "smc_types"] <- "sb-SMCs"
obj.smc@meta.data[WhichCells(obj.smc, idents = c("8")), "smc_types"] <- "int-SMCs"

colors <- setNames(c("#FF0000", "#023858", "#800000"), c("viSMC_progenitors", "sb-SMCs", "int-SMCs"))
DimPlot(obj.smc, group.by = "smc_types", cols = colors) + NoLegend() + NoAxes()

##Get the non-neural crest trajectory of pericyte-3 and vSMCs
DimPlot(obj.smc, label = T)
obj.smc@meta.data$non_ncc <- NA
obj.smc@meta.data[WhichCells(obj.smc, idents = c("11")), "non_ncc"] <- "vSMC_artery_1"
obj.smc@meta.data[WhichCells(obj.smc, idents = c("21")), "non_ncc"] <- "vSMC_artery_2"
obj.smc@meta.data[WhichCells(obj.smc, idents = c("3")), "non_ncc"] <- "pericyte-3"

p <- DimPlot(obj.smc, group.by = "non_ncc")
LabelClusters(p, id = "ident", color = unique(ggplot_build(p)$data[[1]]$colour))

##Check co-expression between cluster 11 and cluster 3
obj <- subset(obj.smc, idents = c("11", "3"))
FeatureScatter(obj, feature1 = "acta2", feature2 = "fgf10a", pt.size = 4.0)

obj.smc <- subset(obj1, idents = c("11", "3"))
FeatureScatter(obj.smc, feature1 = "acta2", feature2 = "foxf2a", pt.size = 4.0)
FeatureScatter(obj.smc, feature1 = "ndufa4l2a", feature2 = "foxf2a", pt.size = 4.0)

##Find some more markers for intestinal SMCs that maybe spatially regionalized
##Look for markers overlapping il13ra2 but not overlapping tesca (tesca localized to mid-posterior intestine)
DimPlot(obj.smc, label = T)
FeaturePlot(obj.smc, features = c("rapgef5b", "cilp", "sptb", "npb", "htr2b", "adamtsl2", "kcnk4a", "pnck", "trpc1", "mpped1", "jph1b", "tesca"))


##Plot visceral SMC featurePlot
pond.with.grey <- c("#CECECE", "#CBDAC2", RColorBrewer::brewer.pal(9, "YlGnBu")[3:9])

FeaturePlot(obj1, c("kcnk18", "smpx", "npb", "il13ra2"), cols = pond.with.grey, pt.size = 1.2)

urd.colors <- defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)
colors_feature_plot <- c("#CECECE80", "#B2B2B280", "#7D9FD180", "#5A90E0", "#307DF0", "#0065FF", "#0078FF", "#008DFF", "#00A1FF", "#00B5FF", "#00CAFF",
  "#00DEFF", "#00F2FF", "#27FFD7", "#8CFF71", "#F1FF0D", "#FFEE00", "#FFDB00", "#FFC900", "#FFB700", "#FFA500", "#FF9200",
  "#FF8000", "#FF6D00", "#FF5B00", "#FF4800", "#FF3600", "#FF2400", "#FF1200", "#FF0000")

FeaturePlot(obj.smc, c("acta2", "tagln", "sfrp2", "elnb"), cols = pond.with.grey, pt.size = 1.2)
fplot <- FeaturePlot(obj.smc, features = c("acta2", "cald1b", "il13ra2", "npb", "fsta", "kcnk18"), pt.size = 2.2)

# Load ggplot2 library
library(ggplot2)

# Customize the plot theme
myTheme <- theme(panel.background = element_rect(fill = "white"),
                 panel.grid = element_blank(),
                 axis.line = element_line(colour = "black"),
                 axis.text = element_text(colour = "black"),
                 axis.title = element_text(colour = "black"),
                 legend.text = element_text(colour = "black"),
                 legend.title = element_text(colour = "black"))

# Modify the plot with the custom theme
fplot + scale_color_gradientn(colours = colors_feature_plot, na.value = "#CECECE")

##Ok, now plot the featureplots with the new URD color scheme
# Define the genes to plot
genes.plot <- c("ctgfa")

# Create a list of plots
fplot_peri <- lapply(genes.plot, function(gene) {
  FeaturePlot(object = obj.smc, features = gene, pt.size = 2.2)
})

# Set the gradient color scale
scale <- scale_color_gradientn(colors = colors_feature_plot, na.value = "#CECECE") 

# Customize the plot theme
theme <- theme(panel.background = element_rect(fill = "white"),
                 panel.grid = element_blank(),
                 axis.line = element_line(colour = "black"),
                 axis.text = element_text(colour = "black"),
                 axis.title = element_text(colour = "black"),
                 legend.text = element_text(colour = "black"),
                 legend.title = element_text(colour = "black"))

# Modify the plots with the custom theme and color scale
fplots <- lapply(fplot_peri, function(plot) {
  plot + scale + theme 
})

dpi <- 300
png(paste0("FeaturePlot_", genes.plot, "_fig.png"), width = 2*dpi, height = 2*dpi)
fplots
dev.off()

##Load the pericyte proportion table
dt <- readxl::read_xlsx("~/Box/Farrell Lab/pericytes_SMC_HCR/2023-03-13_pericyte_proportions.xlsx")
colnames(dt) <- c("animal", "forebrain", "hindbrain", "pharyngeal_arches", "eye")
# Melt the data frame to convert the data from wide format to long format
library(reshape2)
df_long <- melt(dt, id.vars = "animal", variable.name = "region", value.name = "proportion")
# Reshape data from wide to long format
df_long <- tidyr::gather(dt, key = "region", value = "proportion", -animal)


# Calculate mean and standard error of mean
summary_df <- aggregate(proportion ~ region, data = df_long, FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))

# Create a line graph using ggplot2
library(viridisLite)
# calculate means and standard errors
# Jitter plot with mean and error bars
# create jitter plot with mean and error bars
ggplot(df_long, aes(x = region, y = proportion, color = animal)) +
  geom_jitter(width = 0.2, height = 0, size = 2) +
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar", width = 0.2, color = "black") +
  scale_color_viridis_d() +
  labs(x = "Region", y = "Proportion of cells expressing gene") +
  theme_classic()

##Plot the number of animals that had either 1, 2, 3 or 4 epas1a+ pericytes
df <- readxl::read_xlsx("~/Box/Farrell Lab/pericytes_SMC_HCR/2023-03-13_pericyte_2_numbers_table.xlsx")
barplot(df$`number of animals`, names.arg = df$type, horiz = T)




