##This code is to draw a venn diagram showing the number of genes shared between human colon, small intestine - Figure 7C
## zebrafish best4+ cells

##load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(Matrix)
library(BioVenn)


##Load human colonic dataset from Smillie et al., 2019
obj <- readRDS("~/Desktop/human_intestine_data/human_colonic_dataset_Smillie_2019_seurat.rds")

##Load human small intestine object from Burclaff et al., 2022
obj.si <- readRDS("~/Box/zfext/results/2022-11_human_intestine_data/Burclaff_2022_human_small_intestine_downloaded_seurat.rds")

##Load zebrafish endoderm seurat object
obj1 <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/endoderm/obj_seurat/endoderm_seurat.rds")
markers <- readRDS("~/Box/zfext/annotations_celltype_curated_newMama/endoderm/obj_seurat/endoderm_markers.rds")

##Load the gene conversion table to convert between human and zebrafish gene names
gene.conversion.table <- read.table("~/Desktop/human_intestine_data/HGNC_human_ZF_genenamesonly.tsv.gz")

##Get cell ids representing best4 cells from human colon, SI and zebrafish intestine
cells.best4.colon <- WhichCells(obj, idents = c("Best4+ Enterocytes"))
cells.best4.si <- rownames(obj.si@meta.data)[obj.si@meta.data$leiden == "SI_BEST4"]
cells.best4.zf <- WhichCells(obj.intestine, idents = c("best4+_enterocytes"))


##Now get genes that are expressed in human colon and small intestine best4+ cells from the 2 datasets
genes.best4.colon <- unique(rownames(obj@assays$RNA@data[, cells.best4.colon])[which(obj@assays$RNA@data[, cells.best4.colon] > 0)])
genes.best4.si <- unique(rownames(obj.si@assays$RNA@data[, cells.best4.si])[which(obj.si@assays$RNA@data[, cells.best4.si] > 0)])
genes.best4.zf <- unique(rownames(obj.intestine@assays$RNA@data[, cells.best4.zf])[which(obj.intestine@assays$RNA@data[, cells.best4.zf] > 0)])

##Human small intestine genes are as ensembl ids - convert them to gene names
##Gene conversion table - from Ensembl (Biomart) that Gennady generated
gene.conversion.table.full <- read.delim("~/Box/zfext/results/2022-11_human_intestine_data/ZF_human_homologs/zebrafish-to-human_gene_conversion_biomaRt.tsv", sep = "\t")

##Directly use the table with the ensembl ids and the human gene names
gene.ensembl.ids.human.table <- read.table("~/Box/zfext/results/2022-11_human_intestine_data/ZF_human_homologs/gene_orthologs/hs_dr_orthologs_ENSEMBL-SYMBOL.tsv")

##Convert the SI genes in best4+ cells to their gene symbols
gene.best4.si.symbols <- gene.conversion.table.full[which(genes.best4.si %in% gene.conversion.table.full$Gene.stable.ID.1), "UniProtKB.Gene.Name.symbol.1"]
genes.best4.si <- setdiff(gene.best4.si.symbols, "")

##Now convert the human gene names into zebrafish gene names
genes.best4.colon.zf.converted <- unique(gene.conversion.table[gene.conversion.table$ensHS_Gene.name %in% genes.best4.colon, "lawsonZF_LLgeneSymbol"])
genes.best4.si.zf.converted <- unique(gene.conversion.table[gene.conversion.table$ensHS_Gene.name %in% genes.best4.si, "lawsonZF_LLgeneSymbol"])

##Get rid of all NAs
genes.best4.colon.zf.converted <- na.omit(genes.best4.colon.zf.converted)
genes.best4.si.zf.converted <- na.omit(genes.best4.si.zf.converted)
genes.best4.zf <- na.omit(genes.best4.zf)

##Now plot venn diagram
##Draw a venn diagram to compare how many of the genes expressed in best4+ cells are shared between 
#colon, SI and ZF intestine - Figure 6C
library(VennDiagram)
# Chart as Venn Diagram
my_col <- c("#08519c", "#FA8072", "#bd8700")
venn.diagram(
  x = list(genes.best4.colon.zf.converted, genes.best4.si.zf.converted, genes.best4.zf),
  category.names = c("1", "2", "v"),
  filename = 'Venn_human_vs_zf_best4_SI_colon_all.pdf',
  output=TRUE,
  main.cex = 0.15, main.just = c(0.25, 0.75), sub.pos = c(0.25,
                                                          0.75),
  # Output features
  imagetype="png",
  height = 1000, 
  width = 1000, 
  resolution = 500,
  compression = "lzw",
  # Circles
  lwd = 0,
  lty = 'blank',
  fill = my_col,
  
  # Numbers
  cex = .4,
  sub.cex = 0.1,
  fontface = "bold",
  fontfamily = "sans")

##Use Biovenn to generate a better venn diagram
dpi <- 300
pdf("~/Box/zfext/annotations_celltype_curated_newMama/endoderm/best4_markers/Venn_human_vs_zf_best4_SI_colon_all.pdf", width = 5*dpi, height = 5*dpi)
biovenn <- draw.venn(genes.best4.colon.zf.converted, genes.best4.si.zf.converted, genes.best4.zf, 
                     subtitle = "Venn_human_vs_zf_best4_SI_colon_all", nrtype = "abs")
dev.off()



