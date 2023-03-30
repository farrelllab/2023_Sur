# Gene module identification using Fuzzy c-means clustering

This code was designed to catalog gene expression programs (GEPs) that are shared between two or more tissues. In total, this approach resulted in 147 GEPs out of which 57 GEPs (i) represented outliers, (ii) consisted of member genes expressed only within one tissue subset, and (iii) and contained member genes that were mitochondrial or ribosomal. These GEPs were excluded from the analysis. 

This folder contains the code that follows the following steps:
1. smoothed the data from all cells in our dataset using a 5 nearest-neighbor network to reduce effects of technical noise 
2. used fuzzy c-means (FCM) clustering to group genes with similar expression over time and across tissues 
3. filtered out poor quality GEPs (<5 member genes or primarily technical member genes), and 
4. filtered out GEPs that were expressed in a single tissue.

As a complementary approach, we also calculated GEPs on tissue-specific subsets and used correlation with cosine distance to identify modules that were found in multiple tissues; while this approach identified many more cell-type specific GEPs, it did not yield any additional shared GEPs that were not already captured from the global analysis 

Finally, all shared GEPs were plotted as a heatmap using the "pheatmap" package and hierarchically clustered using the "hclust" function. Also for select GEPs included in the manuscript, dotplots for the member genes were plotted. 

Files contained in this folder contains code that was used to generate the following figures: Figure 2A-D, and Figure S4A-B
