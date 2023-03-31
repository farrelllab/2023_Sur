# 2023_Sur_zebrafish_timecourse_scRNAseq

## 489,686 cells across 3.3–120 hpf of zebrafish development

This is the repository which contains the code that was used to generate the results and figures for the manuscript 
*"Single-cell analysis of shared signatures and transcriptional diversity during zebrafish development"*
(https://www.biorxiv.org/content/10.1101/2023.03.20.533545v1.full). 

The data is also available to navigate at our web portal Daniocell: https://daniocell.nichd.nih.gov/

![Single-cell transcriptomes were collected from whole zebrafish embryos at 50 different developmental stages (colored dots) between 14–120 hpf and then merged with our previous dataset encompassing 3.3–12 hpf (Farrell et al., 2018). Size of dots represents the number of cells recovered from each stage](./ZF_timecourse.jpeg)

### Overview

This data was used to search for gene expression programs that are reused across multiple tissues, to determine the duration of transcriptional states during development, and identify dividing populations that are present for long developmental durations. Since global analysis of a dataset of this complexity rarely reveals its full cellular diversity, focused analyses was also performed within selected tissues - the non-skeletal muscle and the endoderm.

### Description

The code is organized into individual folders and numbered in the order they were used to generate the figures in the manuscript. Below is a description of which folders contain the code used to generate the different figure panels in the manuscript.  

 	01_Pre-processing: contains the code used to combine the previously published (Farrell et al., 2018) dataset (3.3-12 hpf) and the newly generated dataset (14-120 hpf). 
 	02_Clustering: contains the code used to generate figure panels Figure 1B, 1C, and Figure S1. 
	03_Persistent: contains the code used to generate the figure panels Figure 1F, Figure S2A-E, and Figure S3. 
		


Languages and packages used are listed below:

1. R >= 4.0
2. Seurat v4.1.0
3. SeuratObject v4.0.4
4. URD v1.1.2
