# 2023_Sur_zebrafish_timecourse_scRNAseq

## 489,686 cells across 3.3–120 hpf of zebrafish development

This is the repository which contains the code that was used to generate the results and figures for the manuscript 
*"Single-cell analysis of shared signatures and transcriptional diversity during zebrafish development"*
(https://www.biorxiv.org/content/10.1101/2023.03.20.533545v1.full). 

The data is also available to navigate at our web portal Daniocell: https://daniocell.nichd.nih.gov/

![Single-cell transcriptomes were collected from whole zebrafish embryos at 50 different developmental stages (colored dots) between 14–120 hpf and then merged with our previous dataset encompassing 3.3–12 hpf (Farrell et al., 2018). Size of dots represents the number of cells recovered from each stage](~/github/2023_Sur/ZF_timecourse.jpeg)

### Overview

The code to generate the results is organized by the different steps takes to get from the raw data to the results. 
Languages and packages used are listed below:

1. R >= 4.0
2. Seurat v4.1.0
3. SeuratObject v4.0.4
4. URD v1.1.2
