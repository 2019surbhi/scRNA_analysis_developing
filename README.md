# Tinglab ureter analysis


## Overview

This repository harbors analysis code and files relevent to scRNA sequencing analysis relevent to our recent work "**Ureter single-cell and spatial mapping reveal cell types, architecture, and signaling networks**". Preprint available [here](https://www.biorxiv.org/content/10.1101/2021.12.22.473889v1)

This code repository utilizes a conventional and popular published scRNA analysis tool- [Seuratv3.2.1](https://www.sciencedirect.com/science/article/pii/S0092867419305598?via%3Dihub), and other available Bioinfirmatics tools, to explore the cellular heterogenity of the human ureters (n=10) using 10X Genomics 3' single cell RNA sequencing platform.

### Seurat-based pipeline 

#### 1. Overview

The pipeline reads aligned data resulting from the cellranger's standard pipelines - [mkfastq](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count#cr-mkfastq) and [count](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count#cr-count). It runs the following set of analysis steps: 

<p> (1) Load data and create Seurat object
<p> (2) Quality Control (Filter low quality cells based on abberent mitochondrial gene%, gene counts and UMI counts)
<p> (3) Data Pre-processing (Normalization and Find Variable genes)
<p> (4) Batch correction and data integration (for sample size>1)
<p> (5) Add metadata (age,sex,etc.) - optional (provide a metadata file in xlsx or csv format with appropriate column names)
<p> (6) Pre-clustering processing - Run PCA 
<p> (7) Clustering optimization - Generate Clustree and Silhouette plots (optional)
<p> (7) Clustering - Number of PCs and resolution can be defined else runs with default value of 50 PCs and clusters at resolution 0.5
<p> (8) Post clustering analysis
       <p> &emsp; (a) Cell proportion table
       <p> &emsp; (b) Differential gene expression
       <p> &emsp; Feature plots (for top 9-20 differential markers per cluster)
       <p> &emsp; (d) Heatmap (for top15 cluster markers)
       <p> &emsp; (e) Dendrogram (again based on top15 differential markers)

#### 2.Scripts and Documentation

* [R](https://github.com/2019surbhi/tinglab_ureter_analysis/tree/main/R): `R` package code
  * [tinglab_scRNA_pipeline_functions.R](https://github.com/2019surbhi/tinglab_ureter_analysis/blob/main/R/tinglab_scRNA_pipeline_functions.R) - user defined functions of customized version of standard functions needed to run the pipeline (Please make sure to download this in the currecnt directory or add correct path to this script in the main script to source these functions). Each function provides a brief description of the required and optional variables with their default values also noted.
  *  [tinglab_scRNA_pipeline.R](https://github.com/2019surbhi/tinglab_ureter_analysis/blob/main/R/tinglab_scRNA_pipeline.R) - main script to run the pipeline from start to finish.
* [Documents](https://github.com/2019surbhi/tinglab_ureter_analysis/tree/main/Documentation): Contains instructions on running this pipeline     

#### 3. Dependencies
              
##### (1) Hardware requirements

* **System Requirements:** For our current dataset we performed this analysis on HPC (using 30 cores and 110 GB memory). The document section of this repository describes how this can be run on a similar High Performance computing systems. We will soon add smaller dataset examples that may be run on a standalone desktop. 
       For larger dataset like ours, the Seurat batch correction and integration step is the memory intensive process but the rest of the analyses post integration can be easily run on standalone laptop/PC with 4+ cores and 16+ GB RAM (Tested on macOS Catalina v10.15.17, 2.5 Ghz, 4 cores and 16GB RAM). 

* **Memory** Please note that Seurat object size ranges from few hundred MBs to several GBs depending on size of dataset (Seurat object for our ureter dataset with 10 samples is 14.33 GB in size), so please ensure there is sufficient storage memory available to save Seurat object(s) and other outputs. The pipeline allows saving either or both *integrated* and *clustered* objects and we advise users to save the *integrated* version to allow flexibilities for using this for clustering with different clustering parameters and other post clustering analysis. 
       
 ##### (2) Software requirements
              
  * This pipeline runs on `R` platform (tested on v3.2 and v4.0)
  * The users will additionally need to install R packages like `Seurat` and others pacakges mentioned in the scRNA functions [script](https://github.com/2019surbhi/tinglab_ureter_analysis/blob/main/R/tinglab_scRNA_pipeline_functions.R)
  * **Version**: To closely match our results, please run the pipeline on Seurat v3.2. However, the pipeline has been also tested on Seurat v4.0 using R v4.0 resulting in no major functional differences in the outputs. The same is true for other R pacakges used in the pipeline.
              
