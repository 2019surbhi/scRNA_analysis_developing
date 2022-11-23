# Linux or HPC dependencies
#library(argparser)
library(future)
# Data wrnaggling dependencies
library(data.table)
library(dplyr)
library(magrittr)
#'omics' dependencies
library(monocle3)
library(Seurat)
library(SeuratWrappers)
#Visualization dependencies
library(ggplot2)
library(patchwork)
library(gridExtra)
library(plotscale)
library(lattice)
library(EBImage)

source('monocle_functions.R')


input_dir<-'/Volumes/tingalab/SingleCell/CURRENT/Uro subset/2021_02_23_ureter10_uro/'

object<-'2021_02_23_ureter10_uro_clustered.rds'
run_name<-'2021_02_23_ureter10_uro_'
out_dir<-'/Users/sonas/Desktop/2021_04_06_monocle_2021_02_23_ureter10_uro/'


#Read Seurat object
obj<-readRDS(paste0(input_dir,object))

#Manually create cds object [different from Seurat wrapper]
cds<-create_cds(obj)

#Pre-process
cds<-preprocess(cds)

# Transfer Seurat UMAP info and cluster
cds<-transfer_seurat_umap_cluster(obj,cds)

# Cluster
cds<- cluster_cells(cds, reduction_method = "UMAP")

bc<-cds
cds2<-cds

# Create Trajectory [for both 'loop' and 'no loop']
# Single partition so running with partition=FALSE
cds<-learn_graph(cds,use_partition = FALSE, close_loop=FALSE)
cds2<- learn_graph(cds2,use_partition = FALSE, close_loop=TRUE)

# Print Trajectory
print_trajectory(cds,out_dir,f_name=base_name1)
print_trajectory(cds2,out_dir,f_name=base_name2)




