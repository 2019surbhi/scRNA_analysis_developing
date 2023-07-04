
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


# Functions related to Monocle cds object processing

## Function to manually create cds from seurat object

create_and_process_cds<-function(obj,verbose=FALSE)
{
#Expression Matrix
exprs<- GetAssayData(obj, slot="counts",assay = "RNA")

## phenodata
pheno.data = obj@meta.data

## feature data
genes <- data.frame(gene_short_name = rownames(exprs))
rownames(genes) <- rownames(exprs)

if(verbose)
{cat('Creating cds object','\n')}
#Create cds object
cds <- new_cell_data_set(exprs,
                         cell_metadata = pheno.data,
                         gene_metadata = genes)

if(verbose)
{cat('Preprocessing cds object','\n')}
# Preprocess - log normalization, scaling,PCA; alignment_group - do i specify it at this stage?
cds<-preprocess_cds(cds, num_dim = 50)

cds<-align_cds(cds, verbose=verbose,alignment_group = "orig.ident")

# Dimentionality reduction; Note Seurat has k=20 as default and here it is k=12, is that optimized for monocle?

cds <-reduce_dimension(
        cds,
        reduction_method="UMAP",
        preprocess_method="PCA",
        verbose=verbose)

return(cds)
}




# Functions to generate plots

## Function to Print trajectories

print_trajectory<-function(cds,out_dir,f_name,dims=c(11,8.5))
{
  t<-plot_cells(cds, 
               color_cells_by="integrated_snn_res.0.5",
               label_cell_groups = FALSE, 
               label_leaves = TRUE,
               label_branch_points = TRUE,
               graph_label_size = 3 )+
  ggtitle("seurat_clusters")

ggsave(paste0(out_dir,f_name,'trajectory.png'), width=dims[2],height=dims[1], units="in",t)
}

