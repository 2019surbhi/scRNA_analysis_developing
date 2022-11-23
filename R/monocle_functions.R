
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



## Function to transfer umap info


transfer_seurat_umap_cluster<-function(obj,cds_seurat)
{

#Set Default Assay
if('integrated' %in% names(obj@assays))
{
DefaultAssay(obj)<-'integrated'
}else
 {
     DefaultAssay(obj)<-'RNA'
 }
 
#Extract umap coordinates
umaps<-Embeddings(obj,reduction="umap")

#identical(rownames(umaps),rownames(reducedDims(cds)$UMAP)) ## TRUE

#Transfer umap coordinates
reducedDims(cds_seurat)$UMAP <- umaps

# Transfer clustering info
cds_seurat@clusters$UMAP_so$clusters <- obj@meta.data$gt_tp_cell_type_integrated_.0.5

cds_seurat <- cluster_cells(cds_seurat, reduction_method = "UMAP")

rownames(cds_seurat@principal_graph_aux$UMAP$dp_mst) <- NULL
colnames(cds_seurat@int_colData@listData$reducedDims$UMAP) <- NULL
return(cds_seurat)

}


# Functions related to finding root

## Function to calculate the root


get_earliest_principal_node<- function(cds,meta='', value='', barcodes='')
{
  if((meta!='') && (value!=''))
  {
  cell_ids <- which(colData(cds)[, meta] == value)
  }else
  {
      cell_ids<-barcodes
  }
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  return(root_pr_nodes)
}




## Function to extract barcodes of cells that coexpress top markers of a selected cluster


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



print_pseudotime<-function(cds,out_dir,f_name,dims=c(11,8.5))
{
  p<-plot_cells(cds, 
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE, 
           label_leaves = TRUE,
           label_branch_points = FALSE,
           graph_label_size = 2)
  ggsave(paste0(out_dir,f_name,'pseudotime.png'),
         width=dims[2],height=dims[1], units="in",t)
}



print_geneplots<-function(cds,out_dir,fname,dims=c(11,8.5))
{
  g<-plot_cells(cds,
           genes=gene_list,
           label_groups_by_cluster=FALSE,
           color_cells_by = "cluster")
  
  ggsave(paste0(out_dir,f_name,'geneplots.png'),
         width=dims[2],height=dims[1], units="in",g)
}

## Function to plot genes in pseudotime


get_genes_in_pseudotime_plot<-function(cds,out_dir,fname,min_exp=0.5)
{
  pg<-plot_genes_in_pseudotime(cds,
                        color_cells_by="seurat_clusters",
                        min_expr=min.exp)
  
  ggsave(paste0(out_dir,fname,'.png'),
         width=11, height=8.5,units="in",pg)
}


