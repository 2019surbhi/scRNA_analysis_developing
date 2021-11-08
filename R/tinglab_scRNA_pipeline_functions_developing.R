   #last edited: 2021_10_05

# scRNA seq Seurat pipeline FUNCTIONS [based on Satija lab's Seurat tool]
# Developers: Surbhi Sona (TING lab) and Jean Clemenceau

#!/usr/bin/env Rscript

#library(argparser)
#library(future)

# Data processing libraries
library(dplyr)
library(data.table)
library(janitor)
library(stringr)

# 'Omics' related packages
library(Seurat)
library(harmony)
library(cluster)
library(factoextra)
library(clustree)
library(ggraph)

# Packages for data visualization
library(ggplot2)
library(lattice)
library(gridExtra)
library(grid)


### Function to create Seurat object ###

## This function creates Seurat object from 10X generated data.User can choose to read the .h5 matrix (default) or read from the 'filtered_feature_bc_matrix' by specifying data_dir=FALSE. ##
# sample_id - name of the sample dir
# input_dir - full path to the 10X output directory for a given sample output is stored
# data_dir - logical argument (TRUE/FALSE)

create_seurat_obj_10X<-function(sample_id,input_dir,data_dir=FALSE,verbose=FALSE)
{
 # Read data from cellranger output 
 path<-paste0(input_dir,sample_id,'/outs/')
 if(data_dir)
 { if(verbose)
     {
         cat(paste0('Loading', sample_id,'from: ',path,'filtered_feature_bc_matrix/'),sep='\n')
     }
   counts.mat <- Read10X(data.dir= paste0(path,'filtered_feature_bc_matrix/'))
 }else{
   if(verbose){
      cat(paste('Loading',sample_id,'from',path,'filtered_feature_bc_matrix.h5',sep=' '),sep='\n')
    }
   counts.mat<-Read10X_h5(filename = paste0(path,'filtered_feature_bc_matrix.h5'))
  }
 
# Create Seurat object
 if(verbose)
    {cat(paste0('Creating Seurat Object for ',sample_id),sep='\n')}

 s.obj<-CreateSeuratObject(counts=counts.mat, project = sample_id)
 s.obj<-add_mt_perc(s.obj,sample_id)

return(s.obj) 

}

### Funtion to create Seurat object from counts data ###

## This function creates Seurat object usig count matrix. ##
# count_dat - name of the sample dir
# sample_id - Sample name

create_seurat_obj_from_counts_data<-function(counts_dat,sample_id,verbose=FALSE)
{
# Create Seurat object
 if(verbose){
    cat(paste0('Creating Seurat Object for sample ',sample_id),sep='\n')
  }

 s.obj<-CreateSeuratObject(counts=counts_dat, project = sample_id)

#Add percentage of mitochondrial genes
s.obj<-add_mt_perc(s.obj,sample_id)

return(s.obj)

}

### Function to add percentage of mitochondrial genes ###

## This function computes % of mitochondrial gene expressed in a given Seurat obj and adds this info back to the Seurat obj. It considers 'MT-' as well as 'mt-' notation defining human and mouse mitochondrial genes, respectively, allowing usage of this pipeline for analysis of both human and mouse scRNA data ##

# sample_id - sample name (must match the sample name defined in the 10X workflow)

add_mt_perc<-function(s.obj,sample_id,verbose=FALSE)
{
  if(verbose)
    {cat(paste0('Calculating MT gene % for ',sample_id),sep='\n')}

if(sum(PercentageFeatureSet(s.obj,pattern = "^MT-"))==0)
 { #In case of mouse dataset, '^mt-' is used
  if(sum(PercentageFeatureSet(s.obj,pattern = "^mt-"))==0)
  { if(verbose)
    {cat("No mitochondrial genes found ", '\n')}
  }else
    {
     s.obj[['perc.mt']]<-PercentageFeatureSet(s.obj,pattern = "^mt-")
    }
 }else
   {
    s.obj[['perc.mt']]<-PercentageFeatureSet(s.obj,pattern = "^MT-")
   }

return(s.obj)

}

### Functions to generate histogram ###

## This function generates histogram for a specified feature in a given seurat object ##
# s.obj - Seurat object to create histogram for
# out_dir - directory path to save histograms at
# prefix - Any prefix (other than sample name) to add to each histogram output
# x.lab - specify which feature to print the histogram for, e.g. LibrarySize or GeneCounts (also X axis label)
# xlim - upper limit of histogram (leave xlim=0 if no upper threshold is to be applied on the final printed output)

get_histogram<-function(s.obj,out_dir,prefix,x.lab,xlim=0,cutoff='none')
{
 sample.mat<-s.obj@assays$RNA@counts

 s_id<-s.obj@project.name
 name<-paste0(prefix,s_id)

 if(x.lab=="LibrarySize")
     {
      x.lim<-c(0,xlim)
      hist_dat<-Matrix::colSums(sample.mat)
     }else if(x.lab=="GeneCounts")
        {
         x.lim<-c(0,xlim)
         hist_dat<-Matrix::colSums(sample.mat>0)
        }else if(x.lab=="MtPerc")
           {
            x.lim<-c(0,xlim)
            hist_dat<-s.obj$perc.mt
           }

# Reset x.lim if no x_lim argument entered
if(xlim==0)
    {x.lim=''}

print_histogram_abline(dat=hist_dat,out_dir,label=name,x_lab=x.lab,cutoff=cutoff,x=x.lim)

return(hist_dat)

}

## This function prints histogram for a given data with a red line specifying the user defined filtering cutoff ##
# dat - data to plot histogram for (gene count/library size/ mitochondrial percentage)
# out_dir - directory path to save histograms at
# label - Any prefix (other than sample name) to add to each histogram output
# x_lab - specify which feature to print the histogram for, e.g. LibrarySize or GeneCounts (also X axis label)
# xlim - upper limit of histogram (leave xlim=0 if no upper threshold is to be applied on the final printed output)

print_histogram_abline<-function(dat,out_dir,label,x_lab,cutoff='none',x='',b=1000)
{
 
 hist.title<-paste(label,x_lab,sep='_')
 options(bitmapType='cairo')
 png(paste0(out_dir,hist.title,".png"))
 
 #Procure xlim if not defined
 if(x=='')
  {   xu<-range(dat)[2]
      x<-c(0,xu)
  }
  

  if(cutoff=='none')
  { cat('printing histogram \n')
    hist(dat,main=hist.title,breaks=b,xlab=x_lab,xlim=x)
  }else{
    cat('printing histogram with abline \n')
    hist(dat,main=hist.title,breaks=b,xlab=x_lab,xlim=x)
    abline(v = cutoff, col=c("red","red"), lwd=c(2,2))
    cat('finished printing histogram with abline \n')
  }
dev.off()

}

### Function to record raw and threshold cell count ###

## This function calculates filtered cell count for each sample after applying each threshold and combination of thresholds to give an idea how each feature based filtering is affecting final cell count ##
# s. obj - seurat obj to calculate cell count table for [Note: the cells are NOT actually filtered by this function]
# out_dir - path where the cell count table will be saved
# mt.thres - Mitchondrial %  (or % mt) threshold (cells with % mt > than the value will be filtered out)
# genecnt.thres - a vector of 2 integers definining lower and higher gene count thresholds for filtering
# libsize.thres - a vector of 2 integers definining lower and higher library size thresholds for filtering

get_cell_table<-function(s.obj,out_dir,mt.thres, genecnt.thres, libsize.thres, verbose=FALSE)
{

if(verbose){
    cat(paste0('Calculating raw and post thresholding cell counts for ',s.obj@project.name),sep='\n')
  }

# Get raw cell count
 raw.cell.cnt<-ncol(s.obj)
 
# Get cell counts after applying one or more thresolds
 mt.thres.cells<-subset(s.obj, subset=perc.mt<mt.thres) %>% ncol()
 
 gene.thres.cells<-subset(s.obj,subset=nFeature_RNA>genecnt.thres[1] & nFeature_RNA<genecnt.thres[2]) %>% ncol()
 
 libsize.thres.cells<-subset(s.obj,subset=nCount_RNA>libsize.thres[1] & nCount_RNA<libsize.thres[2]) %>% ncol()
  
 gene_lib.thres.cells<-subset(s.obj,subset= nFeature_RNA>genecnt.thres[1] & nFeature_RNA<genecnt.thres[2] & nCount_RNA>libsize.thres[1] & nCount_RNA<libsize.thres[2]) %>% ncol()

 gene_mt.thres.cells<-subset(s.obj, subset= nFeature_RNA>genecnt.thres[1] & nFeature_RNA<genecnt.thres[2] & perc.mt<mt.thres) %>% ncol()

 lib_mt.thres.cells<-subset(s.obj,subset= nCount_RNA>libsize.thres[1] & nCount_RNA<libsize.thres[2] & perc.mt<mt.thres) %>% ncol()

 gene_lib_mt.thres.cells<-subset(s.obj, subset= nFeature_RNA>genecnt.thres[1] & nFeature_RNA<genecnt.thres[2] & nCount_RNA>libsize.thres[1] & nCount_RNA<libsize.thres[2] & perc.mt<mt.thres) %>% ncol()

if(verbose){
    cat(paste0(s.obj@project.name,'- Creating cell table'),sep='\n')
  }

cell_table<-data.frame(raw.cell.cnt, mt.thres.cells,gene.thres.cells,libsize.thres.cells,gene_lib.thres.cells, gene_mt.thres.cells, lib_mt.thres.cells, gene_lib_mt.thres.cells)

row.names(cell_table)<-s.obj@project.name

return(cell_table)

}

### Function to filter low quality cells; returns filtered seurat object and a dataframe with cell count pre and post each filtering step ###
# s.obj - seurat obj to apply the filters to
# mt.thres - Mitchondrial %  (or % mt) threshold (cells with % mt > than the value will be filtered out)
# genecnt.thres - a vector of 2 integers definining lower and higher gene count thresholds for filtering
# libsize.thres - a vector of 2 integers definining lower and higher library size thresholds for filtering

filter_cells<-function(s.obj,mt.thres, genecnt.thres, libsize.thres, verbose=FALSE)
{
 if(verbose)
   { 
    cat(paste0(s.obj@project.name,' - Apply Thresholds:'),sep='\n') 
   }
  
 # Get raw  cell and gene metrics
  raw.cell.cnt<-ncol(s.obj)
  raw.min.genes <-min(s.obj[['nFeature_RNA']])
  raw.max.genes <-max(s.obj[['nFeature_RNA']])
  raw.total.genes <- nrow(s.obj)

# Apply threholds

s.obj<-subset(s.obj, subset= nFeature_RNA>genecnt.thres[1] & nFeature_RNA<genecnt.thres[2])
s.obj<-subset(s.obj,subset= nCount_RNA>libsize.thres[1] & nCount_RNA<libsize.thres[2])
s.obj<-subset(s.obj,perc.mt<mt.thres)

 #Get thresholded metrics
  thr.cell.cnt<-ncol(s.obj)
  thr.min.genes <-min(s.obj[['nFeature_RNA']])
  thr.max.genes <-max(s.obj[['nFeature_RNA']])
  thr.total.genes <- nrow(s.obj)

cell.data<-data.frame(
  raw.cell.cnt=raw.cell.cnt,
  thr.cell.cnt=thr.cell.cnt,
  raw.min.genes=raw.min.genes,
  thr.min.genes=thr.min.genes,
  raw.max.genes=raw.max.genes,
  thr.max.genes=thr.max.genes,
  raw.total.genes=raw.total.genes,
  thr.total.genes=thr.total.genes,
  row.names=s.obj@project.name)

obj_list<-list(s.obj,cell.data)

return(obj_list)
}


# Function to plot data before and after applying thresholds -  from Jean Clemenceau
plot_threshold_effects<-function(plot.data,thresholds,title.root)
{  
  sample_ids<-rownames(plot.data)
 #Plot sample cell counts pre/post thresholding
  cell.cnt.thres.plot<-ggplot(plot.data)+
    geom_point(aes(x=rank,y=raw.cell.cnt,color='raw'))+
    geom_point(aes(x=rank,y=thr.cell.cnt,color='thresholded'))+
    theme_light()+
    scale_x_continuous(breaks=plot.data$rank,labels=sample_ids)+
    scale_color_manual(name='Filter Status',values =c('raw'='blue','thresholded'='darkorange'), labels = c('Raw','Thresholded'))+
    theme(axis.text.x = element_text(angle = 30,hjust = 1),legend.position = 'bottom')+
    ggtitle(label=paste0(title.root,' - Sample Cell Counts'),subtitle = paste0('Thresholds: %Mt=',thresholds[1],', gene counts=[',thresholds[2],',',thresholds[3],']' ))+
    labs(x="Ranked Samples", y="Cell Count")

  #Plot sample gene counts pre/post thresholding
  gene.cnt.thres.plot<-ggplot(plot.data)+
    geom_point(aes(x=rank,y=raw.total.genes,color='raw'),alpha=0.5)+
    geom_point(aes(x=rank,y=thr.total.genes,color='thresholded'),alpha=0.5)+
    theme_light()+
    scale_x_continuous(breaks=plot.data$rank,labels=sample_ids)+
    scale_color_manual(name='Filter Status',values =c('raw'='blue','thresholded'='darkorange'), labels = c('Raw','Thresholded'))+
    theme(axis.text.x = element_text(angle = 30,hjust = 1),legend.position = 'bottom')+
    ggtitle(label=paste0(title.root,' - Gene Counts'),subtitle = paste0('Thresholds: %Mt=',thresholds[1],', gene counts=[',thresholds[2],',',thresholds[3],']' ))+
    labs(x="Ranked Samples", y="Total Gene Count")
  #Plot sample gene count ranges pre/post thresholding
  gene.cnt.range.thres.plot<-ggplot(plot.data)+
    geom_point(aes(x=rank,y=raw.min.genes,color='raw',fill='raw'),alpha=0.5,shape=24)+
    geom_point(aes(x=rank,y=thr.min.genes,color='thresholded',fill='thresholded'),alpha=0.5,shape=24)+
    geom_point(aes(x=rank,y=raw.max.genes,color='raw',fill='raw'),alpha=0.5,shape=25)+
    geom_point(aes(x=rank,y=thr.max.genes,color='thresholded',fill='thresholded'),alpha=0.5,shape=25)+
    theme_light()+
    scale_x_continuous(breaks=plot.data$rank,labels=sample_ids)+
    scale_color_manual(name='Filter Status',values =c('raw'='blue','thresholded'='darkorange'), labels = c('Raw','Thresholded'))+
    scale_fill_manual(name='Filter Status',values =c('raw'='blue','thresholded'='darkorange'), labels = c('Raw','Thresholded'))+
    theme(axis.text.x = element_text(angle = 30,hjust = 1),legend.position = 'bottom')+
    ggtitle(label=paste0(title.root,' - Gene Count Range'),subtitle = paste0('Thresholds: %Mt=',thresholds[1],', gene counts=[',thresholds[2],',',thresholds[3],']' ))+
    labs(x="Ranked Samples", y="Gene Count")

  plots<-list(cell.cnt.thres.plot, gene.cnt.thres.plot, gene.cnt.range.thres.plot)
  return(plots)
}



### Funtion to pre-process data ###
## This function processes seurat obj prior to batch correction or other downstream analysis. The seurat obj is normalized and variable genes are calculated for a given sample ##

# s.obj - Seurat obj
# hvg - number of variable genes to be calculated (default: 2000)

pre_process<-function(s.obj,hvg=2000, verbose=FALSE)
{
 
#Log normalize data
if(verbose)
   {
    cat(paste0('Normalizing data for ',s.obj@project.name),sep='\n')
   }
 s.obj<-NormalizeData(s.obj,normalization.method = 'LogNormalize',scale.factor = 10^4, verbose = verbose)

#Select genes that excibit the most variance across cells
 
if(verbose){
      cat(paste0(' Getting Variable Genes for ',s.obj@project.name),sep='\n')
   }

s.obj<-FindVariableFeatures(s.obj,selection.method = 'vst', nfeatures=hvg, verbose = verbose)

return(s.obj)

}

### Function to create variable genes table  ### - from Jean Clemenceau
## This function creates variable gene table for a given seurat obj. The variable genes must be already computed for this seurat object ##
# s.obj - Seurat obj
# out_dir = path to directory to output the tables

get_var_genes<-function(s.obj,out_dir='./',verbose=FALSE)
{
 
  f.name<-paste0(out_dir,s.obj@project.name,'-variableGenes.tsv')
  
  #Rank and sort variable genes
  var.gen.data<-HVFInfo(s.obj)
  var.gen.data<-var.gen.data[order(var.gen.data[,3],decreasing = TRUE),]
  
  #Export to tsv file
  fwrite(var.gen.data,f.name,sep='\t',eol ='\n',col.names = TRUE,row.names = TRUE)
}



### Funtions to integrate multiple samples with batch correction ###
## This function uses a popular and robust Seurat batch correction method (CCA-MNN).
# s.obj.list - list of seurat objects to be integrated after batch correction
# merged.title - project name of the integrated object
# genes - Specify the genes to be integrated (if not specified, only sample anchors included in the final expression matrix of the integrated object). Other options are: 'all'  to integrate all genes or a number to specify the number of top variable genes to include (the specified number of variable genes are calculated on merged obj)

cca_batch_correction<-function(s.obj.list,project.name,anchors=2000,int.genes='',verbose=FALSE)
{
  if(verbose){
    cat("Performing CCA-MNN Pipeline",sep='\n')
  }

# Extract the normalization method [assuming all objects have same normalization method]
norm<-s.obj.list[[1]]@commands$NormalizeData.RNA@params$normalization.method

#Determine min number of neighbors
 mink<-min(200, min(sapply(seq_along(s.obj.list),function(x) ncol(s.obj.list[[x]]) ))  )
 
# Select features
features<-SelectIntegrationFeatures(obj.list,nfeatures=anchors)

# Find Integration anchors
sample.anchors<-FindIntegrationAnchors(s.obj.list,dims = 1:30,k.filter=mink,reduction='cca',anchor.features=features,normalization.method=norm)

# Integrate with specified number of genes
 if(int.genes=='')
{if(verbose)
    {cat("Integrating the sample.anchors only",sep='\n')}
 s.obj.integrated<-IntegrateData(anchorset=sample.anchors, dims=1:30,normalization.method=norm)
}else if(int.genes=='all')
   {
    all.genes <- lapply(s.obj.list, row.names) %>% Reduce(intersect, .)
    
     if(verbose)
       {cat("CCA_MNN batch correction - integrating all genes", sep='\n')
        cat("Total genes being used for integration = ", length(all.genes),'\n')
       }
    s.obj.integrated<-IntegrateData(anchorset=sample.anchors, dims=1:30,features.to.integrate=all.genes,normalization.method=norm)
   }

rm(sample.anchors)

s.obj.integrated@project.name<-project.name

 return(s.obj.integrated)
}

### Function to add metadata to Seurat object ###
## Note: This function needs to be adapted as per your metafile format (add or modify column names, as needed) ##
# s.obj - seurat obj to add the meta data to
# meta_file - file containing metadata (stored in the format: rows = sample and column = meta data)

add_metadata<-function(s.obj,meta_file)
{
  #Import metadata
  
  metadata<- fread(meta_file,header = TRUE,stringsAsFactors = TRUE)
  indx<- match(s.obj@meta.data$orig.ident,metadata$sample_id)
  
  #Add Sex data
  s.obj[['sex']]<- metadata$Sex[indx]
  
  #Add Age data
  s.obj[['age']]<- metadata$Age[indx]
  
  #Add Normal/Tumor data
  s.obj[['normal_tumor']]<- metadata$Normal_Tumor[indx]
  
  return(s.obj)
}

### Function to plot variable genes in samples ### - from Jean Clemenceau
## This function generates variable genes dispersion plot for a given Seurat obj. Note: the variable genes must be already computed for this Seurat obj ##
# s.obj - Seurat obj

get_var_genes_plot<-function(s.obj,verbose=FALSE)
{ 
 #Determine top 10 most variable genes
  top10var<- head(VariableFeatures(s.obj),10)
  
  #Plot variable features
  gene_var_plot <-VariableFeaturePlot(s.obj)+
    labs(title=s.obj@project.name)+
    theme(legend.position = c(0,0.95))
  gene_var_plot <-LabelPoints(plot = gene_var_plot, points = top10var, repel = TRUE)
 
  return(gene_var_plot)
}

### Function to plot symmetric UMAPs ###

## This function plots symmetric UMAPs for a given Seurat obj ##
# s.obj - Seurat obj
# umap_col - if specified, these colors will replace the default ggplot colors for clusters
# label - logical argument (TRUE/FALSE) to specify whether to print labels on UMAP or not
# title - UMAP title (printed on UMAP)
# group - grouping variable. By defualt the cells are colored by clusters (e.g.  user can specify to color the cells by 'sample')
# split -
# dot - size of dot (i.e. cell) on UMAP [Default: 0.3]
# save - logical argument (TRUE/FALSE) to specify whether or not to save the UMAP. By default, the UMAP is returned.
# out - path to output directory
# file_tag - file name prefix
# h = height  of the final output [Default: 8 inch]
# w = width of the final output [Default: 8 inch]
# reduction - dimentionality reduction to use to plot UMAP [Default: 'umap']

customized_umap<-function(s.obj,umap_cols=NULL,label=FALSE,title=NULL, group=NULL,split=NULL,dot=0.3,save=FALSE,out='./',file_prefix,h=8,w=8,reduction='umap')
{
    ns<-s.obj@meta.data$sample %>% unique %>% length()
    nc<-ncol(s.obj)

    if(is.null(title)==TRUE)
    {
      title<-str_sub(s.obj@project.name,end=-2)
    }

    if(ns>1)
    {
      sub<-paste0('n_samples= ',ns,'\n',
                                'n_cells= ', nc)
    }else
    {
      sub<-paste0('n_cells= ', nc)
    }

    umap<-DimPlot(s.obj, group.by=group, split.by=split,pt.size=dot,cols=umap_cols, label=label,reduction=reduction)

    r1<-range(s.obj@reductions[[reduction]]@cell.embeddings[,1])
    r2<-range(s.obj@reductions[[reduction]]@cell.embeddings[,2])

    sc<-c((floor(min(r1[1],r2[1]))),
          ceiling((max(r1[2],r2[2]))))

    umap<-umap+theme(legend.position = 'bottom')+
      scale_x_continuous(limits=sc)+
      scale_y_continuous(limits=sc)+
      ggtitle(label=title,
              subtitle = sub)
if(is.null(group)==FALSE)
{
    file_prefix<-paste0(file_prefix,group,'_')
}
if(save==TRUE)
{
    ggsave(paste0(out,file_prefix,'UMAP.png'),
               width=w, height=h, units="in",umap)

}
    return(umap)
    
}

### Function to tabulate differentially expressed genes for all clusters ###

## This function uses FindMarkers() to generate and save a list of differentially expressed genes  for all clusters in a given Seurat obj. The output tables are filtered for p_val_adj<0.0.5 and sorted by pct.diff and avg_logFC ##
# s.obj - Seurat obj
# clusters - list of clusters to compute the differential genes for
# out_dir - patht to output the differential gene tables
# file_tag - file prefix which will be added to all files

differential_gene_exp<-function(s.obj,clusters,out_dir='./',file_prefix="",verbose=FALSE,save=TRUE)

{if(verbose)

{cat("Finding Differentially expressed cluster markers", '\n')}
  
 DefaultAssay(s.obj)<-"RNA"
 markers<-list()
 
#Create directory to export markers
  if(!dir.exists(paste0(out_dir,'diff_genes/')))
        {dir.create( paste0(out_dir,'diff_genes/'))}
    out_path<-paste0(out_dir,'diff_genes/')

# Find markers and export to tsv file
seurat_clusters<-grep('snn',colnames(s.obj@meta.data))
  for(i in 1:length(clusters))
  {if( sum(s.obj@meta.data[[seurat_clusters]]==clusters[i])<4)
    {
     cat('Cluster',clusters[i],'has too few cells','\n')
    }else{
   	  markers[[i]]<-FindMarkers(s.obj, ident.1= clusters[i])
    	  gene<-rownames(markers[[i]])
    	  markers[[i]]<-cbind(gene,markers[[i]])
    	  markers[[i]]<-markers[[i]] %>% filter(p_val_adj<0.05) %>% mutate(pct.diff=pct.1-pct.2)%>% arrange(desc(pct.diff,avg_logFC))
   #markers[[i]]<-markers[[i]][order(markers[[i]]["avg_logFC"], decreasing=TRUE),]
    	  if(save)
		{fwrite(markers[[i]], paste0(out_path,file_prefix,'cluster',clusters[i],'_diff_markers.tsv') ,append=FALSE,sep='\t',row.names=FALSE,col.names=TRUE )}  
	}
  } 
 return(markers)
}


### Function to generate Features plot ###

## This function prints Feature Plots using specified set of top makers or genes for a given seurat obj ##
# features.list - vector of sorted genes (preferrably differential markers)
# s.obj - Seurat obj
# top - number to specify how many of the sorted genes to generate the Feature Plots for (Default: 20)

get_features_plot<-function(features.list,s.obj, top=20,verbose=FALSE)
{
 if(verbose)
	{cat("Generating FeaturePlots for all clusters",'\n')}
 DefaultAssay(s.obj)<-"RNA"
 features<-head(features.list,top)
 
 plot.list<-list()

 plot.list<-lapply(X=features, FUN=function(x) {FeaturePlot(s.obj,feature=x, min.cutoff='q30')})
 
 return(plot.list)
}

### Funtion to generate dendrogram ###

## This function generates 2 dendrograms for a given Seurat obj: (1) dendrogram generated using all genes and (2) dendrogram generated using specified set of differential features (or genes)  ##
# s.obj - Seurat obj
# file_tag - file name prefix
# features - vector of genes to plot the dendrogram for [This argument is used only to generate the second dendrogram]
# out_dir - path to output directory

plot_dendrogram<-function(s.obj,file_tag,features, out_dir='./',verbose=FALSE)
  {
  if(verbose)
    {cat("Generating dendrogram",'\n')}
  clus<-levels(s.obj)
  collapsed_norm_mat<-vector()

  for(i in 1: length(clus))
    {
      sub<-NULL
      sub<-subset(s.obj,idents = clus[i])
      sub_mat<-as.matrix(sub@assays$RNA@data)
      rm(sub)

      clus_sum<-rowSums(sub_mat)
      rm(sub_mat)
      
      collapsed_norm_mat<-cbind(collapsed_norm_mat,clus_sum)
      rm(clus_sum)
    }


  colnames(collapsed_norm_mat)<-clus

#Filter genes, not expressed in any cluster
  filter<-which(rowSums(collapsed_norm_mat)==0)
  collapsed_norm_mat<-collapsed_norm_mat[-filter,]

#Calculate dissimilarity [default: Euclidean] - all markers
options(bitmapType='cairo')
png(paste0(out_dir,file_tag,"dendrogram.png"),width=11,height = 8.5, units="in",res=300)
  hc_avg<-hclust(as.dist(1-cor(collapsed_norm_mat)),method = "average")
  plot(hc_avg, main=paste0(file_tag,"Cluster Dendrogram"), xlab="")
  dev.off()

# Select diff markers
  m<-which(rownames(collapsed_norm_mat) %in% features==TRUE)
  collapsed_norm_diff_markers<-collapsed_norm_mat[m,]

  options(bitmapType='cairo')
  png(paste0(out_dir,file_tag,"diff_markers_dendrogram.png"),width=11,height = 8.5, units="in",res=300)
  hc2_avg<-hclust(as.dist(1-(cor(collapsed_norm_diff_markers))),method = "average")
  plot(hc2_avg)
  dev.off()
}

### Functions to iteratively cluster seurat obj ###
## This function iteratively clusters seurat obj based on a set of user-defined clustering parameters. Each round of clustering info is stored in the metadata slot of the seurat obj and later used for generating clustree ##
# s.obj - Seurat obj
# res = vector definining the range of resolutions to be used for clustering for the defined set of PCs
# dims - range of PCs to be used for clustering
# reduction - dimentionality reduction to be used (Default: 'pca')
# assay - which assay to use (Default: 'integrated')

iterative_clus_by_res<-function(s.obj,res, dims_use,reduction='pca',assay='integrated',verbose=FALSE)
{ if(verbose)
   {cat("Performing iterative clustering by resolution for PCs 1:",max(dims_use),'\n')}
  DefaultAssay(s.obj)<-assay

  s.obj<-FindNeighbors(s.obj,dims=dims_use,reduction=reduction,assay=assay)
  for(i in 1:length(res))
  {
    s.obj<-FindClusters(s.obj, res=res[i])
  }
  return(s.obj)
}

### Function to generate clustree ###
## This function uses 'clustree' pacakge to generate clustree outputs that illustrate the distribution of cells at varying resolution for when clustered using a specific set of PCs ##
# s.obj - Seurat obj
# prefix - the prefix to look for in the seurat obj metadata that defines the clusters cplumn (e.g. 'integrated_snn_res.')


print_clustree_png<-function(s.obj,prefix,out_dir='./',file_prefix,verbose=FALSE)
{ if(verbose)
   {cat("Generating clustree png(s)", sep='\n')}
   if(!dir.exists(paste0(out_dir,'clustree/')))
  {dir.create(paste0(out_dir,'clustree/'))}
 out_path<-paste0(out_dir,'clustree/')
 options(bitmapType='cairo')
 png(paste0(out_path,file_prefix,"_clustree.png"), width=8.5, height=11, units="in",res=300)
  obj_clustree<-clustree(s.obj,prefix=prefix)
  print(obj_clustree)
  dev.off()
}

## Function to generate clustree geneplots ###
## This function uses 'clustree' pacakge to generate clustree geneplots that illustrate the disctibution of cells at varying resolution for when clustered using a specific set of PCs with the mean/median expression of genes plotted on top of this information. ##

# s.obj - Seurat obj
# prefix - the prefix to look for in the seurat obj metadata that defines the clusters cplumn (e.g. 'integrated_snn_res.')
# assay - which assay to look into for plotting gene expression (Default: 'integrated')
# fun_use - specify whether to plot mean/median expression of gene (Default: 'median')
# out_dir - output directory path to store the plots
# file_prefix - file prefix name

print_geneplots_on_clustree<-function(s.obj,genes, prefix='integrated_snn_res.',assay='integrated', fun_use='median', out_dir, file_prefix, verbose=FALSE)
{
 validated_genes<-genes[which(genes %in% rownames(s.obj)==TRUE)]
 if(length(validated_genes)>0)
{ if(verbose)
   {cat("Generating clustree geneplots",sep='\n')}
  
 if(!dir.exists(paste0(out_dir,'ƒgeneplots/',assay,'/')))
  {dir.create(paste0(out_dir,file.path("clustree_geneplots",assay)),recursive=TRUE)}
 out_path<-paste0(out_dir,'clustree_geneplots/',assay,'/')
 
 DefaultAssay(s.obj)<-assay
 cat('validated genes: \n')
 head(validated_genes)
 
 for(i in 1:length(validated_genes))
 {
  clustree(s.obj, prefix=prefix,node_colour = validated_genes[i], node_colour_aggr = fun_use)
  ggsave(paste0(file_prefix,'_',assay,'_',validated_genes[i],".png"),path=out_path, width=8.5, height=11,units="in")
 }
}else
  {cat("None of the requested gene(s) found!",'\n')}
}


### Function to generate Silhouette plots ###
## This function generates Silhouette plots for clusters of a give seurat obj using 'cluster' and 'factoExtra' packages ##
# s.obj - Seurat obj
# reduction - which dimentionality reduction to use (Default: 'pca')
# dims - range of PCs that were used for clustering the seurat obj (Default: 1:50)
# out_dir - output path where plots are saved
# file_prefix - output file prefix

get_silhouette_plot<-function(s.obj,reduction='pca',dims=1:50,out_dir='./',file_prefix,verbose=FALSE)
{
 if(verbose)
  {cat('Generating Silhouette plots','\n')} 
 if(!dir.exists(paste0(out_dir,'Silhouette/')))
  {dir.create(paste0(out_dir,'Silhouette/'))}
 out_path<-paste0(out_dir,'Silhouette/')
  clusters<-s.obj@meta.data[,grep('snn',colnames(s.obj@meta.data))]
  
  dist.matrix<-dist(x=Embeddings(object=s.obj[[reduction]])[,dims])
  sil<-silhouette(x=as.numeric(as.factor(clusters)),dist=dist.matrix)
  
  #Printing the Silhouette plot
  fviz_silhouette(sil)
  ggsave(filename=paste0(file_prefix,"_silhouette.png"),path=out_path, width=33,height=10)
}
