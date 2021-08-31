#last edited: 2021_08_24

# scRNA seq Seurat pipeline FUNCTIONS [developer: Surbhi Sona, TING lab]
# Objective: Processing scRNA seq data on HPC using Seurat workflow

###################
#SCRIPT begins here
###################

#!/usr/bin/env Rscript

library(argparser)
library(future)

# Data processing libraries
library(dplyr)
library(data.table)
library(janitor)
library(stringr)

# 'Omics' related packages
library(Seurat)
library(cluster)
library(factoextra)
library(clustree)
library(ggraph)

# Packages for data visualization
library(ggplot2)
library(lattice)
library(gridExtra)
library(grid)


### FUNCTIONS ###

# Function to create Seurat object #

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

#Funtion to create Seurat object from counts data #

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

# Function to add percentage of mitochondrial genes #

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


# Function to record raw and threshold cell count # - from Jean Clemenceau
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

# Functions to generate histogram #
print_histogram<-function(dat,out_dir,label,x_lab,x='',b=1000)
{
 #cat("xlim= ",x,'\n')
 #c.xlim=class(x)
 #cat("class of xlim= ",c.xlim)
 hist.title<-paste(label,x_lab,sep='_')
 options(bitmapType='cairo')
 png(paste0(out_dir,hist.title,".png"))
 if(x=='')
  {
   hist(dat,main=hist.title,breaks=b,xlab=x_lab)
  }else
    {
     hist(dat,main=hist.title,breaks=b,xlab=x_lab,xlim=x)
    }
dev.off()
}

print_histogram_abline<-function(dat,out_dir,label,x_lab,cutoff,x='',b=1000)
{
 #cat("xlim= ",x,'\n')
 #c.xlim=class(x)
 #cat("class of xlim= ",c.xlim)
 hist.title<-paste(label,x_lab,sep='_')
 options(bitmapType='cairo')
 png(paste0(out_dir,hist.title,".png"))
 if(x=='')
  {
   hist(dat,main=hist.title,breaks=b,xlab=x_lab)
   abline(v = cutoff, col=c("red","red"), lwd=c(2,2))
   # text(x=cutoff, y=ylim[1], labels=cutoff, adj=c(1.1,1), col='blue')
   }else
    {
     hist(dat,main=hist.title,breaks=b,xlab=x_lab,xlim=x)
     abline(v = cutoff, col=c("red","red"), lwd=c(2,2))
     #text(x=cutoff, y=ylim[1], labels=cutoff, adj=c(1.1,1), col='blue')
    }
dev.off()
}

get_histogram<-function(s.obj,out_dir,prefix,x.lab,xlim=0)
{
 sample.mat<-s.obj@assays$RNA@counts

 s_id<-s.obj@project.name
 name<-paste0(prefix,s_id)

 if(x.lab=="LibrarySize")
     {
      x.lim<-c(0,xlim)
      hist_dat<-Matrix::colSums(sample.mat)
      #print_histogram(hist_dat,out_dir,label=name,x.lab,x.lim)
     }else if(x.lab=="GeneCounts")
        {
         x.lim<-c(0,xlim)
         hist_dat<-Matrix::colSums(sample.mat>0)
         #print_histogram(hist_dat,out_dir,label=name,x.lab,x.lim)
        }else if(x.lab=="MtPerc")
           {
            x.lim<-c(0,xlim)
            hist_dat<-s.obj$perc.mt
            #print_histogram(hist_dat,out_dir,label=name,x.lab,x.lim)
           }

# Reset x.lim if no x_lim argument entered
if(xlim==0)
    {x.lim=''}

print_histogram_abline(hist_dat,out_dir,label=name,x.lab,x.lim)

return(hist_dat)

}

get_histogram_all<-function(s.obj,out_dir,run_tag,xlab)
{
 sample.mat<-s.obj@assays$RNA@counts
 
 s_id<-s.obj@project.name
 name<-paste0(run_tag,s_id)
 
 if(xlab=='')
   {xlab=c('LibrarySize','GeneCounts','MtPerc')}

 for(j in 1:length(xlab))
   {
    if(xlab[i]=="LibrarySize")
     {
      #x_lim<-c(0,50000)
      counts.dat<-Matrix::colSums(sample.mat)
      print_histogram(counts.dat,out_dir,label=name,xlab)
     }else if(feature[j]=="gene")
        {xlab<-"GeneCounts"
         #x_lim<-c(0,2000)
         gene.dat<-Matrix::colSums(sample.mat>0)
         print_histogram(gene.dat,out_dir,label=name,xlab)
	}else if(feature[j]=="mito")
           {xlab<-"Mt_perc"
            #x_lim<-c(0,50)
            mt.dat<-s.obj$perc.mt
 	    print_histogram(mt.dat,out_dir,label=name,xlab)
           }
     }

hist_list<-list(counts.dat,gene.dat,mt.dat)

return(hist_list)

}

# Funtion to generate histograms of aggregate data #
get_agg_histograms<-function(hist_list,run_tag,out_dir)
{
agg_counts<-do.call(c,(lapply(hist_list, `[[`,1)))
agg_genes<-do.call(c,(lapply(hist_list, `[[`,2)))
agg_mt<-do.call(c,(lapply(hist_list, `[[` ,3)))

agg_list<-list(agg_counts,agg_genes,agg_mt)
x_lab<-c("LibrarySize","GeneCounts","Mt_perc")
#x_lim<-c(50000,10000,50)
for(i in 1:length(agg_list))
{
 f.name<-paste0(run_tag,'Aggregate_Histogram_',x_lab[i])
 print_histogram(agg_list[[i]],out_dir,x_lab[i],file.name=f.name)
}
 #rm(agg_list,agg_counts,agg_genes,agg_mt)
}

# Function to filter low quality cells #
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

# Note: The following function helps preserve R friendly data structure but generates very different results from Seurat subset() so decided not to use this!

# vars<-c('perc.mt','nCount_RNA','nFeature_RNA')
# low<-c(0 ,genecnt.thres[1],libsize.thres[1])
# high<-c(mt.thres, genecnt.thres[2], libsize.thres[2])
# for(i in length(vars))
#   {
#    expr<-FetchData(s.obj, vars=vars[i])
#    s.obj<-s.obj[,which(x=expr>low[i] & expr<high[i])]
#   }   

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

# Print ViolinPlots pre and post normalization


# Funtion to pre-process data #

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

# Scaling all genes
#if(verbose)
#{cat('Scaling all features ',sep='\n')}
#all.genes<-row.names(s.obj)
#s.obj<-ScaleData(s.obj,features=all.genes)

return(s.obj)

}

# Function to plot variable genes in samples # - from Jean Clemenceau
get_var_genes<-function(s.obj,out_dir='./',verbose=FALSE)
{
 
  f.name<-paste0(out_dir,s.obj@project.name,'-variableGenes.tsv')
  
  #Rank and sort variable genes
  var.gen.data<-HVFInfo(s.obj)
  var.gen.data<-var.gen.data[order(var.gen.data[,3],decreasing = TRUE),]
  
  #Export to tsv file
  fwrite(var.gen.data,f.name,sep='\t',eol ='\n',col.names = TRUE,row.names = TRUE)
}

#PCA plot function # - from Jean Clemenceau
plot_pca<-function(sample,cores=argv$cores,verbose=argv$verbose,plotTitle=argv$run_tag,theDims=argv$pca_dimensions){
  # Split PCs for rendering
  pca_split<-split(theDims, ceiling(seq_along(theDims)/6))

  # Visualize which genes are associated with reduced dimensions
  pca_gene_plots<-mclapply(pca_split,function(x){ VizDimLoadings(sample,dims=x,ncol=3,reduction='pca') },mc.cores=cores)

  #Rank PCs by % of variance explained by each, then find elbow
  pca_elbow_plot <-ElbowPlot(sample,ndims=50)+
    labs(title=plotTitle)
  
return(list(pca_gene_plots,pca_elbow_plot))
}



# Funtions to integrate multiple samples with batch correction #

cca_batch_correction<-function(s.obj.list,merged.title,genes='',reduction='cca',verbose=FALSE)
{
  if(verbose){
    cat("Performing CCA-MNN Pipeline",sep='\n')
  }
#Determine min number of neighbors
 mink<-min(200, min(sapply(seq_along(s.obj.list),function(x) ncol(s.obj.list[[x]]) ))  )

 sample.anchors<-FindIntegrationAnchors(s.obj.list,dims = 1:30,k.filter=mink,reduction=reduction)
 if(genes=='')
{if(verbose)
    {cat("Integration the sample.anchors only",sep='\n')}
 s.obj.integrated<-IntegrateData(anchorset=sample.anchors, dims=1:30)
}else if(genes=='all')
   {
    all.genes <- lapply(s.obj.list, row.names) %>% Reduce(intersect, .)
    
     if(verbose)
       {cat("CCA_MNN batch correction - integrating all genes", sep='\n')
        cat("Total genes being used for integration = ", length(all.genes),'\n')
       }
    s.obj.integrated<-IntegrateData(anchorset=sample.anchors, dims=1:30,features.to.integrate=all.genes)   
   }else{
      if(verbose)
       {cat("CCA_MNN batch correction - integrating top ", genes, " variable genes", sep='\n')}
      m.obj<-merge(s.obj.list[[1]],s.obj.list[2:length(s.obj.list)])
      m.obj<-FindVariableFeatures(m.obj, nfeatures=genes)
      
      int.genes<- m.obj@assays$RNA@var.features
      #int.genes<-union(sample.anchors,m.obj@assays$RNA@var.features)
      s.obj.integrated<-IntegrateData(anchorset=sample.anchors, dims=1:30,features.to.integrate=int.genes)
    rm(m.obj)
    }

#saveRDS(sample.anchors,"sample.anchors.rds")
rm(sample.anchors)

s.obj.integrated@project.name<-merged.title

 return(s.obj.integrated)
}

# Function to add metadata to Seurat object # - adapt function as per metafile
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

# Function to retrieve markers from file(s) #

get_diff_gene_list<-function(diff_dir,verbose=FALSE)
{
if(verbose)
{cat("Extracting diff genes from files",sep='\n')}
 setwd(diff_dir)
 f<-list.files()
 markers<-lapply(f,read.delim,header=TRUE)

 return(markers)

}

get_cluster_markers<-function(diff_dir,n=10)
{
  f.name<-list.files(diff_dir)
  f.name<-f.name[grep(".tsv",f.name)]
  cluster<-vector()
  markers<-vector()
  
  for(i in 1:length(f.name))
  {
    f<-read.delim(paste0(diff_dir,f.name[i]), header=TRUE)
    m<-f[(1:n),1]
    clus<-unlist(strsplit(unlist(strsplit(f.name[i], split='-'))[2],split='_'))[1]
    cluster<-append(cluster,rep(clus, times=n))
    markers<-append(markers,m)
  }
  
  cluster<-as.numeric(cluster)
  marker.table<-data.frame("Cluster"= cluster,"Marker"=markers)
  #write.csv(marker.table, paste0(out_dir,run_tag,"Markers.table.csv"),row.names=FALSE)
  return(marker.table)
}


### Plotting Functions ###

#Function to plot variable genes in samples # - from Jean Clemenceau

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

# Function to plot consistent UMAPs #

customized_umap<-function(s.obj,umap_cols=NULL,label=FALSE,title=NULL, group=NULL,split=NULL,dot=0.3,save=FALSE,out='./',run_tag,h=8,w=8)
{
    ns<-s.obj@meta.data$sample %>% unique %>% length()
    nc<-ncol(obj)

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

    umap<-DimPlot(s.obj, group.by=group, split.by=split,pt.size=dot,cols=umap_cols, label=label)

    r1<-range(s.obj@reductions$umap@cell.embeddings[,1])
    r2<-range(s.obj@reductions$umap@cell.embeddings[,2])

    sc<-c((floor(min(r1[1],r2[1]))),
          ceiling((max(r1[2],r2[2]))))

    umap<-umap+theme(legend.position = 'bottom')+
      scale_x_continuous(limits=sc)+
      scale_y_continuous(limits=sc)+
      ggtitle(label=title,
              subtitle = sub)
if(is.null(group)==FALSE)
{
    run_tag<-paste0(run_tag,group,'_')
}
if(save==TRUE)
{
    ggsave(paste0(out,run_tag,'UMAP.png'),
               width=w, height=h, units="in",umap)

}
    return(umap)
    
}

# Function to tabulate differentially expressed genes 
differential_gene_exp<-function(s.obj,clusters,out_dir='./',run_file_name="",verbose=FALSE,save=TRUE)

{if(verbose)

{cat("Finding Differentially expressed cluster markers", '\n')}
  
 DefaultAssay(s.obj)<-"RNA"
 markers<-list()
 
#Create directory to export markers
  if(!dir.exists(paste0(out_dir,'diff_genes/')))
        {dir.create( paste0(out_dir,'diff_genes/'))}
    out_path<-paste0(out_dir,'diff_genes/')

# Find markers and export to tsv file
  for(i in 1:length(clusters))
  {if( sum(s.obj@meta.data$seurat_clusters==clusters[i])<4)
    {
     cat('Cluster',clusters[i],'has too few cells','\n')
    }else{
   	  markers[[i]]<-FindMarkers(s.obj, ident.1= clusters[i])
    	  gene<-rownames(markers[[i]])
    	  markers[[i]]<-cbind(gene,markers[[i]])
    	  markers[[i]]<-markers[[i]] %>% filter(p_val_adj<0.05) %>% mutate(pct.diff=pct.1-pct.2)%>% arrange(desc(pct.diff,avg_logFC))
   #markers[[i]]<-markers[[i]][order(markers[[i]]["avg_logFC"], decreasing=TRUE),]
    	  if(save)
		{fwrite(markers[[i]], paste0(out_path,run_file_name,'cluster',clusters[i],'_diff_markers.tsv') ,append=FALSE,sep='\t',row.names=FALSE,col.names=TRUE )}  
	}
  } 
 return(markers)
}

# Function to tabulate differentially expressed genes - between 2 or more clusters
differential_gene_exp2<-function(ident1,ident2,s.obj,out_dir='./',run_tag="",verbose=FALSE)

{if(verbose)

{cat("Finding Differentially expressed cluster markers", '\n')}

 DefaultAssay(s.obj)<-"RNA" 

# Find markers and export to tsv file
          markers<-FindMarkers(s.obj, ident.1=ident1,ident.2=ident2)
          gene<-rownames(markers)
          markers<-cbind(gene,markers)
          markers<-markers %>% filter(p_val_adj<0.05) %>% mutate(pct.diff=pct.1-pct.2)%>% arrange(desc(pct.diff,avg_logFC))
   #markers[[i]]<-markers[[i]][order(markers[[i]]["avg_logFC"], decreasing=TRUE),]

if(length(ident1)>1)
{
 ident1<-paste0(ident1, collapse="_")
}

if(length(ident2)>1)
{
 ident2<-paste0(ident2,collapse="_")
}

          fwrite(markers,paste0(out_dir,run_tag,'cluster_',ident1,'_vs_cluster_',ident2,'_diff_markers.tsv') ,append=FALSE,sep='\t',row.names=FALSE,col.names=TRUE )
        
}

# Function to generate Features plot #

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

# Funtion to generate dendrogram #
plot_dendrogram<-function(s.obj,run_file_name,features, out_dir='./',verbose=FALSE)
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
png(paste0(out_dir,run_file_name,"dendrogram.png"),width=11,height = 8.5, units="in",res=300)
  hc_avg<-hclust(as.dist(1-cor(collapsed_norm_mat)),method = "average")
  plot(hc_avg, main=paste0(run_file_name,"Cluster Dendrogram"), xlab="")
  dev.off()

  # ggsave(paste0(run_file_name,"dendrogram.png"), path=out_dir,width=11,height = 8.5, units="in",d)
   # Select diff markers
  m<-which(rownames(collapsed_norm_mat) %in% features==TRUE)
  collapsed_norm_diff_markers<-collapsed_norm_mat[m,]

  options(bitmapType='cairo')
  png(paste0(out_dir,run_file_name,"diff_markers_dendrogram.png"),width=11,height = 8.5, units="in",res=300)
  hc2_avg<-hclust(as.dist(1-(cor(collapsed_norm_diff_markers))),method = "average")
  plot(hc2_avg)
  dev.off()
}

# Functions to implement clustree workflow # 
iterative_clus_by_res<-function(s.obj,res,dims_use,reduction='PCA', verbose=FALSE)
{ if(verbose)
   {cat("Performing iterative clustering by resolution for PCs 1:",max(dims_use),'\n')}

  s.obj<-FindNeighbors(s.obj,dims=dims_use,reduction=reduction)
  for(i in 1:length(res))
  {
    s.obj<-FindClusters(s.obj, res=res[i])
  }
  return(s.obj)
}

print_clustree_png<-function(s.obj,prefix,out_dir='./',run_tag,verbose=FALSE)
{ if(verbose)
   {cat("Generating clustree png(s)", sep='\n')}
   if(!dir.exists(paste0(out_dir,'clustree/')))
  {dir.create(paste0(out_dir,'clustree/'))}
 out_path<-paste0(out_dir,'clustree/')
 options(bitmapType='cairo')
 png(paste0(out_path,run_tag,"_clustree.png"), width=8.5, height=11, units="in",res=300)
  obj_clustree<-clustree(s.obj,prefix=prefix)
  print(obj_clustree)
  dev.off()
}

print_clustree_pdf<-function(s.obj,prefix,out_dir='./',run_tag, verbose=FALSE)
{ if(verbose)
   {cat("Generating clustree pdf(s)", sep='\n')}
   if(!dir.exists(paste0(out_dir,'clustree/')))
  {dir.create(paste0(out_dir,'clustree/'))}
  out_path<-paste0(out_dir,'clustree/')

  obj_clustree<-clustree(s.obj,prefix=prefix)
  ggsave(paste0(run_tag,"_clustree.pdf"),path=out_path, width=11,height=8.5,units="in",clus)
  #dev.off()
}

print_geneplots_on_clustree_RNA<-function(s.obj,genes, prefix, fun_use='median', out_dir, run_tag, verbose=FALSE)
{
 validated_genes<-genes[which(genes %in% rownames(s.obj)==TRUE)]
 if(length(validated_genes)>0)
{ if(verbose)
   {cat("Generating clustree geneplots",sep='\n')}
  
 if(!dir.exists(paste0(out_dir,'clustree_geneplots/RNA_assay')))
  {dir.create(paste0(out_dir,file.path("clustree_geneplots","RNA_assay")),recursive=TRUE)}
 out_path<-paste0(out_dir,'clustree_geneplots/RNA_assay/')
 
 DefaultAssay(s.obj)<-'RNA'
 
 for(i in 1:length(validated_genes))
 {
  clustree(s.obj, prefix=prefix,node_colour = validated_genes[i], node_colour_aggr = fun_use)
  ggsave(paste0(run_tag,'_RNAassay_',validated_genes[i],".png"),path=out_path, width=8.5, height=11,units="in")
 }
}else
  {cat("None of the requested gene(s) found!",'\n')}
}

print_geneplots_on_clustree_integrated<-function(s.obj,genes,  prefix = "integrated_snn_res.", fun_use='median', out_dir, run_tag, verbose=FALSE)
{
 validated_genes<-genes[which(genes %in% rownames(s.obj)==TRUE)]
 if(length(validated_genes)>0)
{ if(verbose)
   {cat("Generating clustree geneplots",sep='\n')}
  
 if(!dir.exists(paste0(out_dir,'clustree_geneplots/integrated_assay')))
  {dir.create(paste0(out_dir,file.path("clustree_geneplots","integrated_assay")),recursive=TRUE)}
 out_path<-paste0(out_dir,'clustree_geneplots/integrated_assay/')
 
 DefaultAssay(s.obj)<-'integrated'
 
 for(i in 1:length(validated_genes))
 {
  clustree(s.obj, prefix=prefix ,node_colour = validated_genes[i], node_colour_aggr = fun_use)
  ggsave(paste0(run_tag,'_',validated_genes[i],".png"),path=out_path, width=8.5, height=11,units="in")
 }
}else
  {cat("None of the requested gene(s) found!",'\n')}
}

# Function to generate Silhouette plots #
get_silhouette_plot<-function(s.obj,reduction='pca',dims=1:50,out_dir='./',run_tag,verbose=FALSE)
{
 if(verbose)
  {cat('Generating Silhouette plots','\n')} 
 if(!dir.exists(paste0(out_dir,'Silhouette/')))
  {dir.create(paste0(out_dir,'Silhouette/'))}
 out_path<-paste0(out_dir,'Silhouette/')
  clusters<-s.obj$seurat_clusters
  dist.matrix<-dist(x=Embeddings(object=s.obj[[reduction]])[,dims])
  sil<-silhouette(x=as.numeric(as.factor(clusters)),dist=dist.matrix)
  
  #Printing the Silhouette plot
  fviz_silhouette(sil)
  ggsave(paste0(run_tag,"_silhouette.png"),path=out_path, width=33,height=10)
}
