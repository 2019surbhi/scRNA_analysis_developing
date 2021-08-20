#!/usr/bin/env Rscript

# This is a script that facilitates running the Ting Lab scRNA seq pipeline on HPC for processing data generated on 10XGenomics platform 

source('./scRNA_pipeline_functions_final_HPC.R')

library(argparser)

### User inputs and parameters ###

parser<-arg_parser(name="scRNA_pipepine_run_final_HPC.R",description="Ting lab Seurat pipeline for scRNA and subset analysis")

parser<-add_argument(
  parser,
  arg='--input_dir',
  short = '-i',
  type="character",
  default='./',
  help="Enter input directory path in the format /path/to/input_dir/ Default=./ ")

parser<-add_argument(
  parser,
  arg='--samples',
  short = '-s',
  type="character",
  default='all',
  help="Enter sample names separated by(, for samples from same flowcell and : for samples from different flowcell)  Default is 'all' samples within the specified input directory")

parser<-add_argument(
  parser,
  arg='--data_dir',
  short = '-z',
  flag=TRUE,
  help="Flag set to read cellranger output from directory instead of hdf5 file")

parser<-add_argument(
  parser,
  arg='--object',
  short='-b',
  type="character",
  default='',
  help="Enter Seurat object with full path - for subset analysis")

parser<-add_argument(
  parser,
  arg='--clusters',
  short='-l',
  type="character",
  default='all',
  help="Enter which clusters (separated by ,) to be used for subset analysis. Use this argument only when performing subset analysis")

parser<-add_argument(
  parser,
  arg='--out_dir',
  short = '-o',
  type="character",
  default='./',
  help="Enter output directory path. Format: path/to/output/dir/ Default= ./")

parser<-add_argument(
  parser,
  arg='--run_tag',
  short = '-r',
  type="character",
  default='',
  help="Enter prefix to mark files for current run")

parser<-add_argument(
  parser,
  arg='--meta_file',
  short = '-f',
  type="character",
  default='',
  help="Enter the meta data file in csv format along with the path")

parser<-add_argument(
  parser,
  arg='--thresholds',
  short = '-t',
  type="character",
  default='',
  help="Enter QC thresholds for mt%, gene count (lower and upper cutoffs), and library size cutoff (lower and upper cutiffs) Format: <% MT genes>,<gene count lower cutoff>,<gene count upper cutoff>,<libsize low cutoff>,<libsize high cutoff>'")

parser<-add_argument(
  parser,
  arg='--hvg',
  short = '-y',
  default=2000,
  type='numeric',
  help="Enter the number of highly variable genes to be used for clustering. Defulat = 2000")

parser<-add_argument(
  parser,
  arg='--pca_dimensions',
  short = '-d',
  type="character",
  default='',
  help="Set of PCA dimensions (components) to be applied for clustering and dimensionality reduction. Format:'<pc1>:<pc_end>' e.g. '1:25' to select first 25 PCs OR '<pc1>,<pc2>,<pc3>' e.g. '1,2,3' to select first 3 PCs. Default=50 PCs")

parser<-add_argument(
  parser,
  arg='--cluster_resolution',
  short = '-e',
  type='numeric',
  default=0.5,
  help="Enter clustering resolution. Default =0.5. Needs to be optimized for each dataset [a helpful tool: clustree]")

parser<-add_argument(
  parser,
  arg='--batch_genes',
  short='-g',
  default='',
  type='character',
  help="Enter number of genes to integrate: Leave blank to integrate sample.anchors, enter 'all' to integrate all genes, or a numeric value to  enter the number of genes to integrate - the top variable genes calculated on merged object will be used for integration for this last option")

parser<-add_argument(
 parser,
 arg='--qc_plots',
 short='-q',
 flag=TRUE,
 help="Set flag if you wish to generate qc plots (pre and post filtering)")

parser<-add_argument(
  parser,
  arg='--clustering_optimization',
  short = '-u',
  type="character",
  default='none',
  help="specify which clustering optimizations to run, Enter 'clustree' for running clustree or/and 'sil' for running silhouette'. Separate multiple entries with ','. Default is set to run neither")

parser<-add_argument(
parser,
arg='--gene_list',
short = '-n',
type='character',
default="KRT5,KRT14,KRT15,KRT17,KRT20,UPK2,UPK3A,UPK3B,UPK1A,UPK1B,EPCAM,PTPRC,CD8A,CD79A,CD79B,LST1,TRAC,NKG7,KLRD1,TPSAB,PECAM1,FLT1,DCN,POSTN,ACTA2,MHY11,MKI67,TP63",
help="Enter the list of genes to be plotted on clustree [separated by ,]")

parser<-add_argument(
 parser,
 arg='--save',
 short='-a',
 type="character",
 default='both',
 help="Specify choices for saving Seurat object. Choices: 'none' , 'integrated' (to save object prior to clustering), 'final' (to save object post clustering), or 'both' (default) ")

parser<-add_argument(
  parser,
  arg='--cores',
  short = '-c',
  type="numeric",
  default=4,
  help= "Number of cores to use when processing samples in parallel.")

parser<-add_argument(
  parser,
  arg='--mem',
  short = '-m',
  type="numeric",
  default=40,
  help= "Memory allocated for gloabl objects in GB. Default is 40")

parser<-add_argument(
  parser,
  arg='--verbose',
  short = '-v',
  flag=TRUE,
  help="Set flag for printing output messages during the run")

argv <- parse_args(parser)

#Print run parameters
cat("List of arguments:",'\n')
argv
cat('\n')

### Process arguments ###

# Prepare output directory
if(! dir.exists(argv$out_dir) ){
  dir.create(argv$out_dir)
}

#Save run parameters - find better ways to save this in txt format
write.table(argv,paste0(argv$out_dir,"run_parameters.txt"),sep='\t',col.names=FALSE)

# Add separater to run name
if(argv$run_tag!=''){
  argv$run_tag<- gsub(' ','_',argv$run_tag)
 }
argv$run_tag<-paste0(argv$run_tag,'_')

#Process thresholds
  if(argv$thresholds!='')
   {
    argv$thresholds= as.numeric(unlist(strsplit(argv$thresholds,split = ',',fixed = TRUE)))
   }

#Process PCA components
if(argv$pca_dimensions==''){
  argv$pca_dimensions<- c(1:50)
}else{
  dim.range<- as.numeric(unlist(strsplit(argv$pca_dimensions,split = ':',fixed = TRUE)))
 if(length(dim.range)>1) 
    {
     argv$pca_dimensions=as.numeric(dim.range[1]):as.numeric(dim.range[2])
    }else
     {
      argv$pca_dimensions<-as.numeric(unlist(strsplit(argv$pca_dimensions, split=',', fixed=TRUE)))
     }  
}

# Process batch genes arguments
if(argv$batch_genes!='')
  {if(argv$batch_genes!='all')
	{argv$batch_genes<-as.numeric(argv$batch_genes)}
  }

# Process cluster numbers
if(argv$clusters!='all')
{
argv$clusters=as.numeric(unlist(strsplit(argv$clusters,split = ',',fixed = TRUE)))

}

# Process clustering optimization arguments
argv$clustering_optimization<-unlist(strsplit(argv$clustering_optimization, split=','))

clustree<-ifelse('clustree' %in% argv$clustering_optimization,TRUE,FALSE)
sil<-ifelse('sil' %in% argv$clustering_optimization,TRUE,FALSE)
ikap<-ifelse('ikap' %in% argv$clustering_optimization,TRUE,FALSE)

#Generate a list of genes for clustree from input
argv$gene_list<-unlist(strsplit(argv$gene_list,split=','))

# Print processed arguments
cat("List of processed arguments:",'\n')
argv
cat('\n')

#Save arguments
saveRDS(argv,paste0(argv$out_dir,argv$run_tag,"arguments.rds"))


#Facilitate parallel processing

plan("multicore", workers=argv$cores)
options(future.globals.maxSize=((argv$mem*1000)*1024^2))
options(future.rng.onMisue = "ignore")
#n<-length(argv$samples)

setwd(argv$out_dir)

obj.list<-list()

#Load object from object_list (if added)
if(argv$object=='')
{

# Determine if the samples are to be loaded from single or multiple directories

paths.list<-unlist(strsplit(argv$input_dir,split=':'))

if(length(paths.list)>1)
{
 if(argv$verbose)
  {cat("Loading data from multiple datasets/paths" ,'\n')}

 s.list<-list()

if(argv$samples=='all')
 {if(argv$verbose)
  {cat('You did not enter any samples, Running pipeline on all samples from the specified input directories',sep='\n')}
   s<-lapply(paths.list,list.files)
   }else{
      #Create a list of samples corresponding to each path
      s<-strsplit(argv$samples,split=':') %>% unlist %>% lapply(FUN=function(x) {strsplit(x,split=',') %>%  unlist})
         }
# Create a list of seurat object
for(i in 1:length(paths.list))
 {
  n<-length(s[[i]])

   s.list[[i]]<-lapply(s[[i]][1:n], create_seurat_obj_10X,input_dir=paths.list[i],argv$data_dir,verbose=argv$verbose)

  }

obj.list<-unlist(s.list)
rm(s.list)
rm(s)

}else{
  # Extract samples list
   if(argv$verbose)
  {cat('Loading data from single path', '\n')}
   if(argv$samples=='all')
     {if(argv$verbose)
       {cat('you did not enter any samples, Running pipeline on all samples in the input directory',sep='\n')}
      argv$samples<-list.files(argv$input_dir)
     }else{
           argv$samples<-unlist(strsplit(argv$samples,split = ',',fixed = TRUE))
          }
  n<-length(argv$samples)
  #Create Seurat object from single directory
   obj.list<-lapply(argv$samples[1:n], create_seurat_obj_10X,input_dir=argv$input_dir,argv$data_dir,verbose=argv$verbose)

 }
}else{
  if(argv$verbose)
{cat("Begin sub-clustering",'\n')}

#Load Seurat object (containing all samples)
obj<-readRDS(argv$object)

if(argv$clusters=='all')
{sub<-obj}else{
		#Subset Seurat object to include only specific clusters defined by user
		sub<-subset(obj,idents=argv$clusters)
	 	}

sub.list<-list()
sub.counts<-list()

#Split the subsetted object by sample
if(argv$verbose)
{ cat("Splitting object by sample ID",'\n')}

sub.list<-SplitObject(sub, split.by = "orig.ident")
sample.id<-names(sub.list)

#Extract counts data
sub.counts<-lapply(X=1:length(sub.list), function(x) {return(sub.list[[x]]@assays$RNA@counts)})

#Create Seurat object
obj.list<-mapply(create_seurat_obj_from_counts_data,sub.counts,sample.id,verbose=argv$verbose)

rm(sub)
rm(sub.list)
rm(sub.counts)
rm(sample.id)

}

if(argv$qc_plots)
{

#Get pre-filter histogram [sample-wise + aggregate]
if(!dir.exists(paste0(argv$out_dir,'qc_plot/histograms/raw')))
        {dir.create(file.path("qc_plots","histograms","raw"), recursive = TRUE)}
    out_path<-paste0(argv$out_dir,'qc_plots/histograms/raw/')

x_lab<-c('LibrarySize','GeneCounts','MtPerc')
#x_lim<-c(120000,10000,100)
x_lim<-c(50000,2000,50)

cutoff_low<-c(argv$thresholds[4],argv$thresholds[2],argv$thresholds[1])
cutoff_high<-c(argv$thresholds[4],argv$thresholds[3],argv$thresholds[1])

for(i in 1:length(x_lab))
{

 pre<-paste0(argv$run_tag,'Histogram_')
 # Get Histogram for each sample
 hist_list<-lapply(obj.list,get_histogram,out_path,prefix=pre,x_lab[i])

 #Get aggregate Histogram
 lab<-paste0(argv$run_tag,'AggregateHistogram')
 agg<-unlist(hist_list)
 print_histogram_abline(agg,out_path,label=lab,x_lab[i],cutoff=c(cutoff_low[i],cutoff_high[i]))

 #Get aggregate Histogram low xlim
 lab<-paste0(argv$run_tag,'AggregateHistogram_low')
 print_histogram_abline(agg,out_path,label=lab,x_lab[i],cutoff=c(cutoff_low[i],cutoff_high[i]),x=c(0,x_lim[i]))

}

rm(hist_list)
rm(agg)
rm(out_path)
rm(lab)
rm(x_lab)
rm(x_lim)

#hist_list<-lapply(X=1:length(obj.list),FUN=function(x){get_histogram(obj.list[[x]],out_path,argv$run_tag)})
#get_agg_histograms(hist_list,argv$run_tag,out_path)

}

if(length(argv$thresholds)>1)
{
# Get cell table
tab_list<-lapply(obj.list[1:length(obj.list)],get_cell_table,mt.thres=argv$thresholds[1], genecnt.thres=argv$thresholds[2:3], libsize.thres=argv$thresholds[4:5], verbose=argv$verbose)
tab<-data.frame()
for(i in 1:length(tab_list))
{
tab<-rbind(tab,tab_list[[i]])
}

total<-colSums(tab)
tab<-rbind(tab,"Total"=total)

write.csv(tab,paste0(argv$out_dir,argv$run_tag,"cell_qc_table.csv"))
rm(tab_list)
rm(tab)

# Filter low quality cells
obj_tab_list<-lapply(X=(1:length(obj.list)), FUN=function(x){filter_cells(obj.list[[x]],mt.thres=argv$thresholds[1], genecnt.thres=argv$thresholds[2:3], libsize.thres=argv$thresholds[4:5], verbose=argv$verbose)})

obj.list<-do.call(c,(lapply(obj_tab_list,`[[`,1)))
cell.data<-do.call(rbind,(lapply(obj_tab_list, `[[`, 2)))
#cell.data2<-adorn_totals(cell.data,"col")

write.csv(cell.data,paste0(argv$out_dir,argv$run_tag,"cell_threhold_table.csv"))

rm(obj_tab_list)
#cell.data<-do.call(c,lapply(cell_tab_list))

# Get threhold plots
cell.data$rank[order(cell.data$thr.cell.cnt,decreasing = TRUE)]<- c(1:nrow(cell.data))
if(argv$verbose){
    cat("Exporting cell count plots",sep='\n')}
threshold_plots<-plot_threshold_effects(cell.data,argv$thresholds,argv$run_tag)
pdf(file=paste0(argv$out_dir,argv$run_tag,"thresholds.pdf"),paper='a4',width=8)
lapply(threshold_plots,print)
dev.off()

rm(threshold_plots)
rm(cell.data)
}

if(argv$qc_plots)
{
# Get post threshold histograms
if(!dir.exists(paste0(argv$out_dir,'qc_plot/histograms/post_filter')))
        {dir.create(file.path("qc_plots","histograms","post_filter"), recursive = TRUE)}
    out_path<-paste0(argv$out_dir,'qc_plots/histograms/post_filter/')

x_lab<-c('LibrarySize','GeneCounts','MtPerc')
x_lim<-c(50000,2000,50)

for(i in 1:length(x_lab))
{
 pre<-paste0(argv$run_tag,'post-filter_Histogram_')
 # Get Histogram for each sample
 hist_list<-lapply(obj.list,get_histogram,out_path,prefix=pre,x_lab[i])

 #Get aggregate Histogram
 lab<-paste0(argv$run_tag,'post-filter_AggregateHistogram')
 agg<-unlist(hist_list)
 print_histogram(agg,out_path,label=lab,x_lab[i])

 #Get aggregate Histogram low xlim
 lab<-paste0(argv$run_tag,'post-filter_AggregateHistogram_low')
 print_histogram(agg,out_path,label=lab,x_lab[i],x=c(0,x_lim[i]))

}

rm(hist_list)
rm(pre)
rm(lab)
rm(agg)
rm(x_lab)
rm(x_lim)

#lapply(X=1:length(obj.list),FUN=function(x){get_histogram(obj.list[[x]],out_path,paste0(argv$run_tag,'post_filter_'))})

}

### Data pre-processing ###

obj.list<-lapply(obj.list[1:length(obj.list)],pre_process, hvg=argv$hvg, verbose=argv$verbose)

#Get variable genes table and plots for each sample
# Prepare output file
  if(!dir.exists(paste0(argv$out_dir,'variable_genes/'))){
    dir.create( paste0(argv$out_dir,'variable_genes/'))
  }
out_path<-paste0(argv$out_dir,'variable_genes/')
lapply(obj.list[1:length(obj.list)],get_var_genes,out_dir=out_path,verbose=argv$vebose)

 varplots<-lapply(obj.list[1:length(obj.list)],get_var_genes_plot,verbose=argv$verbose)
 pdf(file=paste0(out_path,argv$run_tag,'VariableGenePlots.pdf'),paper='a4')
 print(varplots)
 dev.off()

rm(out_path)
rm(varplots)

#saveRDS(obj.list, paste0(argv$out_dir,argv$run_tag,'obj_list.rds'))

if(length(obj.list)>1)
 {
 # Integrate data with batch correction
 
obj.integrated<-cca_batch_correction(obj.list,merged.title=argv$run_tag,genes=argv$batch_genes,verbose=argv$verbose)

}else{
   obj.integrated=obj.list[[1]]
   }

### Add meta data ###

if(argv$meta_file!='')
{
obj.integrated<-add_metadata(obj.integrated,argv$meta_file)
}

rm(obj.list)

if((argv$save=='integrated')||(argv$save=='both'))
 {
  #Saving object pre-clustering for clustering optimization (if needed)   
  saveRDS(obj.integrated,file=paste0(argv$out_dir,argv$run_tag,"integrated_only.rds"))
 }

### Pre-clustering ###

# Scale data

  #obj.integrated<-FindVariableFeatures(obj.integrated)
  all.features<-rownames(obj.integrated)
  obj.integrated<-ScaleData(obj.integrated,features=all.features, verbose=argv$verbose)
 
# Run PCA

  obj.integrated<-RunPCA(obj.integrated, npcs=50,ndims.print = 1:15, verbose=argv$verbose)
  obj.integrated[['sample']]<-obj.integrated[['orig.ident']]


### PCA Plots ###

# Elbow Plot

options(bitmapType='cairo')
png(file=paste0(argv$out_dir,argv$run_tag,"PCA_elbow_plots.png"),width = 11,height = 8.5, units='in', res=300)
  ElbowPlot(obj.integrated,ndims=50)+labs(paste0(argv$run_tag,'ElbowPlot'))
  dev.off()

# PCA gene plot (print 4 PCs per page)

  pc_genes_plot_list<-lapply(argv$pca_dimensions,function(x){ VizDimLoadings(obj.integrated,dims=x,ncol=1,reduction='pca') })

pca_plots<-marrangeGrob(pc_genes_plot_list, nrow=2, ncol=2)
ggsave(paste0(argv$out_dir,argv$run_tag,"PC_gene_plots.pdf"), width=8.5, height=11, units = "in", pca_plots)


### Clustering Optimization ###

# This section is optional to generate a set of plots for clustering optimization

if(clustree)
{
# Generating clustree and silhouette plots using batch integrated object
res<-seq(0.1,1.2,by=0.1)
pc<-c(15,20,25,30,35,40,50)

for(i in 1:length(pc))
{
#dims<-seq(1,pc[i],by=1)
obj_clustree<-NULL
clus_run=paste0(argv$run_tag,'PC',pc[i])
obj_clustree<-iterative_clus_by_res(obj.integrated, res=res,dims_use=1:pc[i],verbose=argv$verbose)
print_clustree_png(obj_clustree,prefix="integrated_snn_res.",out_dir=argv$out_dir,run_tag=clus_run,verbose=argv$verbose)

print_geneplots_on_clustree_integrated(obj_clustree,prefix="integrated_snn_res.", genes=argv$gene_list, fun_use='median',out_dir= argv$out_dir, run_tag=clus_run, verbose=FALSE)

 print_geneplots_on_clustree_RNA(obj_clustree,genes=argv$gene_list, fun_use='median',prefix='integrated_snn_res.', out_dir=argv$out_dir,run_tag=clus_run, verbose=FALSE)

}

rm(res)
rm(pc)
rm(clus_run)
rm(obj_clustree)

}

if(sil)
{

res<-seq(0.1,1.2,by=0.1)
pc<-c(15,20,25,30,35,40,50)

DefaultAssay(obj.integrated)<-'integrated'

for(i in 1:length(pc))
{
for(j in 1:length(res))
   {
   sil_run<-paste0(argv$run_tag,'PC',pc[i],'_res',res[j])
   obj_sil<-iterative_clus_by_res(obj.integrated,res=res[j],dims_use=1:pc[i],verbose=argv$verbose)
    get_silhouette_plot(obj_sil,reduction='pca',dims=1:pc[i],out_dir=argv$out_dir,run_tag=sil_run,verbose=argv$verbose)
   }
}

rm(pc)
rm(res)
rm(sil_run)
rm(obj_sil)

}

DefaultAssay(obj.integrated)<-'integrated'

### CLUSTERING ###

# In the pipeline the clusters are generated by default using 50 PCs at resolution 0.5

#Find Clusters
  obj.integrated<-FindNeighbors(obj.integrated, dims=argv$pca_dimensions)
  obj.integrated<-FindClusters(obj.integrated, res=as.numeric(argv$cluster_resolution))

#Generate UMAP
  if(argv$verbose)
  {cat("Genetrating UMAP plots", sep='\n')}
  obj.integrated<-RunUMAP(obj.integrated, reduction="pca", dims=argv$pca_dimensions)

#Print UMAPs

if(!dir.exists(paste0(argv$out_dir,'umaps/')))
        {dir.create( paste0(argv$out_dir,'umaps/'))}
    out_dir<-paste0(argv$out_dir,'umaps/')

customized_umap(obj.integrated, umap_cols=NULL,label=TRUE,title=NULL, group=NULL,split=NULL,dot=0.3,save=TRUE,out=out_dir,run_tag=argv$run_tag)
customized_umap(obj.integrated, umap_cols=NULL,label=TRUE,title=NULL, group='orig.ident',split=NULL,dot=0.3,save=TRUE,out=out_dir,run_tag=argv$run_tag)

if(argv$meta_file!='')
{
 customized_umap(obj.integrated, umap_cols=NULL,label=TRUE,title=NULL, group='sex',split=NULL,dot=0.3,save=TRUE,out=out_dir,run_tag=argv$run_tag)
 customized_umap(obj.integrated, umap_cols=NULL,label=TRUE,title=NULL, group='age',split=NULL,dot=0.3,save=TRUE,out=out_dir,run_tag=argv$run_tag)
 customized_umap(obj.integrated, umap_cols=NULL,label=TRUE,title=NULL, group='normal_tumor',split=NULL,dot=0.3,save=TRUE,out=out_dir,run_tag=argv$run_tag)
}
 
### Export seurat objects ###

# obj.integrated is an object with normalized batch corrected counts and allows clustering with optimal parameters while obj.clustered contains clustering info and can be used for further downstream analyses

if((argv$save=='final')||(argv$save=='both'))
{
 
 if(argv$verbose)
  {cat("Exporting Seurat integrated object",sep='\n')}
 
 saveRDS(obj.integrated,file=paste0(argv$out_dir,argv$run_tag,"clustered.rds"))
}

### Cell proportion table ###

tab<-table(Idents(obj.integrated),obj.integrated@meta.data$orig.ident)

#Convert to data frame
tab<-as.data.frame.matrix(tab)

#total<-table(obj.integrated@meta.data$orig.ident)
#cell.tab<-rbind(tab,total)

#Add totals

tab<-adorn_totals(tab, c("row","col"))

write.csv (tab,paste0(argv$out_dir,argv$run_tag,"cells_by_cluster_by_sample.csv"))

rm(tab)

### Differential gene expression ###

clusters_num<-levels(obj.integrated$seurat_clusters)

marker.list<-differential_gene_exp(obj.integrated, clusters=clusters_num, out_dir=argv$out_dir,run_file_name=argv$run_tag,verbose=argv$verbose)

### Feature Plots ###

if(!dir.exists(paste0(argv$out_dir,'feature_plots/')))
        {dir.create( paste0(argv$out_dir,'feature_plots/'))}
    out<-paste0(argv$out_dir,'feature_plots/')

# Extract feature list from marker list file

features.list<-list()

for(i in 1:length(marker.list))
{
 features.list[[i]]<-as.vector(marker.list[[i]][,'gene'])
}

# Generating pdf Plots

if(argv$verbose)
        {cat("Creating Feature plots (pdf format)", sep='\n')}

plot.list<-list()

for(i in 1:length(marker.list))
{
plot.list<-get_features_plot(features.list[[i]],obj.integrated, top=20)
p1<-marrangeGrob(plot.list, nrow=3, ncol=3)
ggsave(paste0(out,argv$run_tag,"FeaturePlot_cluster",clusters_num[i],".pdf"),width=11, height=8.5, units = "in", p1)
dev.off()
}

# Saving plots in png format as well
if(argv$verbose)
        {cat("Creating Feature plots (png format)", sep='\n')}

 DefaultAssay(obj.integrated)<-'RNA'
 n<-9 # Num of features to print per cluster

for(i in 1:length(features.list))
{

#png(paste0(out,run_file_name,"FeaturePlot_cluster",clusters_num[i],".png"),width=11, height=8.5, units = "in",res=600)
features<-features.list[[i]]
p2<-FeaturePlot(obj.integrated,features[1:n],min.cutoff='q30')
ggsave(paste0(out,argv$run_tag,"FeaturePlot_cluster",clusters_num[i],".png"),width=11, height=8.5, units = "in",p2)
}

rm(clusters_num)
rm(marker.list)
rm(plot.list)
rm(p1)
rm(p2)
rm(out)

### Heatmaps ###

# Generating heatmap for top 15 markers per cluster

DefaultAssay(obj.integrated)<-'integrated'
heatmap_features<-vector()
top_h<-15
for(i in 1:length(features.list))
{

heatmap_features<-append(heatmap_features,features.list[[i]][1:top_h])

}

rm(features.list)

#png(paste0(run_file_name,"Heatmap.png"),width=55,height=34,units="in",res=300)
h1<-DoHeatmap(obj.integrated,features=heatmap_features)
ggsave(paste0(argv$run_tag,"Heatmap_integrated.pdf"),path=argv$out_dir,width=22,height=17,units="in",h1)

rm(h1)
 
# Plot Heatmap for RNA assay

DefaultAssay(obj.integrated)<-'RNA'
obj.integrated<-ScaleData(obj.integrated, features=rownames(obj.integrated))
h2<-DoHeatmap(obj.integrated,features=heatmap_features)
ggsave(paste0(argv$run_tag,"Heatmap_RNA.pdf"),path=argv$out_dir,width=22,height=17,units="in",h2)

rm(h2)

### Plot Dendrogram ###

plot_dendrogram(obj.integrated,argv$run_tag,features=heatmap_features,out_dir=argv$out_dir)

rm(heatmap_features)
rm(obj.integrated)

sessionInfo()

cat('\n')

cat('End of script ','\n')

