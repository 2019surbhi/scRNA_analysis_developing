#!/usr/bin/env Rscript


# Developers: Surbhi Sona and Jean Clemenceu

source('/home/sonas/beegfs/scripts/R/tinglab_scRNA_pipeline_functions.R')

library(argparser)

### User inputs and parameters ###

parser<-arg_parser(name="tinglab_scRNA_pipeline.R",description="Ting lab Seurat pipeline for scRNA and subset analysis")

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
  help="Enter which clusters (separated by ,) to be used for subset analysis. Use this argument only when performing subset analysis to specify which clusters to use to subset the specified Seurat obj. Specify this option as 'none' to skip subset and batch correction")

parser<-add_argument(
  parser,
  arg='--output_dir',
  short = '-o',
  type="character",
  default='./',
  help="Enter output directory path. Format: path/to/output/dir/ Default= ./")

parser<-add_argument(
  parser,
  arg='--file_prefix',
  short = '-f',
  type="character",
  default='',
  help="Enter file prefix to mark files for current run")

parser<-add_argument(
  parser,
  arg='--meta_file',
  short = '-F',
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
  arg='--batch_correction',
  short='-C',
  default='none',
  type='character',
  help="Enter the batch correction method to be used. Options are: harmony and cca-mnn. Default is none and if not changed the rest of the piepline will run on the single input object without batch correction")

parser<-add_argument(
  parser,
  arg='--batch_genes',
  short='-g',
  default='',
  type='character',
  help="Enter number of genes to integrate: Leave blank to integrate sample.anchors, enter 'all' to integrate all genes, or a numeric value to  enter the number of genes to integrate (should be same as number of variable genes computed) ")

parser<-add_argument(
 parser,
 arg='--qc_only',
 short='-Q',
 flag=TRUE,
 help="Set flag to run only QC")

parser<-add_argument(
 parser,
 arg='--qc_plots',
 short='-q',
 flag=TRUE,
 help="Set flag to generate qc plots (pre and post filtering)")

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
 help="Specify choices for saving Seurat object. Choices: 'none' , 'integrated' (to save batch corrected object prior to clustering), 'clustered' (to save object post clustering), or 'both' (default) ")

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

args <- parse_args(parser)

#Print run parameters
cat("List of arguments:",'\n')
args
cat('\n')

### Process arguments ###

# Prepare output directory
if(! dir.exists(args$output_dir))
{
  dir.create(args$output_dir)
}

#Save run parameters - find better ways to save this in txt format
user_input<-unlist(args)
names(user_input)<-NULL
args_tab<-cbind(names(args),user_input)
colnames(args_tab)<-c('arguments','user_input')
write.csv(args_tab,paste0(args$output_dir,args$file_prefix,"run_parameters.csv"))

# Add separator to file name
if(args$file_prefix!='')
{
  args$file_prefix<- gsub(' ','_',args$file_prefix)
 }

args$file_prefix<-paste0(args$file_prefix,'_')

#Process cell filtering thresholds
  if(args$thresholds!='')
   {
    args$thresholds= as.numeric(unlist(strsplit(args$thresholds,split = ',',fixed = TRUE)))
   }

#Process PCA components

if(args$pca_dimensions=='')
{# If user didn't choose PCs then use all 50 PCs
  args$pca_dimensions<- c(1:50)
}else
  {
   # Create a vector of user specified PCs
   dim.range<- as.numeric(unlist(strsplit(args$pca_dimensions,split = ':',fixed = TRUE)))
   if(length(dim.range)>1)
    {
     # Create PCs vector if user specified range of PCs (e.g. 1:15)
     args$pca_dimensions=as.numeric(dim.range[1]):as.numeric(dim.range[2])
    }else
      {
       # Create PCs vector if user entered specific sets of PCs (e.g. 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
       args$pca_dimensions<-as.numeric(unlist(strsplit(args$pca_dimensions, split=',', fixed=TRUE)))
      }
  }

# Process resolution
args$cluster_resolution<-as.numeric(args$cluster_resolution)

# Process batch genes arguments
if(args$batch_genes!='')
  {
   if(args$batch_genes!='all')
	{
     args$batch_genes<-as.numeric(args$batch_genes)
    }
  }

# Process cluster numbers
if(args$clusters!='all')
{
 args$clusters=as.numeric(unlist(strsplit(args$clusters,split = ',',fixed = TRUE)))
}

# Process clustering optimization arguments
args$clustering_optimization<-unlist(strsplit(args$clustering_optimization, split=','))

clustree<-ifelse('clustree' %in% args$clustering_optimization,TRUE,FALSE)
sil<-ifelse('sil' %in% args$clustering_optimization,TRUE,FALSE)

#Generate a list of genes for clustree from input
args$gene_list<-unlist(strsplit(args$gene_list,split=','))

# Print processed arguments
cat("List of processed arguments:",'\n')
args
cat('\n')

#Save arguments
saveRDS(args,paste0(args$output_dir,args$file_prefix,"arguments.rds"))


###(1) Load data or prepare subset analysis data ###

#Facilitate parallel processing
plan("multicore", workers=args$cores)
options(future.globals.maxSize=((args$mem*1000)*1024^2))
options(future.rng.onMisue = "ignore")
#n<-length(args$samples)

setwd(args$output_dir)

#Initialize skip to FALSE to avoid skipping batch correction
skip<-FALSE
obj.list<-list()

# [if-else block 1]: If no object specified, then the script looks for expression matrix in the input directory
if(args$object=='')
{
 # Determine if the samples are to be loaded from single or multiple directories
 paths.list<-unlist(strsplit(args$input_dir,split=':'))
  
  ## [if-else block 1a]: If multiple paths are specified by the user
  if(length(paths.list)>1)
   {
    if(args$verbose)
      {cat("Loading data from multiple datasets/paths" ,'\n')}

    s.list<-list()
       
    ### [if-else block 1b]: If user didn't specify sample names, all samples from input directory will be loaded
    if(args$samples=='all')
     {if(args$verbose)
        {cat('You did not enter any samples, Running pipeline on all samples from the specified input directories',sep='\n')}
      s<-lapply(paths.list,list.files)
     }else ### [if-else block 1b]: If user specified sample names
       {
        # Create a list of samples corresponding to each path
        s<-strsplit(args$samples,split=':') %>% unlist %>% lapply(FUN=function(x) {strsplit(x,split=',') %>%  unlist})
       }
       
    # Now Create a list of seurat object [all samples or user specified samples]
    for(i in 1:length(paths.list))
    {
     n<-length(s[[i]])

     s.list[[i]]<-lapply(s[[i]][1:n], create_seurat_obj_10X,input_dir=paths.list[i],args$data_dir,verbose=args$verbose)
    }

   obj.list<-unlist(s.list)
   rm(s.list)
   rm(s)

  }else
    {## [if-else block 1a]: If single path specified

     # Extract samples list
     if(args$verbose)
        {cat('Loading data from single path', '\n')}
    
     ### [if-else block 1c]: If user didn't specify sample names, all samples from single input directory will be loaded
     if(args$samples=='all')
        {
         if(args$verbose)
          {cat('you did not enter any samples, Running pipeline on all samples in the input directory',sep='\n')}
         args$samples<-list.files(args$input_dir)
        }else ### [if-else block 1c]: User specified samples from single input directory will be loaded
          {
           args$samples<-unlist(strsplit(args$samples,split = ',',fixed = TRUE))
          }
     
    n<-length(args$samples)
    #Create Seurat object from single directory
    obj.list<-lapply(args$samples[1:n],create_seurat_obj_10X,input_dir=args$input_dir,args$data_dir,verbose=args$verbose)
    }
    
}else
  { # [if-else block 1]: If object is specified then run subset analysis
  
   if(args$verbose)
     {cat("Begin sub-clustering",'\n')}
   
   #Load Seurat object (containing all samples)
   obj<-readRDS(args$object)

   if(args$clusters=='none')
     { ## [if-else block 1d]: if clusters='none' then the object is entered to run rest of the pipeline skipping batch correction and the steps prior to it

          obj.integrted<-obj
          skip<-TRUE
      }else if(args$clusters=='all')
        {## [if-else block 1d]: if clusters='all' then no subsetting is done (useful for running analysis on manually merged obj)
		 sub<-obj
         
	 	}else
            {## [if-else block 1d]: Subset Seurat object to include only specific clusters defined by user
                sub<-subset(obj,idents=args$clusters)
            }

  if(skip!=TRUE)
   {
       
    sub.list<-list()
    sub.counts<-list()

   #Split the subsetted object by sample
   if(args$verbose)
    { cat("Splitting object by sample ID",'\n')}
   sub.list<-SplitObject(sub, split.by = "orig.ident")
   sample.id<-names(sub.list)

   #Extract counts data
   sub.counts<-lapply(X=1:length(sub.list), function(x) {return(sub.list[[x]]@assays$RNA@counts)})

   #Create Seurat object
   obj.list<-mapply(create_seurat_obj_from_counts_data,sub.counts,sample.id,verbose=args$verbose)

   rm(sub)
   rm(sub.list)
   rm(sub.counts)
   rm(sample.id)
  
  }
}
   
###(2) QC ###

##(2a) Generate raw count qc plots ## (optional)
if(skip!=TRUE)
{  if(args$qc_plots)
  {

#Get pre-filter histogram [sample-wise + aggregate]
if(!dir.exists(paste0(args$output_dir,'qc_plot/histograms/raw')))
        {dir.create(file.path("qc_plots","histograms","raw"), recursive = TRUE)}
out_path<-paste0(args$output_dir,'qc_plots/histograms/raw/')

x_lab<-c('LibrarySize','GeneCounts','MtPerc')
x_lim<-c(50000,2000,50)

if(args$thresholds=='')
{
    cutoff_low<-rep('none',3)
    cutoff_high<-rep('none',3)
}else
  {
   cutoff_low<-c(args$thresholds[4],args$thresholds[2],args$thresholds[1])
   cutoff_high<-c(args$thresholds[4],args$thresholds[3],args$thresholds[1])
  }

  
for(i in 1:length(x_lab))
{

 pre<-paste0(args$file_prefix,'Histogram_')
 # Get Histogram for each sample
 hist_list<-lapply(obj.list,get_histogram,out_path,prefix=pre,x_lab[i])

 #Get aggregate Histogram
 lab<-paste0(args$file_prefix,'AggregateHistogram')
 agg<-unlist(hist_list)
 print_histogram_abline(agg,out_path,label=lab,x_lab=x_lab[i],cutoff=c(cutoff_low[i],cutoff_high[i]))

 #Get aggregate Histogram low xlim
 lab<-paste0(args$file_prefix,'AggregateHistogram_low')
 print_histogram_abline(agg,out_path,label=lab,x_lab=x_lab[i],cutoff=c(cutoff_low[i],cutoff_high[i]),x=c(0,x_lim[i]))

}

rm(hist_list)
rm(agg)
rm(out_path)
rm(lab)
rm(x_lab)
rm(x_lim)

# Get Vlnplot and Scatterplots as well

obj_merged<-merge(obj.list[[1]],y=obj.list[2:length(obj.list)])

options(bitmapType='cairo')
png(paste0(args$output_dir,'qc_plot/',args$file_prefix,'genecount_VlnPlot.png'), width=16,height=8,units='in',res=300)
VlnPlot(obj_merged, features = "nFeature_RNA")
dev.off()

png(paste0(args$output_dir,'qc_plot/',args$file_prefix,'UMIcount_VlnPlot.png'), width=16,height=8,units='in',res=300)
VlnPlot(obj_merged, features = "nCount_RNA")
dev.off()

png(paste0(args$output_dir,'qc_plot/',args$file_prefix,'mt_perc_VlnPlot.png'), width=16,height=8,units='in',res=300)
VlnPlot(obj_merged, features = "perc.mt")
dev.off()

png(paste0(args$output_dir,'qc_plot/',args$file_prefix,'scatter_plots.png'),width=16,height=8,units='in',res=300)
plot1 <- FeatureScatter(obj_merged, feature1 = "nCount_RNA", feature2 = "perc.mt")
plot2 <- FeatureScatter(obj_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

rm(obj_merged)

}


##(2b) If user specified, thresholds, calculate post filtering metrics ##
if(length(args$thresholds)>1)
{
 # Get cell table
 tab_list<-lapply(obj.list[1:length(obj.list)],get_cell_table,mt.thres=args$thresholds[1], genecnt.thres=args$thresholds[2:3], libsize.thres=args$thresholds[4:5], verbose=args$verbose)
 tab<-data.frame()
 for(i in 1:length(tab_list))
  {
   tab<-rbind(tab,tab_list[[i]])
  }

 total<-colSums(tab)
 tab<-rbind(tab,"Total"=total)

 write.csv(tab,paste0(args$output_dir,args$file_prefix,"cell_qc_table.csv"))
 rm(tab_list)
 rm(tab)

 # Filter low quality cells
 
 obj_tab_list<-lapply(X=(1:length(obj.list)), FUN=function(x){filter_cells(obj.list[[x]],mt.thres=args$thresholds[1], genecnt.thres=args$thresholds[2:3], libsize.thres=args$thresholds[4:5], verbose=args$verbose)})
 
# Extract filtered data
 obj.list<-do.call(c,(lapply(obj_tab_list,`[[`,1)))
 cell.data<-do.call(rbind,(lapply(obj_tab_list, `[[`, 2)))
 cell.data2<-adorn_totals(cell.data,"col")

 write.csv(cell.data2,paste0(args$output_dir,args$file_prefix,"cell_threhold_table.csv"))

 rm(obj_tab_list)
 #cell.data<-do.call(c,lapply(cell_tab_list))

 # Get threhold plots
 cell.data$rank[order(cell.data$thr.cell.cnt,decreasing = TRUE)]<- c(1:nrow(cell.data))
 if(args$verbose){
    cat("Exporting cell count plots",sep='\n')}
 threshold_plots<-plot_threshold_effects(cell.data,args$thresholds,args$file_prefix)
 pdf(file=paste0(args$output_dir,args$file_prefix,"thresholds.pdf"),paper='a4',width=8)
 lapply(threshold_plots,print)
 dev.off()

 rm(threshold_plots)
 rm(cell.data)

 ## Get post filtering plots ## (optional)
 if(args$qc_plots)
 {
  # Get post threshold histograms
  if(!dir.exists(paste0(args$output_dir,'qc_plot/histograms/post_filter')))
        {dir.create(file.path("qc_plots","histograms","post_filter"), recursive = TRUE)}
  out_path<-paste0(args$output_dir,'qc_plots/histograms/post_filter/')

  x_lab<-c('LibrarySize','GeneCounts','MtPerc')
  x_lim<-c(50000,2000,50)
  
  cutoff_low<-c(args$thresholds[4],args$thresholds[2],args$thresholds[1])
  cutoff_high<-c(args$thresholds[4],args$thresholds[3],args$thresholds[1])

  for(i in 1:length(x_lab))
  {
   pre<-paste0(args$file_prefix,'post-filter_Histogram_')
   # Get Histogram for each sample
   hist_list<-lapply(obj.list,get_histogram,out_path,prefix=pre,x_lab[i])

   #Get aggregate Histogram
   lab<-paste0(args$file_prefix,'post-filter_AggregateHistogram')
   agg<-unlist(hist_list)
   print_histogram_abline(agg,out_path,label=lab,x_lab[i])

   #Get aggregate Histogram low xlim
   lab<-paste0(args$file_prefix,'post-filter_AggregateHistogram_low')
   print_histogram_abline(agg,out_path,label=lab,x_lab[i],x=c(0,x_lim[i]))

   }

  rm(hist_list)
  rm(pre)
  rm(lab)
  rm(agg)
  rm(x_lab)
  rm(x_lim)
  
  # Get Vlnplot and Scatterplots as well
  obj_merged<-merge(obj.list[[1]],y=obj.list[2:length(obj.list)])
  
  options(bitmapType='cairo')

  png(paste0(args$output_dir,'qc_plot/',args$file_prefix,'post_filter_genecount_VlnPlot.png'), width=16,height=8,units='in',res=300)
  VlnPlot(obj_merged, features = "nFeature_RNA")
  dev.off()

  png(paste0(args$output_dir,'qc_plot/',args$file_prefix,'post_filter_UMIcount_VlnPlot.png'), width=16,height=8,units='in',res=300)
  VlnPlot(obj_merged, features = "nCount_RNA")
  dev.off()

  png(paste0(args$output_dir,'qc_plot/',args$file_prefix,'post_filter_mt_perc_VlnPlot.png'), width=16,height=8,units='in',res=300)
  VlnPlot(obj_merged, features = "perc.mt")
  dev.off()

  png(paste0(args$output_dir,'qc_plot/',args$file_prefix,'post_filter_scatter_plots.png'),width=16,height=8,units='in',res=300)
  plot1 <- FeatureScatter(obj_merged, feature1 = "nCount_RNA", feature2 = "perc.mt")
  plot2 <- FeatureScatter(obj_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2
  dev.off()

  rm(obj_merged)

  }
}

# End script here if user specified only qc outputs

if(args$qc_only)
{
    quit(save='no')
}

###(3) Data pre-processing ###

##(a) Log normalize ##

if(args$verbose)
{cat('Performing Log normalization \n')}

n<-length(obj.list)
obj.list<-lapply(X=1:n,function(x){NormalizeData(obj.list[[x]],normalization.method = 'LogNormalize',scale.factor = 10^4, verbose =args$verbose)})

if(args$verbose)
{cat('Preparing batch correction using ',args$batch_correction, '\n')}

if(args$batch_correction!='none')
    {if(args$batch_correction=='cca-mnn')
        {
            assay<-'integrated'
            reduction<-'pca'
    
        }else if(args$batch_correction=='harmony')
            {
                assay<-'RNA'
                reduction<-'harmony'
                obj.merged<-merge(x=obj.list[[1]],y=obj.list[2:length(obj.list)])
                DefaultAssay(obj.merged)<-'RNA'
            }else
                  {
                    cat('You did not enter correct batch correction method, setting it to cca-mnn')
                    args$batch_correction<-'cca-mnn'
                    assay<-'integrated'
                    reduction<-'pca'
                }
  }else
    {
      obj.list<-obj    
    }

##(a) Find variable genes ##

if(args$verbose)
{cat('Finding Variable genes per sample \n')}
    
obj.list<-lapply(X=1:length(obj.list),function(x){FindVariableFeatures(obj.list[[x]],selection.method = 'vst', nfeatures=args$hvg, verbose = args$verbose)})

#Get variable genes table and plots for each sample
# Prepare output file
if(!dir.exists(paste0(args$output_dir,'variable_genes/')))
  {dir.create( paste0(args$output_dir,'variable_genes/'))}
out_path<-paste0(args$output_dir,'variable_genes/')

lapply(X=1:length(obj.list),function(x){get_var_genes(obj.list[[x]],out_dir=out_path,verbose=args$vebose)})

varplots<-lapply(obj.list[1:length(obj.list)],get_var_genes_plot,verbose=args$verbose)

 pdf(file=paste0(out_path,args$file_prefix,'VariableGenePlots_all_samples.pdf'),paper='a4')
 print(varplots)
 dev.off()

rm(out_path)
rm(varplots)

cat('Obj list size: ', length(obj.list), '\n')

###(4) Batch correction and integration ###

if(length(obj.list)>1)
 {
 # Integrate data with batch correction
 
 if(args$batch_correction=='cca-mnn')
 {
  if(args$verbose)
  {cat('Performing cca-mnn batch integration \n')}
  
  obj.integrated<-cca_batch_correction(obj.list,project.name=args$file_prefix, anchors=args$hvg, int.genes=args$batch_genes, verbose=args$verbose)

}else if(args$batch_correction=='harmony')
 {
     if(args$verbose)
     {cat('Performing harmony batch integration \n')}
  # Variable features for merged obj
  obj.merged<-FindVariableFeatures(obj.merged,selection.method = 'vst', nfeatures=args$hvg, verbose = args$verbose)
  
  if(!dir.exists(paste0(args$output_dir,'variable_genes/')))
    {dir.create( paste0(args$output_dir,'variable_genes/'))}

  out_path<-paste0(args$output_dir,'variable_genes/')

  #Var genes table
  get_var_genes(obj.merged,out_dir=out_path,verbose=args$vebose)
  
  # Var genes plot
  varplots<-get_var_genes_plot(obj.merged,verbose=args$verbose)

   pdf(file=paste0(out_path,args$file_prefix,'VariableGenePlots.pdf'),paper='a4')
   print(varplots)
   dev.off()

  # Scaling
  all.features<-rownames(obj.merged)
  obj.merged<-ScaleData(obj.merged,features=all.features, verbose=args$verbose,assay = 'RNA')
  
   # PCA
   obj.merged<-RunPCA(obj.merged, npcs=50, verbose=args$verbose,assay = 'RNA')
   
   #Run Harmony
   obj.integrated<-RunHarmony(obj.merged,
                  group.by.vars = "orig.ident",
                  plot_convergence = FALSE,
                  assay.use = 'RNA')
 }
}else
  {
    # This allows skipping batch correction in case user wants to run the rest of the pipeline on an already batch corrected object (most likely batch corrected by another method)
    
   obj.integrated<-obj.list
  }


###(5) Add meta data ### (optional)

if(args$verbose)
{cat('Adding metdata \n')}

# Add sample column
obj.integrated[['sample']]<-obj.integrated[['orig.ident']]

# Add metadata from file (if specified)
if(args$meta_file!='')
{
obj.integrated<-add_metadata(obj.integrated,args$meta_file)
}

rm(obj.list)

### Save batch corrected obj ### (optional)
if((args$save=='integrated')||(args$save=='both'))
 {
  #Saving object pre-clustering for clustering optimization (if needed)   
  saveRDS(obj.integrated,file=paste0(args$output_dir,args$file_prefix,"integrated_only.rds"))
 }

 } # skip until this point if only running post processing
 
###(6) Pre-clustering processing ###

if(clustree)
{
 #Save a copy of obj before running PCA to generate clustree geneplots on RNA assay (optional)
 obj.integrated_RNA<-obj.integrated
}

##(6a) Run PCA ## (skip if ran harmony or if reduction already exists)

if(is.null(obj.integrated@reductions$pca)==TRUE)
{
 if(is.null(obj.integrated@reductions$harmony)==TRUE)
    {
       # Perform PCA since no reduction assay found
       DefaultAssay(obj.integrated)<-'integrated'

        if(args$verbose)
        {cat('Running PCA on integrated data')}
          all.features<-rownames(obj.integrated)
          obj.integrated<-ScaleData(obj.integrated,features=all.features, verbose=args$verbose)
         
        # Run PCA
        obj.integrated<-RunPCA(obj.integrated, npcs=50,ndims.print = 1:15, verbose=args$verbose)

    }
}


##(6b) PCA Plots ##

# Elbow Plot

options(bitmapType='cairo')
png(file=paste0(args$output_dir,args$file_prefix,"PCA_elbow_plots.png"),width = 11,height = 8.5, units='in', res=300)
  ElbowPlot(obj.integrated,ndims=50,reduction=reduction)+ggtitle(paste0(args$file_prefix,'ElbowPlot'))
  dev.off()

# PCA gene plot (print 4 PCs per page)
pc_genes_plot_list<-lapply(args$pca_dimensions,function(x){ VizDimLoadings(obj.integrated,dims=x,ncol=1,reduction=reduction) })

pca_plots<-marrangeGrob(pc_genes_plot_list, nrow=2, ncol=2)
ggsave(paste0(args$output_dir,args$file_prefix,"PC_gene_plots.pdf"), width=8.5, height=11, units = "in", pca_plots)


###(7) Clustering Optimization ### (optional)

##(7a) Generating clustree plots ## (optional)
if(clustree)
{

if(args$verbose)
{cat('Running clustree \n')}

res<-seq(0.1,1.2,by=0.1)
pc<-c(15,20,25,30,35,40,50)

for(i in 1:length(pc))
{
#dims<-seq(1,pc[i],by=1)
obj_clustree<-NULL
clus_run=paste0(args$file_prefix,'PC',pc[i])
obj_clustree<-iterative_clus_by_res(obj.integrated, res=res,dims_use=1:pc[i],reduction=reduction,assay=assay,verbose=args$verbose)

col<-grep('snn',colnames(obj_clustree@meta.data))
prefix<-gsub('[0-9].+','',colnames(obj_clustree@meta.data))[col[1]]

print_clustree_png(obj_clustree, prefix=prefix, out_dir=args$output_dir, file_prefix=clus_run, verbose=args$verbose)

# Generate clustree geneplots on integrated assay
print_geneplots_on_clustree(obj_clustree,genes=args$gene_list,prefix=prefix,assay=assay , fun_use='median',out_dir= args$output_dir, file_prefix=clus_run, verbose=FALSE)
}

# Generate clustree geneplots on RNA assay using a copy of batch corrected object

# Skip this step if harmony integration was run

if(args$batch_correction!='harmony')
{
 DefaultAssay(obj.integrated_RNA)<-'RNA'

 cat('Scaling \n')

 all.genes<-rownames(obj.integrated_RNA)
 obj.integrated_RNA<-ScaleData(obj.integrated_RNA,features=all.genes)

cat('Running PCA on RNA assay \n')

# Get var features from RNA assay if integrated assay is null or has 0 var.features

if((is.null(obj.integrated@assays$integrated@var.features))==FALSE)
{
 v<-obj.integrated@assays$integrated@var.features %>% length()
 
 if(v==0) #Integrated assay exists but has 0 var.features
     {
         var_genes<-obj.integrated_RNA@assays$RNA@var.features
     }else #Integrated assay exists and has >0 var.features
        {
          var_genes<-obj.integrated_RNA@assays$integrated@var.features
        }
}else # Integrated assay doesn't exist (in case of harmony integration)
  {
      var_genes<-obj.integrated_RNA@assays$RNA@var.features
  }

obj.integrated_RNA<-RunPCA(obj.integrated_RNA,assay='RNA',features=var_genes)

cat('Generating clustree geneplots on RNA assay \n')
for(i in 1:length(pc))
{
obj_clustree<-NULL
clus_run=paste0(args$file_prefix,'PC',pc[i])
obj_clustree<-iterative_clus_by_res(obj.integrated_RNA, res=res,dims_use=1:pc[i],verbose=args$verbose,assay='RNA')

print_geneplots_on_clustree(obj_clustree,genes=args$gene_list, fun_use='median',assay='RNA',prefix='RNA_snn_res.', out_dir=args$output_dir,file_prefix=clus_run, verbose=FALSE)
}

rm(clus_run)
rm(obj_clustree)
rm(obj.integrated_RNA)

}

}

##(7b) Generate silhouette plots ## (optional)
if(sil)
{

if(args$verbose)
 {cat('Generating Silhouette plots \n')}


res<-seq(0.1,1.2,by=0.1)
pc<-c(15,20,25,30,35,40,50)

DefaultAssay(obj.integrated)<-assay

for(i in 1:length(pc))
{
for(j in 1:length(res))
   {
   sil_run<-paste0(args$file_prefix,'PC',pc[i],'_res',res[j])
   obj_sil<-iterative_clus_by_res(obj.integrated,dims_use=1:pc[i],res=res[j],verbose=args$verbose,reduction=reduction,assay=assay)
    get_silhouette_plot(s.obj=obj_sil,reduction=reduction,dims=1:pc[i],out_dir=args$output_dir,file_prefix=sil_run,verbose=args$verbose)
   }
}

rm(pc)
rm(res)
rm(sil_run)
rm(obj_sil)

}


###(8) CLUSTERING ###

# In the pipeline the clusters are generated by default using 50 PCs at resolution 0.5. These parameters can be modified by the user

DefaultAssay(obj.integrated)<-assay

if(args$verbose)
{cat('Clustering \n')}

#Find Clusters
  obj.integrated<-FindNeighbors(obj.integrated, dims=args$pca_dimensions,assay=assay,reduction=reduction)
  obj.integrated<-FindClusters(obj.integrated, res=args$cluster_resolution)

#Generate UMAP
  if(args$verbose)
  {cat("Genetrating UMAP plots", sep='\n')}
  obj.integrated<-RunUMAP(obj.integrated, reduction=reduction,assay=assay, dims=args$pca_dimensions)

#Print UMAPs
if(!dir.exists(paste0(args$output_dir,'umaps/')))
        {dir.create( paste0(args$output_dir,'umaps/'))}
    out_dir<-paste0(args$output_dir,'umaps/')

customized_umap(obj.integrated, umap_cols=NULL,label=TRUE,title=NULL, group=NULL,split=NULL,dot=0.3,save=TRUE,out=out_dir,file_prefix=args$file_prefix)

customized_umap(obj.integrated, umap_cols=NULL,label=TRUE,title=NULL, group='orig.ident',split=NULL,dot=0.3,save=TRUE,out=out_dir,file_prefix=args$file_prefix)

if(args$meta_file!='')
{
 customized_umap(obj.integrated, umap_cols=NULL,label=TRUE,title=NULL, group='sex',split=NULL,dot=0.3,save=TRUE,out=out_dir,file_prefix=args$file_prefix)
 
 customized_umap(obj.integrated, umap_cols=NULL,label=TRUE,title=NULL, group='age',split=NULL,dot=0.3,save=TRUE,out=out_dir,file_prefix=args$file_prefix)
 
 customized_umap(obj.integrated, umap_cols=NULL,label=TRUE,title=NULL, group='normal_tumor',split=NULL,dot=0.3,save=TRUE,out=out_dir,file_prefix=args$file_prefix)
}
 
## Save clustered seurat objects ##

if((args$save=='clustered')||(args$save=='both'))
{
 
 if(args$verbose)
  {cat("Exporting Seurat integrated object",sep='\n')}
 
 saveRDS(obj.integrated,file=paste0(args$output_dir,args$file_prefix,"clustered.rds"))
}

###(9) Post-clustering downstream analysis ###

##(9a) Cell proportion table ##
tab<-table(Idents(obj.integrated),obj.integrated@meta.data$orig.ident)

#Convert to data frame
tab<-as.data.frame.matrix(tab)
tab<-cbind(rownames(tab),tab)
colnames(tab)[1]<-'cluster'

#Add totals
tab2<-adorn_totals(tab, c("row","col"))
write.csv (tab2,paste0(args$output_dir,args$file_prefix,"cells_by_cluster_by_sample.csv"),row.names = FALSE)


##(9b) Differential gene expression ##

clusters_num<-levels(obj.integrated)
DefaultAssay(obj.integrated)<-"RNA"
marker.list<-differential_gene_exp(obj.integrated, clusters=clusters_num, out_dir=args$output_dir,file_prefix=args$file_prefix,verbose=args$verbose)

##(9c) Feature Plots ##

if(!dir.exists(paste0(args$output_dir,'feature_plots/')))
        {dir.create( paste0(args$output_dir,'feature_plots/'))}
    out<-paste0(args$output_dir,'feature_plots/')

# Extract feature list from marker list file

features.list<-list()

for(i in 1:length(marker.list))
{
 features.list[[i]]<-as.vector(marker.list[[i]][,'gene'])
}

# Generating pdf Plots
if(args$verbose)
        {cat("Creating Feature plots (pdf format)", sep='\n')}

plot.list<-list()

for(i in 1:length(marker.list))
{
plot.list<-get_features_plot(features.list[[i]],obj.integrated, top=20)
p1<-marrangeGrob(plot.list, nrow=3, ncol=3)
ggsave(paste0(out,args$file_prefix,"FeaturePlot_cluster",clusters_num[i],".pdf"),width=11, height=8.5, units = "in", p1)
dev.off()
}

# Saving plots in png format as well
if(args$verbose)
        {cat("Creating Feature plots (png format)", sep='\n')}

 DefaultAssay(obj.integrated)<-'RNA'
 n<-9 # Num of features to print per cluster

for(i in 1:length(features.list))
{

#png(paste0(out,file_name,"FeaturePlot_cluster",clusters_num[i],".png"),width=11, height=8.5, units = "in",res=600)
features<-features.list[[i]]
p2<-FeaturePlot(obj.integrated,features[1:n],min.cutoff='q30')
ggsave(paste0(out,args$file_prefix,"FeaturePlot_cluster",clusters_num[i],".png"),width=11, height=8.5, units = "in",p2)
}

rm(clusters_num)
rm(marker.list)
rm(plot.list)
rm(p1)
rm(p2)
rm(out)

##(9d) Heatmaps ##

# Generating heatmap for top 15 markers per cluster

DefaultAssay(obj.integrated)<-assay
heatmap_features<-vector()
top_h<-15
for(i in 1:length(features.list))
{

heatmap_features<-append(heatmap_features,features.list[[i]][1:top_h])

}

rm(features.list)

#png(paste0(file_name,"Heatmap.png"),width=55,height=34,units="in",res=300)
h1<-DoHeatmap(obj.integrated,features=heatmap_features)
ggsave(paste0(args$file_prefix,"Heatmap_integrated.pdf"),path=args$output_dir,width=22,height=17,units="in",h1)

rm(h1)
 
# Plot Heatmap for RNA assay

DefaultAssay(obj.integrated)<-'RNA'
obj.integrated<-ScaleData(obj.integrated, features=rownames(obj.integrated))
h2<-DoHeatmap(obj.integrated,features=heatmap_features)
ggsave(paste0(args$file_prefix,"Heatmap_RNA.pdf"),path=args$output_dir,width=22,height=17,units="in",h2)

rm(h2)

##(9e) Dendrogram ###

plot_dendrogram(obj.integrated,args$file_prefix,features=heatmap_features,out_dir=args$output_dir)

rm(heatmap_features)
rm(obj.integrated)

# Print session info
sessionInfo()
cat('\n')

cat('End of script ','\n')

