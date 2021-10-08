#!/usr/bin/env Rscript

# This is a script that facilitates running the Ting Lab scRNA seq pipeline on HPC for processing data generated on 10XGenomics platform 

source('./scRNA_pipeline_functions_final_HPC.R')

library(argparser)

### User defined arguments ###
parser<-arg_parser(name="scRNA_pipepine_qc.R",description="Processing multiple scRNA seq sample using Seurat pipeline")


parser<-add_argument(
  parser,
  arg='--input_dir',
  short = '-i',
  type="character",
  default="./",
  help="Enter the path(s) for the dataset(s). If multiple separate by : .Must have trailing '/'. Default: [./]")

parser<-add_argument(
  parser,
  arg='--out_dir',
  short = '-o',
  type="character",
  default="./",
  help="location+name for output directory. Must contain trailing '/'. Default: [./data/cluster_plots/]")

parser<-add_argument(
  parser,
  arg='--samples',
  short = '-s',
  type="character",
  default="all",
  help="Enter the sample names to be included from each datasets separated by , while datasets separated by : Default is 'all' samples")

parser<-add_argument(
 parser,
 arg='--object',
 short='-b',
 type='character',
 default='',
 help='Enter Seurat object with full path - from multiple dataset or subcluster analysis')

parser<-add_argument(
  parser,
  arg='--matrix_name',
  short = '-n',
  type="character",
  default="filtered_feature_bc_matrix",
  help="Name of the matrix file/directory with no file extension. Eg: 'raw_feature_bc_matrix'. Default: [filtered_feature_bc_matrix]")

parser<-add_argument(
  parser,
  arg='--run_name',
  short = '-r',
  type="character",
  default="",
  help="Title used to identify all plots that correspond to the current run instance. Default: None")

parser<-add_argument(
  parser,
  arg='--thresholds',
  short = '-t',
  type="character",
  default='',
  help="Set of thresholds to be applied to all samples for QC. Format:'<% MT genes>,<gene count minimum>,<gene count maximum>,<libsize low cutoff>,<libsize high cutoff>'")

parser<-add_argument(
  parser,
  arg='--data_dir',
  short = '-z',
  flag=TRUE,
  help="Flag set to read cellranger output from directory instead of hdf5 file")


parser<-add_argument(
  parser,
  arg='--cores',
  short = '-c',
  type='numeric',
  default=4,
  help= "Number of cores to use when processing samples in parallel.")

parser<-add_argument(
  parser,
  arg='--mem',
  short = '-m',
  type='numeric',
  default=40,
  help= "Memory allocated for gloabls in GB Default is 10")

parser<-add_argument(
  parser,
  arg='--verbose',
  short = '-v',
  flag=TRUE,
  help="Set flag for printing output messages during the run")

parser<-add_argument(
  parser,
  arg='--clusters',
  short='-l',
  type="character",
  default='',
  help="Enter which clusters (cluster numbers separated by ,) for subcluster analysis")

argv <- parse_args(parser)

#Print run parameters
cat("List of arguments:",'\n')
argv
cat('\n')

#Save run parameters
write.table(argv,"run_parameters.txt",sep='\t',col.names=FALSE)

### Argument validation and processing ### 

#Process thresholds
if(argv$thresholds!='')
{
 argv$thresholds= as.numeric(unlist(strsplit(argv$thresholds,split = ',',fixed = TRUE)))
}

# Prepare output directory
if(! dir.exists(argv$out_dir) ){
  dir.create(argv$out_dir)
}

#Validate Sublcuster numbers
if(argv$clusters!='')
{
argv$clusters=as.numeric(unlist(strsplit(argv$clusters,split = ',',fixed = TRUE)))}

# Add separater to run name
if(argv$run_name!=''){
  argv$run_name<- gsub(' ','_',argv$run_name)
 }
argv$run_name<-paste0(argv$run_name,'_')

# Print processed arguments
cat("List of processed arguments:",'\n')
argv
cat('\n')

#Save arguments
saveRDS(argv,paste0(argv$out_dir,argv$run_name,"arguments.rds"))

### Main Script ###

#Facilitate parallel processing

plan("multicore", workers=argv$cores)
options(future.globals.maxSize=((argv$mem*1000)*1024^2))
#n<-length(argv$samples)
obj.list<-list()

#Create timer
timer<-createTimer(verbose=FALSE)

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
  {cat("you did not enter any samples, Running pipeline on all samples from the specified input directories",sep='\n')}
   s<-lapply(paths.list,list.files)
   }else{
      #Create a list of samples corresponding to each path
      s<-strsplit(argv$samples,split=':') %>% unlist %>% lapply(FUN=function(x) {strsplit(x,split=',') %>%  unlist})
         }
# Create a list of seurat object
for(i in 1:length(paths.list))
 {
  n<-length(s[[i]])

   s.list[[i]]<-lapply(s[[i]][1:n], create_seurat_obj_10X,input_dir=paths.list[i],argv$matrix_name,argv$data_dir,verbose=argv$verbose)

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
       {cat("you did not enter any samples, Running pipeline on all samples in the input directory",sep='\n')}
      argv$samples<-list.files(argv$input_dir)
     }else{
           argv$samples<-unlist(strsplit(argv$samples,split = ',',fixed = TRUE))
          }
  n<-length(argv$samples)
  #Create Seurat object from single directory
   obj.list<-lapply(argv$samples[1:n], create_seurat_obj_10X,input_dir=argv$input_dir,argv$matrix_name,argv$data_dir,verbose=argv$verbose)

 }
}else{
  if(argv$verbose)
{cat("Begin sub-clustering",'\n')}

#Load Seurat object (containing all samples)
obj<-readRDS(argv$object)

sub.list<-list()
sub.counts<-list()

#Subset Seurat object to include only specific clutsers defined by user
sub<-subset(obj,idents=argv$clusters)

#Split the subsetted object by sample
if(argv$verbose)
{ cat("Splitting object by sample ID",'\n')}

sub.list<-SplitObject(sub, split.by = "orig.ident")
sample.id<-unique(sub@meta.data$orig.ident)

#Extract counts data
for(i in 1:length(sub.list))
{sub.counts[[i]]<-sub.list[[i]]@assays$RNA@counts}

#Create Seurat object
obj.list<-mapply(create_seurat_obj_from_counts_data,sub.counts,sample.id,verbose=argv$verbose)

rm(sub)
rm(sub.list)
rm(sub.counts)
rm(sample.id)

}
 
#Get pre-filter histogram [sample-wise + aggregate]

timer$start("qc plots")
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

 pre<-paste0(argv$run_name,'Histogram_')
 # Get Histogram for each sample
 hist_list<-lapply(obj.list,get_histogram,out_path,prefix=pre,x_lab[i])

 #Get aggregate Histogram
 lab<-paste0(argv$run_name,'AggregateHistogram')
 agg<-unlist(hist_list)
 print_histogram_abline(agg,out_path,label=lab,x_lab[i],cutoff=c(cutoff_low[i],cutoff_high[i]))

 #Get aggregate Histogram low xlim
 lab<-paste0(argv$run_name,'AggregateHistogram_low')
 print_histogram_abline(agg,out_path,label=lab,x_lab[i],cutoff=c(cutoff_low[i],cutoff_high[i]),x=c(0,x_lim[i]))

}

rm(hist_list)
rm(agg)
rm(out_path)
rm(lab)
rm(x_lab)
rm(x_lim)

#hist_list<-lapply(X=1:length(obj.list),FUN=function(x){get_histogram(obj.list[[x]],out_path,argv$run_name)})
#get_agg_histograms(hist_list,argv$run_name,out_path)

# Generate  violin plot
n<-length(obj.list)
obj<-merge(obj.list[[1]],obj.list[2:n])

VlnPlot(obj,features = "nFeature_RNA",pt.size=0)
ggsave(paste0(argv$out_dir,"qc_plots/",argv$run_name,"Violin plot_gencnt.png"),width=11,height=8.5)

VlnPlot(obj,features = "nCount_RNA",pt.size=0)
ggsave(paste0(argv$out_dir,"qc_plots/",argv$run_name,"Violin plot_libsize.png"),width=11,height=8.5)

VlnPlot(obj,features = "perc.mt",pt.size=0)
ggsave(paste0(argv$out_dir,"qc_plots/",argv$run_name,"Violin plot_mt.png"),width=11,height=8.5)

#Generate scatter plot

#plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "perc.mt")
#plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#ggsave(paste0(argv$out_dir,"qc_plots/",argv$run_name,"scatter1.png"),plot1)
#ggsave(paste0(argv$out_dir,"qc_plots/",argv$run_name,"scatter2.png"),plot2)

rm(obj)

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
write.csv(tab,paste0(argv$out_dir,"cell_qc_table.csv"))
rm(tab_list)
rm(tab)

# Filter low quality cells
obj_tab_list<-lapply(X=(1:length(obj.list)), FUN=function(x){filter_cells(obj.list[[x]],mt.thres=argv$thresholds[1], genecnt.thres=argv$thresholds[2:3], libsize.thres=argv$thresholds[4:5], verbose=argv$verbose)})

obj.list<-do.call(c,(lapply(obj_tab_list,`[[`,1)))
cell.data<-do.call(rbind,(lapply(obj_tab_list, `[[`, 2)))

write.csv(cell.data,paste0(argv$out_dir,"cell_threhold_table.csv"))

rm(obj_tab_list)

#cell.data<-do.call(c,lapply(cell_tab_list))

# Get threhold plots
cell.data$rank[order(cell.data$thr.cell.cnt,decreasing = TRUE)]<- c(1:nrow(cell.data))
if(argv$verbose){
    cat("Exporting cell count plots",sep='\n')}
threshold_plots<-plot_threshold_effects(cell.data,argv$thresholds,argv$run_name)
pdf(file=paste0(argv$out_dir,argv$run_name,"thresholds.pdf"),paper='a4',width=8)
lapply(threshold_plots,print)
dev.off()

rm(threshold_plots)
rm(cell.data)


# Get post threshold histograms
if(!dir.exists(paste0(argv$out_dir,'qc_plot/histograms/post_filter')))
        {dir.create(file.path("qc_plots","histograms","post_filter"), recursive = TRUE)}
    out_path<-paste0(argv$out_dir,'qc_plots/histograms/post_filter/')

x_lab<-c('LibrarySize','GeneCounts','MtPerc')
x_lim<-c(50000,2000,50)

for(i in 1:length(x_lab))
{
 pre<-paste0(argv$run_name,'post-filter_Histogram_')
 # Get Histogram for each sample
 hist_list<-lapply(obj.list,get_histogram,out_path,prefix=pre,x_lab[i])

 #Get aggregate Histogram
 lab<-paste0(argv$run_name,'post-filter_AggregateHistogram')
 agg<-unlist(hist_list)
 print_histogram(agg,out_path,label=lab,x_lab[i])

 #Get aggregate Histogram low xlim
 lab<-paste0(argv$run_name,'post-filter_AggregateHistogram_low')
 print_histogram(agg,out_path,label=lab,x_lab[i],x=c(0,x_lim[i]))

}

rm(hist_list)
rm(pre)
rm(lab)
rm(agg)
rm(x_lab)
rm(x_lim)

#lapply(X=1:length(obj.list),FUN=function(x){get_histogram(obj.list[[x]],out_path,paste0(argv$run_name,'post_filter_'))})
}

timer$stop("qc plots", comment="generate qc plots and tables (pre+post filtering")

# Extract and save timer components
time_dat<-getTimer(timer)
time_dat
write.csv(time_dat,paste0(argv$out_dir,argv$run_name,"timer.csv"))

cat('End of script ','\n')
