Here are examples of how to run the script

**1) Using raw data from cellranger output**

  **(a) Inputs from Single flowcell**

/path/to/script/scRNA_pipeline_run_final_HPC.R -i /path/to/cellranger/ouput/flowcell1/ -r run_tag -o path/to/output/dir/ -s sample1,sample2,sample3,sample4 -t mt,gene_low,gene_high,lib_low,lib_high -c 20 -v -q -u clustree -g all -d 1:30 -e 0.2

*Notes:*
* Leave out -s argument to load all samples
* Leave out -q and -u to not generate plot qc and clustree (saves time)
* -c is number of cores and change it according to dataset
* -m can be used to specify memory for globals (I use Seurat's suggested parallel processing in this pipeline so need to specify size global object passed to functions; usually set it to 50 or 100GB but may be higher depending on the dataset)
* -g is used to decide whether sample anchors, or all genes are integrated post batch correction. Integrating all genes is problematic for larger datasets.
* By default clustering is done using 1st 50 PCs at resolution 0.5; a different resolution can be specified used -e and different number of PCs using -d like illustrated in the example.

 **(b) Inputs from Multiple flowcells**

/path/to/script/scRNA_pipeline_run_final_HPC.R -i /path/to/cellranger/ouput/flowcell1/:path/to/cellranger/output/flowcells2/ -r run_tag -o path/to/output/dir/ -s sample1,sample2,sample3:sample11,sample12,sample13 -t mt,gene_low,gene_high,lib_low,lib_high -c 20 -v -q -u clustree -g all

*Notes:*
* again leave out -s argument if all samples are to be loaded from both (or more) flowcells

**2) Using Seurat object (Subset analysis)**

/path/to/script/scRNA_pipeline_run_final_HPC.R -b /path/to/seurat/object/seruat_object.rds -r run_tag -o path/to/output/dir/ -l 0,1,2,5,6 -c 20 -v -u clustree -g all

*Notes:*
* -l specifies clusters to be inlcuded in the subset analysis
* Use -t if additional filtering needs to be applied (cannot be lower than the threshold values used for parent object)
