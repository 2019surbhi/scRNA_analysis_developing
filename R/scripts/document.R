#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

x <- readRDS(args[1])
if(length(x$thresholds) > 1){
	th <- x$thresholds
	names(th) <- c("perc.mt", "low.gene.cutoff", "high.gene.cutoff", "low.lib.size", "high.lib.size")
	x <- c(x,th)
	x$thresholds <- NULL
}
x$gene_list <- paste(x$gene_list, collapse=",")
if(length(x$pca_dimensions) > 1) x$pca_dimensions <- paste0(min(x$pca_dimensions),":", max(x$pca_dimensions))
x[1:2] <- NULL
write.csv(x)
