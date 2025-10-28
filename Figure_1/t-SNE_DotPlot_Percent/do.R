
### Deal arguements
args <- commandArgs(T)

file    <- 'parameter.yaml'
outdir  <- './'

handlers <- list("bool#no" = function(x){if ( x %in% c("false", "FALSE") ) FALSE else x}, "bool#yes" = function(x){if ( x %in% c("true", "TRUE") ) TRUE else x})
parameter <- yaml::yaml.load_file( file, handlers = handlers)

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)

library(future)
options(future.globals.maxSize = 100 * 1024 * 1024^2)
plan("multiprocess", workers = 4)
if ( ! is.na(outdir) ) setwd(outdir)

obj <- Load('obj.Rda')

## t-SNE
source('plotting.R')
PlotCluster(obj, reduction = 'umap', outpref = "UMAP" )
PlotCluster(obj, reduction = 'tsne', outpref = "tSNE" )

## Percent
source('cluster.R')
StatCluster(obj)

## DotPlot
features <- readLines('gene.list')
PlotDotPlot(obj, features = features, group.by = "seurat_clusters", outfile = "DotPlot.pdf", is.use.name = TRUE)


