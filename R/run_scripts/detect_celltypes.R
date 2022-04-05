## This script finds cell beloning to celltypes using AUC as input to the evaluation framework

##########
### Load parameters and packages
##########
message(" Load parameters and packages ")

library(magrittr)
library(scUtils)
library(Seurat)
library(AUCell)
library(tidyverse)

# get params-filename from commandline
command_args<-commandArgs(TRUE)
param_file = command_args[1]
# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# read features to exclude
features_exclude_list= jsonlite::read_json(parameter_list$genes_to_exclude_file)
features_exclude_list = lapply(features_exclude_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# read signatures
signaturelist = jsonlite::read_json(path = parameter_list$genes_to_exclude_file)
signaturelist = lapply(signaturelist,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# load seurat
seurat_merged = readRDS(paste0(parameter_list$data_path,parameter_list$merged_file))

##########
### Run AUCEll based function
##########

auc_mat_celltypes = map_celltype_signatures2(exprMatrix=seurat_object_raw@assays$RNA@counts,block_size=10000,aucMaxRank_n=parameter_list$auc_max_rank,
                                             gene_set_list=gene_set_list,min_rowSum=20,global_seed =parameter_list$global_seed)
# save
data.table::fwrite(data.table::as.data.table(auc_mat_celltypes),file = paste0(parameter_list$mapped_celltypes_auc_file),sep="\t")

##########
### Run AUCEll based function
##########
