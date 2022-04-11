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
# source functions
source("R/evaluation_functions.R")

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
signaturelist = jsonlite::read_json(path = parameter_list$celltype_signature_file)
signaturelist = lapply(signaturelist,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# load seurat
seurat_merged = readRDS(paste0(parameter_list$merged_file))

##########
### Run AUCEll based function
##########

auc_mat_celltypes = map_celltype_signatures2(exprMatrix=seurat_merged@assays$RNA@counts,block_size=10000,aucMaxRank_n=parameter_list$auc_max_rank,
                                             gene_set_list=signaturelist,min_rowSum=10,global_seed =parameter_list$global_seed)
# save
data.table::fwrite(data.table::as.data.table(auc_mat_celltypes),file = paste0(parameter_list$auc_backup_file),sep="\t")

##########
### Use matrix and subset to cells based on threshold
##########
if(parameter_list$alpha<1){parameter_list$alpha=1} # factor to manually increase thresholds because the default AUCell Approach is not stringent enough in most cases

## calculate thresholds for each mapped celltype
threshold_fun = function(celltype_values,thrP,smallestPopPercent,auc_min_pos_thresh=0.05,auc_max_pos_thresh=0.15,alpha=1){
  glProb <- 1-(thrP/length(celltype_values) + smallestPopPercent)
  pos_thresh1 = as.numeric(qnorm(glProb,mean=mean(celltype_values),sd=sd(celltype_values)) * alpha)
  pos_thresh1= max(auc_min_pos_thresh,min(pos_thresh1,auc_max_pos_thresh))
  return(pos_thresh1)
}
thresholds = apply(auc_mat_celltypes,2,threshold_fun,thrP = parameter_list$thrP, smallestPopPercent =parameter_list$smallestPopPercent,
                   alpha=parameter_list$alpha,auc_max_pos_thresh = parameter_list$auc_max_pos_thresh,auc_min_pos_thresh= parameter_list$auc_min_pos_thresh)

## multiple auc matrix with gene set lengths (short signatures often have high scores which can distort comparison across signatures)
auc_mat_celltypes_mult = t(t(auc_mat_celltypes)*sapply(signaturelist,length))
# for each cells: which celltype scored highest (including the length factor):
max_val = apply(auc_mat_celltypes_mult,1,function(x){max(x)[1]})
# build mask with TRUE where max
max_mask=auc_mat_celltypes_mult==max_val
# build second mask with TRUE where above threshold
thresh_mask = t(t(auc_mat_celltypes)>thresholds)
# add masks and take 'intersection'
merged_masks = max_mask+thresh_mask
merged_masks[merged_masks==1]=0
merged_masks[merged_masks==2]=1
# there are a few cells with two max values --> we just take the first
#length(which(rowSums(merged_masks)>1))
# get the max celltype for each cell (that is also above threshold)
max_idx = data.frame(Cell_ID = rownames(merged_masks),Signature_CellType = as.character(apply(merged_masks,1,function(x){names(x)[x==1][1]}))) %>% na.omit()
# put into a list per cluster
all_mapped_cells = base::split(max_idx$Cell_ID,f=max_idx$Signature_CellType) # or: avg_log2fc or fc_mast

##########
### Save
##########

# save as separate json
writeList_to_JSON(list_with_rows = all_mapped_cells,filename = paste0(parameter_list$integration_folder_path,parameter_list$detected_cells_filename))

## add to seurat:
tmp_meta = dplyr::left_join(seurat_merged@meta.data,max_idx,by=c("Cell_ID"="Cell_ID"))
rownames(tmp_meta) = tmp_meta$Cell_ID


message(Sys.time(),": Finalized celltype detection")
