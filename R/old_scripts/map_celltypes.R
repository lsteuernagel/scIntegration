# run processing script
# this script runs some runtime-intensive stpes sequentially (SCtransform on different feature set sizes). For very large datasets this could be reworked.

message(Sys.time(),": Starting mapping of celltypes." )

##########
### Load parameters and packages
##########

# get params-filename from commandline
command_args<-commandArgs(TRUE)
paramfile_name = command_args[1]
# read all parameters and filepaths
message("Reading parameters from: ",paramfile_name)
parameter_list = jsonlite::read_json(paramfile_name)


message("Reading parameters ...")
if(length(parameter_list)==0){
  stop("Please provide a valid parameter list!")  
}
#message(parameter_list)

# source functions and load libs
require(tidyverse)
require(Seurat)
require(Matrix)
require(SeuratDisk)
require(AUCell)
source("evaluation/evaluation_functions.R")
source("utils.R")
##########
### Map celltypes and save matrix with values
##########
message("Load data: ")
# read seurat object
#seurat_object_raw <-SeuratDisk::LoadH5Seurat(paste0(parameter_list$scHarmonize_path,parameter_list$key,"/",parameter_list$seurat_data,"/",parameter_list$key,"_raw.h5Seurat"))
seurat_object_raw <-SeuratDisk::LoadH5Seurat(paste0(parameter_list$seurat_object_raw_path))

# load gene sets
gene_set_list=readRDS(parameter_list$gene_set_list_file)
# if(is.null(parameter_list$gene_set_list)){
#   gene_set_list=readRDS(parameter_list$gene_set_list_file)
# }else{
#   gene_set_list=parameter_list$gene_set_list
# }
message("Run mapping: ")
# get mappings
auc_mat_celltypes = map_celltype_signatures2(exprMatrix=seurat_object_raw@assays$RNA@counts,block_size=10000,aucMaxRank_n=parameter_list$auc_max_rank,
                                             gene_set_list=gene_set_list,min_rowSum=20,global_seed =parameter_list$global_seed)
# save
data.table::fwrite(data.table::as.data.table(auc_mat_celltypes),file = paste0(parameter_list$mapped_celltypes_auc_file),sep="\t")

##########
### Use matrix and subset to cells based on threshold
##########
# 
# metadata = seurat_object_raw@meta.data
# 
# all_celltypes=colnames(auc_mat_celltypes)
# all_mapped_cells = list()
# for(i in 1:length(all_celltypes)){
#   ## subset to cells above threhsold T ? ---> Probably makes sense
#   celltype_values = auc_mat_celltypes[,i]
#   #print(hist(celltype_values))
#   current_celltype = all_celltypes[i]
#   
#   # subset
#   # get threshold 
#   if(parameter_list$alpha<1){parameter_list$alpha=1}
#   glProb <- 1-(parameter_list$thrP/length(celltype_values) + parameter_list$smallestPopPercent)
#   pos_thresh = qnorm(glProb,mean=mean(celltype_values),sd=sd(celltype_values)) * parameter_list$alpha
#   if(parameter_list$auc_min_pos_thresh>pos_thresh){pos_thresh = parameter_list$auc_min_pos_thresh}
#   if(parameter_list$auc_max_pos_thresh<pos_thresh){pos_thresh = parameter_list$auc_max_pos_thresh}
#   
#   # get cells:
#   cells_incelltype = metadata$Cell_ID[which(celltype_values>pos_thresh)]
#   all_mapped_cells[[current_celltype]] = cells_incelltype
#   
# }
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

## multiple auc matrix with geen set lengths (short signatures often ahve high scores which can distort comparison across signatures)
auc_mat_celltypes_mult = t(t(auc_mat_celltypes)*sapply(gene_set_list,length))
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
max_idx = data.frame(cell_id = rownames(merged_masks),cluster_id = as.character(apply(merged_masks,1,function(x){names(x)[x==1][1]}))) %>% na.omit()
# put into a list per cluster
all_mapped_cells = base::split(max_idx$cell_id,f=max_idx$cluster_id) # or: avg_log2fc or fc_mast

##########
### Save
##########

mapped_cells_file = parameter_list$mapped_celltypes_file
saveRDS(all_mapped_cells,file = mapped_cells_file)
mapped_cells_file = paste0(gsub(".txt|.rds|.json","",mapped_cells_file),".json")
writeList_to_JSON(list_with_rows = all_mapped_cells,filename = mapped_cells_file)

message(Sys.time(),": Finalized celltype mapping")


