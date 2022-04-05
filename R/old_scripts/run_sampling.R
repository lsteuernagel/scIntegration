# run evaluation script
message(Sys.time(),": Starting subsampling" )

# source functions
require(tidyverse)
require(Seurat)
require(Matrix)
require(SeuratDisk)
source("utils.R")

# get params from commandline
command_args<-commandArgs(TRUE)
paramfile_name = command_args[1]
print(paramfile_name)

# read all parameters and filepaths
message("Reading parameters from: ",paramfile_name)
parameter_list = jsonlite::read_json(paramfile_name)

# read metadata
seurat_metadata <- data.table::fread(parameter_list$seurat_metadata_path,data.table = F)
rownames(seurat_metadata) = seurat_metadata[,"Cell_ID"]
message(nrow(seurat_metadata))
# extract target sample size and set steps (or take from params)
target_sub_sample = parameter_list$target_sub_sample
message(target_sub_sample)
if(is.null(parameter_list$stepsize)){
  message("hi")
  stepsize = floor((nrow(seurat_metadata)-target_sub_sample) / 100)
}else{
  stepsize = parameter_list$stepsize
}
message("Running with step-size: ",stepsize)
# run sampling
tmp_res = downsample_balanced_iterative(seurat_metadata,global_seed=parameter_list$global_seed,
                                        predictor_var = parameter_list$predictor_var,stepsize = stepsize,target_sample_size = target_sub_sample)  
# get ids
sampling_ids= tmp_res[,"Cell_ID"]

#save in file
data.table::fwrite(data.table::as.data.table(sampling_ids),file = parameter_list$id_file_name,col.names = F)

# 
message("Finalized sampling")



