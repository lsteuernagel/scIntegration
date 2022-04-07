
##########
### Load parameters and packages
##########
message(" Load parameters and packages ")


library(magrittr)
library(scUtils)

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

# load seurat
seurat_merged = readRDS(paste0(parameter_list$data_path,parameter_list$merged_file))
seurat_metadata = seurat_merged@meta.data

##########
### Down sample
##########

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
tmp_res = scUtils::downsample_balanced_iterative(seurat_metadata,global_seed=parameter_list$global_seed,
                                        predictor_var = parameter_list$batch_var,stepsize = stepsize,target_sample_size = target_sub_sample)
# get ids
sampling_ids = list(cell_ids = tmp_res[,"Cell_ID"])

##########
### Save
##########

#save in file
writeList_to_JSON(list_with_rows = sampling_ids,filename = paste0(parameter_list$integration_folder_path,parameter_list$id_file_name))

#finalized
message("Finalized downsampling")



