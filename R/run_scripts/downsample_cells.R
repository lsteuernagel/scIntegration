# subsample_param_file = list(seurat_metadata_path=seurat_metadata_path,
#                             id_file_name=id_file_name,
#                             stepsize=100,
#                             target_sub_sample=target_sub_sample,
#                             predictor_var="Dataset",#batch_var,
#                             global_seed=global_seed)
# hash = digest::digest(subsample_param_file)
# ## save parameters
# subsample_param_path = paste0(scHarmonize_path,key,"/",script_params,"/","subsample_parameter_file",".json")
# writeList_to_JSON(list_with_rows = subsample_param_file,filename = subsample_param_path)
# ## execute bash script or slurm bash scipt to run the job
# script_path = "evaluation/run_sampling.R"
# # run
# if(run_slurm_job){
#   system(paste0("sbatch bash/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",subsample_param_path))
# }else{
#   job_path = paste0(scHarmonize_path,key,"/",log_files,"/","subsample_",hash,".log")
#   system(paste0("bash/run_Rscript.sh ",singularity_path," ",script_path," ",subsample_param_path," ",job_path))
# }
#
# ## check:
# meta_data = data.table::fread(seurat_metadata_path,data.table = F)
# subsample_ids = data.table::fread(id_file_name,data.table = F,header = F)[,1]
# table(meta_data$Dataset[meta_data$Cell_ID %in% subsample_ids])

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
writeList_to_JSON(list_with_rows = sampling_ids,filename = paste0(integration_folder_path,parameter_list$detected_cells_filename))

#finalized
message("Finalized downsampling")



