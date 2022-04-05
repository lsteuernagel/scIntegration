# run harmony script
message(Sys.time(),": Starting harmony execution." )

# source functions
require(tidyverse)
require(Seurat)
require(Matrix)
require(SeuratDisk)
source("evaluation/evaluation_functions.R")
#source("processing/export_seurat.R")
source("integration/harmony_functions.R")
source("utils.R")


# get params from commandline
command_args<-commandArgs(TRUE)
paramfile_name = command_args[1]
message(paramfile_name)
#paramfile_name = paste0(scHarmonize_path,key,"/",script_params,"/","harmony_parameter_file_",harmony_param_file$job_id,".RDS")

# read all parameters and filepaths
parameter_list = readRDS(file=paramfile_name)
#parameter_list = jsonlite::read_json(paramfile_name)
# with json: unlist some entries
# which_transform = c("")
# parameter_list[names(parameter_list) %in% which_transform] = lapply(parameter_list[names(parameter_list) %in% which_transform],function(x){return(unlist(x))})

# make sure output dir exists
system(paste0("mkdir -p ",paste0(parameter_list$result_path)))

#I only load a subset of the object, by default trying counts from RNA
# but if I don't provide RNA i can specify another assay in the param slot
if(is.null(parameter_list$assay_slot)){parameter_list$assay_slot=list(RNA=c("counts"))}

# read seurat object
# change to metadata only ? --> for now just load a lighter object but with all reductions, then here we don't need handling to load all PCAs from file
message("Loading from: ",parameter_list$seurat_object_path)
seurat_object_harmony <- SeuratDisk::LoadH5Seurat(parameter_list$seurat_object_path,assays=parameter_list$assay_slot)

# add reduction
if(parameter_list$reduction_name %in% names(seurat_object_harmony@reductions)){
  message("Found reduction in seurat object reductions slot. Will use this one")  
}else{
  message("Did not find reduction in seurat object. Loading from: ",parameter_list$pca_path)  
  seurat_object_harmony = add_reduction_seurat(seurat_object_harmony,integration_name=parameter_list$reduction_name ,integration_path=parameter_list$pca_path,
                                               max_dim=parameter_list$dims_to_use,calc_umap=FALSE,overwrite =TRUE,overwrite2=FALSE)
}

# I am creating ids in the param dataframe
#full_param_df = data.table::fread(parameter_list$paramdf_path,data.table = F)
full_param_df = parameter_list$param_df
for(i in 1:nrow(full_param_df)){
  # get param set into list
  arguments_as_list = as.list(full_param_df[i,])
  arguments_as_list$sigma = unlist(as.data.frame(arguments_as_list$sigma))
  arguments_as_list$theta = unlist(as.data.frame(arguments_as_list$theta))
  arguments_as_list$vars_use = as.character(unlist(as.data.frame(arguments_as_list$vars_use)))
  full_param_df$full_id[i] = paste0("harmony_",i,"_",paste0(arguments_as_list[["vars_use"]],collapse = "."),"_",paste0(arguments_as_list[["theta"]],collapse = "."),"_",paste0(arguments_as_list[["sigma"]],collapse = "."),"_",parameter_list$nclust_K,"_",parameter_list$reduction_name,"_",parameter_list$job_id)
}

######## run harmony ######
run_harmony(seurat_object=seurat_object_harmony,reduction_name=parameter_list$reduction_name,param_df=full_param_df,dims_to_use = parameter_list$dims_to_use,
            kmeans_nstart=parameter_list$kmeans_nstart,nclust_K=parameter_list$nclust_K,kmeans_celltypeMin=parameter_list$kmeans_celltypeMin,
            filepath=parameter_list$result_path,global_seed=parameter_list$global_seed,random_id=parameter_list$job_id)

# add some additional info to param df before saving
full_param_df$reduction = parameter_list$dims_to_use
full_param_df$kmeans_nstart = parameter_list$kmeans_nstart
full_param_df$global_seed = parameter_list$global_seed


## Save parameters
param_file_full_path = paste0(parameter_list$result_path,"seuratIntegration_",parameter_list$assay,"_",parameter_list$random_id,"_parameters.rds")
# save updated files
message("Saving parameters to RDS file. Finished execution.")
saveRDS(full_param_df,param_file_full_path)


#### some old code I might need to build a function that assemble all parameter files

# check if file exists, filename contains length(cells_sets) ntrees and the seed used, so if this changes a new file will be created
# if not create empty dataframe(?) RDS and potentially path
# if(!file.exists(param_file_full)){
#   message("Did not find parameter file @ ",param_file_full," || Creating new RDS and directoy if necessary")  
#   system(paste0("mkdir -p ",paste0(filepath)))
#   all_parameters_harmony = data.frame()
#   saveRDS(all_parameters_harmony, param_file_full)
# }
# # load current results and extract names (only based on normal prob)
# all_parameters_harmony = readRDS(param_file_full)
# 
# # add details and input names to param file
# if(nrow(all_parameters_harmony)==0){
#   all_parameters_harmony = param_df_print
# }else{
#   all_parameters_harmony = rbind(all_parameters_harmony,param_df_print)
# }
# # save updated files
# message("Saving parameters to updated RDS file.")
# saveRDS(all_parameters_harmony,param_file_full)




