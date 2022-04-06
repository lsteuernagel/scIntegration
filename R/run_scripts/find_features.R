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
### Normalized and run feature detection
##########

message(" Add variable features ")

# normalize data
merged_seurat <- Seurat::NormalizeData(object = merged_seurat,  verbose = F, assay = "RNA")

# find HVGs
merged_seurat = scUtils::identify_variable_features(merged_seurat,
                                                    n_hvgs_sizes = parameter_list$feature_set_sizes,
                                                    batch_var = parameter_list$sample_column,
                                                    assay_name = "RNA",
                                                    method = "vst",
                                                    ignore_genes_vector = features_exclude_list,
                                                    returnSeurat = TRUE,
                                                    seed = parameter_list$global_seed)

feature_sets = merged_seurat@misc$var_features

##########
### Save
##########

feature_set_path = paste0(parameter_list$integration_folder_path,"features/")
system(paste0("mkdir -p ",paste0(feature_set_path)))

scUtils::writeList_to_JSON(feature_sets,filename = paste0(feature_set_path,"feature_sets.json"))



