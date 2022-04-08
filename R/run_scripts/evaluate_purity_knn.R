##########
### Load parameters and packages
##########

message(" Load parameters and packages ")

library(magrittr)
library(scUtils)

source("R/evaluation_functions.R")
source("R/annoy_from_seurat_functions.R")

# get params-filename from commandline
command_args<-commandArgs(TRUE)
param_file = command_args[1]
# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# load seurat meta
seurat_metadata_evaluation <- data.table::fread(parameter_list$seurat_merged_metadata,data.table = F)
rownames(seurat_metadata_evaluation) = seurat_metadata_evaluation[,"Cell_ID"]

# read cell sets
cells_sets = jsonlite::read_json(paste0(parameter_list$integration_folder_path,parameter_list$detected_cells_filename))
# if some fields are lists --> unlist
cells_sets = lapply(cells_sets,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})


##########
### Get all evaluation files
##########

#run mixing evaluation

#eval_name = paste0(parameter_list$evaluation_file_path,"evaluation_entropy_knn.",parameter_list$eval_id,".",length(loaded_cell_sets),".",parameter_list$k_param,".",parameter_list$global_seed,".txt")

evaluate_purity_knn(seurat_object_metadata = seurat_metadata_evaluation,
                    integration_names = parameter_list$integration_names,
                    integration_path = parameter_list$integration_res_path,
                    evaluation_file=parameter_list$evaluation_file,
                    cells_sets = cells_sets,
                    k_param=parameter_list$k_param,
                    ncores=parameter_list$ncores,
                    dist_type=parameter_list$dist_type,
                    global_seed=parameter_list$global_seed)

message(" Finalized ")


