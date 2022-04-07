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

# load seurat
seurat_metadata_evaluation <- data.table::fread(parameter_list$seurat_merged_metadata,data.table = F)
rownames(seurat_metadata_evaluation) = seurat_metadata_evaluation[,"Cell_ID"]

##########
### Get all evaluation files
##########


# read subsample ids
subset_cell_ids = data.table::fread(parameter_list$id_file_name,data.table = F,header = F)[,1]

#run mixing evaluation

evaluate_mixing(seurat_object_metadata = seurat_metadata_evaluation,
                integration_files = parameter_list$files_for_evaluation,
                integration_path = parameter_list$integration_res_path,
                evaluation_file=parameter_list$evaluation_file,
                batch_var=parameter_list$batch_var,
                subset_cell_ids=subset_cell_ids,
                ntrees=parameter_list$ntrees_mixing,
                sampsize_pct=parameter_list$sampsize_pct,
                ncores=parameter_list$ncores,
                max_dim = parameter_list$max_dim,
                scale_to_max=TRUE,
                returnNormalized=TRUE,
                returnCollapsed =TRUE,
                max_for_norm=parameter_list$max_for_norm,
                global_seed=parameter_list$global_seed)

message(" Finalized ")


