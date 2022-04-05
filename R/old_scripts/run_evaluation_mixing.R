# run evaluation script
message(Sys.time(),": Starting evaluation" )

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
print(paramfile_name)

# read all parameters and filepaths
message("Reading parameters from: ",paramfile_name)
parameter_list = jsonlite::read_json(paramfile_name)

# read seurat object --> changed to metadata only
# seurat_object_evaluation <- SeuratDisk::LoadH5Seurat(parameter_list$seurat_object_path)
seurat_metadata_evaluation <- data.table::fread(parameter_list$seurat_object_path,data.table = F)
rownames(seurat_metadata_evaluation) = seurat_metadata_evaluation[,"Cell_ID"]

# read subsample ids
subset_cell_ids = data.table::fread(parameter_list$id_file_name,data.table = F,header = F)[,1]

#run mixing evaluation
eval_name = paste0(parameter_list$evaluation_file_path,"evaluation_mixing_prob.",parameter_list$eval_id,".",parameter_list$batch_var,".",parameter_list$ntrees_mixing,".",parameter_list$max_dim,".",parameter_list$global_seed,".txt")
evaluate_mixing(seurat_object_metadata = seurat_metadata_evaluation,
                integration_files = parameter_list$files_for_evaluation,
                integration_path = parameter_list$integration_res_path,
                evaluation_file=eval_name,
                batch_var=parameter_list$batch_var,
                subset_cell_ids=subset_cell_ids,
                ntrees=parameter_list$ntrees_mixing,
                ncores=parameter_list$ncores,
                max_dim = parameter_list$max_dim,
                scale_to_max=TRUE,
                max_for_norm=parameter_list$max_for_norm,
                global_seed=parameter_list$global_seed)

message(Sys.time(),": Finished evaluation" )


