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

#run purity evaluation
if(file.exists(parameter_list$perBatch_celltypes_file)){
  loaded_cell_sets = readRDS(parameter_list$perBatch_celltypes_file)
  eval_name = paste0(parameter_list$evaluation_file_path,"evaluation_purity_perBatch_cluster.",parameter_list$eval_id,".",length(loaded_cell_sets),".",parameter_list$global_seed,".txt")
  evaluate_purity_perBatch(seurat_object_metadata = seurat_metadata_evaluation,
                           cell_sets = loaded_cell_sets,
                           integration_files = parameter_list$files_for_evaluation,
                           integration_path = parameter_list$integration_res_path,
                           evaluation_file=eval_name,
                           ncores=parameter_list$ncores,
                           max_dim = parameter_list$max_dim,
                           global_seed=parameter_list$global_seed)
}
message(Sys.time(),": Finished evaluation" )



