# run evaluation script
message(Sys.time(),": Starting evaluation" )

# source functions
require(tidyverse)
require(Seurat)
require(Matrix)
require(SeuratDisk)
source("evaluation/evaluation_functions.R")
source("harmonization/seurat_annoy_functions.R")
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
# subset_cell_ids = data.table::fread(parameter_list$id_file_name,data.table = F,header = F)[,1]

# run purity evaluation
if(file.exists(parameter_list$mapped_celltypes_file)){
  loaded_cell_sets = readRDS(parameter_list$mapped_celltypes_file)
  # if(!is.null(parameter_list$max_cells_per_celltype)){
  #   mapped_celltypes = lapply(mapped_celltypes,)
  # }
  eval_name = paste0(parameter_list$evaluation_file_path,"evaluation_purity_knn.",parameter_list$eval_id,".",length(loaded_cell_sets),".",parameter_list$k_param,".",parameter_list$global_seed,".txt")
  evaluate_purity_knn(seurat_object_metadata = seurat_metadata_evaluation,
                      cells_sets = loaded_cell_sets,
                      k_param=parameter_list$k_param,
                      integration_files= parameter_list$files_for_evaluation,
                      integration_path=parameter_list$integration_res_path,
                      evaluation_file=eval_name,
                      ncores=parameter_list$ncores,
                      dist_type=parameter_list$dist_type,
                      max_dim = parameter_list$max_dim,
                      global_seed=parameter_list$global_seed)
  # seurat_object_metadata,cells_sets,k_param=20,integration_files,integration_path,evaluation_file,ncores=1,dist_type="cosine",max_dim=50,global_seed=123){
  
}else{
  warning("Cannot find file with mapped celltypes. Please provide valid path via 'mapped_celltypes_file'.")
}

message(Sys.time(),": Finished evaluation" )


