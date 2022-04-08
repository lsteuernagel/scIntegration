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
seurat_merged = readRDS(paste0(parameter_list$merged_file))

#get metadata and save
seurat_metadata = seurat_merged@meta.data
data.table::fwrite(x = seurat_metadata,file = paste0(parameter_list$integration_folder_path,gsub(".rds","_metadata.txt",parameter_list$merged_file)),sep = "\t")

##########
### Export to anndata
##########
dummy=matrix(data = as.numeric())
seurat_merged@assays[["RNA"]]@var.features = character()
seurat_merged@assays[["RNA"]]@scale.data <- dummy[,-1] # error is okay

# make file name
merged_file_name = paste0(parameter_list$integration_folder_path,gsub(".rds","",parameter_list$merged_file))

# save h5seurat
SeuratDisk::SaveH5Seurat(object = seurat_merged,filename = paste0(merged_file_name,".h5seurat"), overwrite = TRUE, verbose = TRUE)

# save to anndata
SeuratDisk::Convert( paste0(merged_file_name,".h5seurat"), dest =  paste0(merged_file_name,".h5ad"),assay="RNA",verbose=TRUE,overwrite=TRUE)

message(" Complete ")




