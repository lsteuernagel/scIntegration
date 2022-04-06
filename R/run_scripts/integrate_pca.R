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
seurat_merged = readRDS(paste0(parameter_list$data_path,parameter_list$merged_file))

# load feature list
feature_set_list = jsonlite::read_json(feature_set_file)
# if some fields are lists --> unlist
feature_set_list = lapply(feature_set_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

##########
### Run PCA
##########

# which feature set to use:
features_to_use = feature_set_list[[paste0(parameter_list$assay_name,".log.", method, ".split_", batch_var, ".features.",parameter_list$feature_set_size)]]
messsage("Features to use: ",length(features_to_use))

# run PCA
seurat_merged = Seurat::RunPCA(seurat_merged,
                               assay = parameter_list$assay_name,
                               features = features_to_use,
                               npcs = max(as.numeric(parameter_list$latent_space_sizes)),
                               seed.use = parameter_list$global_seed)

# extract DimReduc
matrix_pca = seurat_merged@reductions[["PCA"]]@cell.embeddings
# make list for all desired sizes
list_with_PCAs = list()
for(size in as.numeric(parameter_list$latent_space_sizes)){
  list_with_PCAs[[size]] = matrix_pca[,1:size]
}

##########
### Save result
##########

# save all sizes
for(size in as.numeric(parameter_list$latent_space_sizes)){
  current_matrix_pca = list_with_PCAs[[size]]
  current_matrix_pca = cbind(Cell_ID=rownames(seurat_merged@meta.data),current_matrix_pca)
  filename = paste0(parameter_list$output_folder,"PCA_",param_list$batch_var,".",size,".",features_to_use,".txt")
  data.table::fwrite(data.table::as.data.table(current_matrix_pca),file=filename,sep="\t",col.names = TRUE)
}

message("Finalized PCA")





