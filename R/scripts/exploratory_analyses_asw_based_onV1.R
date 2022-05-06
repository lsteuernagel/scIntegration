
##########
### Export
##########

## folder_for_v1_integrations
folder_for_v1_integrations = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_integration/integration/"
system(paste0("mkdir -p ",paste0(folder_for_v1_integrations)))

# load hypomap v1 neurons
hypoMap_seurat_neurons_v1 = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_largeFiles/hypothalamus_neurons_map.rds")

# save embedding:
current_matrix =hypoMap_seurat_neurons_v1@reductions$scvi@cell.embeddings
current_matrix = cbind(Cell_ID=rownames(hypoMap_seurat_neurons_v1@meta.data),current_matrix)
data.table::fwrite(data.table::as.data.table(current_matrix),file=paste0(folder_for_v1_integrations,"scvi_hypoMap_neurons_v1.txt"),sep="\t",col.names = TRUE)


# load hypomap v1
hypoMap_seurat_v1 = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_largeFiles/hypothalamus_full_map.rds")

# save embedding:
current_matrix =hypoMap_seurat_v1@reductions$scvi@cell.embeddings
current_matrix = cbind(Cell_ID=rownames(hypoMap_seurat_v1@meta.data),current_matrix)
data.table::fwrite(data.table::as.data.table(current_matrix),file=paste0(folder_for_v1_integrations,"scvi_hypoMap_v1.txt"),sep="\t",col.names = TRUE)


## export label file:
# this file is custom made: col1: Cell_ID, col2: label
label_file_name = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_integration/cell_labels_K169_v1_1.csv")
cell_label_df = hypoMap_seurat_neurons_v1@meta.data[,c("Cell_ID","K169_pruned")]
# ensure that allc ell ids will work:
intersection_cellids=intersect(hypoMap_seurat_neurons_v1@meta.data$Cell_ID,hypoMap_seurat_v1@meta.data$Cell_ID) %>% na.omit()
cell_label_df = cell_label_df[cell_label_df$Cell_ID %in% intersect(hypoMap_seurat_neurons_v1@meta.data$Cell_ID,hypoMap_seurat_v1@meta.data$Cell_ID),]
data.table::fwrite(cell_label_df,file=label_file_name,col.names = TRUE)

### export new
#
# merged_file_h5ad = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_integration/hypothalamus_neurons_map.h5ad"
# # save h5seurat
# SeuratDisk::SaveH5Seurat(object = seurat_merged,filename = gsub("h5ad","h5seurat",merged_file_h5ad), overwrite = TRUE, verbose = TRUE)
# # save to anndata
# SeuratDisk::Convert( paste0(merged_file_name,".h5seurat"), dest =  merged_file_h5ad,assay="RNA",verbose=TRUE,overwrite=TRUE)


##########
### Define params
##########

# image
singularity_path = "~/Documents/r_scvi_015.simg"
# where to store temporary json's with params for jobs:
param_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/hypoMap_v2_params/"
# where to save log files --> use this path in the slurm.sh files!
log_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/hypoMap_v2_slurmlogs/"
# define helper function
writeList_to_JSON = function (list_with_rows, filename){
  jsonfile = jsonlite::toJSON(list_with_rows, pretty = TRUE, auto_unbox = TRUE, digits = NA)
  writeLines(jsonfile, con = paste0(filename))
}

# define integrations based on file names used above!:
available_integrations = c("scvi_hypoMap_neurons_v1","scvi_hypoMap_v1")

## folder_for_v1_integrations
folder_for_v1_integrations = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_integration/integration/"
system(paste0("mkdir -p ",paste0(folder_for_v1_integrations)))
folder_for_v1_evaluation = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_integration/evaluation/"
system(paste0("mkdir -p ",paste0(folder_for_v1_evaluation)))

# define file names
evaluation_knownLabel_asw_file = paste0(folder_for_v1_evaluation,"all_K169_asw.txt")
evaluation_knownLabel_asw_grouped_file = paste0(folder_for_v1_evaluation,"grouped_K169_asw.txt")

# make param list
param_set =list()
param_set$integration_folder_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/"#"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_integration/"
param_set$merged_file_h5ad = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_objects/hypothalamus_neurons_reference_RNA.h5ad"
param_set$integration_res_path = folder_for_v1_integrations
param_set$evaluation_file = evaluation_knownLabel_asw_file
param_set$evaluation_file_grouped = evaluation_knownLabel_asw_grouped_file
param_set$global_seed = 123456
param_set$integration_names = available_integrations
param_set$label_file_name = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_integration/cell_labels_K169_v1_1.csv") # see alos above --> hardcoded twice!!

##########
### Execute eval asw
##########

# make unique id:
job_id=digest::digest(param_set)
param_set$job_id = job_id
# write to JSON as transfer file
param_file = paste0(param_path,"evaluation_knownLabels_asw_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# execute job
script_path = "python/run_asw_evaluation_knownLabels.py"
# set sbatch params:
jobname = paste0("evaluation_knownLabels_asw","_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --kill-on-invalid-dep=yes python/run_Python_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
print(output_message)


##########
### Execute eval asw on author cell types
##########

# this file is custom made: col1: Cell_ID, col2: label
label_folder ="/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_integration/cell_labels/"
label_file_names = list.files(label_folder)

for(label_file in label_file_names){

  label_name = gsub("_AuthorCellTypes.csv","",label_file)

  param_set$label_file_name = paste0(paste0(label_folder,label_file)) # see alos above --> hardcoded twice!!
  param_set$label_name = label_name
  param_set$evaluation_file = paste0(folder_for_v1_evaluation,label_name,"_AuthorCellType_all_asw.txt")
  param_set$evaluation_file_grouped = paste0(folder_for_v1_evaluation,label_name,"_AuthorCellType_grouped_asw.txt")
  # make unique id:
  job_id=digest::digest(param_set)
  param_set$job_id = job_id
  # write to JSON as transfer file
  param_file = paste0(param_path,"evaluation_knownLabels_asw_params_",job_id,".json")
  writeList_to_JSON(list_with_rows = param_set,filename = param_file)
  # execute job
  script_path = "python/run_asw_evaluation_knownLabels.py"
  # set sbatch params:
  jobname = paste0("evaluation_knownLabels_asw","_",job_id)
  outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
  errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
  output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --kill-on-invalid-dep=yes python/run_Python_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
  print(output_message)

}

##########
### Read results
##########

all_eval_files = list.files(folder_for_v1_evaluation,pattern="all")
all_eval_files = all_eval_files[grepl("_AuthorCellTypes",all_eval_files)]
table_list = list()
for(i in 1:length(all_eval_files)){
  label = gsub("_AuthorCellTypes.csv","",all_eval_files[[i]])
  table_list[[label]] = data.table::fread(paste0(folder_for_v1_evaluation,all_eval_files[[i]]))
  table_list[[label]]$Dataset = label
}

a1 = do.call(rbind,table_list)

all_eval_files = list.files(folder_for_v1_evaluation,pattern="grouped")
all_eval_files = all_eval_files[grepl("_AuthorCellTypes",all_eval_files)]
table_list = list()
for(i in 1:length(all_eval_files)){
  label = gsub("_AuthorCellTypes.csv","",all_eval_files[[i]])
  table_list[[label]] = data.table::fread(paste0(folder_for_v1_evaluation,all_eval_files[[i]]))
  table_list[[label]]$Dataset = label
}

a2 = do.call(rbind,table_list)
colnames(a2) = c("reduction","celltype","asw_eulid","asw_cosine","dataset")

a2_wide = a2 %>% dplyr::select(-asw_cosine) %>% tidyr::spread(key="reduction",value="asw_eulid")
colnames(a2_wide)
ggplot(a2_wide,aes(x=scvi_hypoMap_v1,y=scvi_hypoMap_neurons_v1,color=dataset))+geom_point()+geom_abline(slope = 1)





