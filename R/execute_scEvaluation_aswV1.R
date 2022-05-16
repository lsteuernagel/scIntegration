##########
### Load
##########

# source("R/integration_functions.R")
# source("R/evaluation_functions.R")
singularity_path = "~/Documents/r_scvi_015.simg"

# direct output and logs to some files on the local filesystem:
# where to store temporary json's with params for jobs:
param_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/hypoMap_v2_params/"
# where to save log files --> use this path in the slurm.sh files!
log_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/hypoMap_v2_slurmlogs/"

# load json file with all other information
params_integration = jsonlite::read_json("data/parameters_integration_v2_neurons_1.json")
# if some fields are lists --> unlist
params_integration = lapply(params_integration,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

### try to creat dir if necessary:
system(paste0("mkdir -p ",paste0(param_path)))
system(paste0("mkdir -p ",paste0(log_path)))

# define helper function
writeList_to_JSON = function (list_with_rows, filename){
  jsonfile = jsonlite::toJSON(list_with_rows, pretty = TRUE, auto_unbox = TRUE, digits = NA)
  writeLines(jsonfile, con = paste0(filename))
}

# make folder:
evaluation_folder = paste0(params_integration$integration_folder_path,"evaluation/")
system(paste0("mkdir -p ",paste0(evaluation_folder)))

require(magrittr)

##########
### Run evaluation: asw
##########

# this file is custom made: col1: Cell_ID, col2: label
label_file_name = paste0(param_set$integration_folder_path,"cell_labels_K169_v1.csv")

# define file names
evaluation_knownLabel_asw_file = paste0(evaluation_folder,"all_knownLabel_K169_asw.txt")
evaluation_knownLabel_asw_grouped_file = paste0(evaluation_folder,"grouped_knownLabel_K169_asw.txt")

# A: Find which results to evaluate
# Read all files with integration results
all_integration_files=list.files(paste0(params_integration$integration_folder_path,"integration/"),recursive = TRUE,pattern = ".txt")
available_integrations= gsub(".txt","",all_integration_files)
available_integrations = as.character(sapply(available_integrations,function(x){ strsplit(x,split = "/")[[1]][length(strsplit(x,split = "/")[[1]])]}))

# set parameters for current dataset
param_set = params_integration
# update with info for this evaluation
param_set$integration_names = available_integrations
param_set$integration_res_path = paste0(param_set$integration_folder_path,"integration/")
param_set$evaluation_file = evaluation_knownLabel_asw_file
param_set$evaluation_file_grouped = evaluation_knownLabel_asw_grouped_file
param_set$merged_file_h5ad = paste0(param_set$integration_folder_path,param_set$new_name_suffix,".h5ad")
param_set$label_file_name = label_file_name
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
### Run evaluation: asw
##########

# this file is custom made: col1: Cell_ID, col2: label
#label_folder ="/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/cell_labels/"
label_folder ="/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_neurons_integration/cell_labels/"
label_file_names = list.files(label_folder)

# A: Find which results to evaluate
# Read all files with integration results
all_integration_files=list.files(paste0(params_integration$integration_folder_path,"integration/"),recursive = TRUE,pattern = ".txt")
available_integrations= gsub(".txt","",all_integration_files)
available_integrations = as.character(sapply(available_integrations,function(x){ strsplit(x,split = "/")[[1]][length(strsplit(x,split = "/")[[1]])]}))

for(label_file in label_file_names){

  label_name = gsub("_AuthorCellTypes.tsv","",label_file)

  # define file names
  evaluation_knownLabel_asw_file = paste0(evaluation_folder,label_name,"_AuthorCellType_all_asw.txt")
  evaluation_knownLabel_asw_grouped_file = paste0(evaluation_folder,label_name,"_AuthorCellType_grouped_asw.txt")

  # set parameters for current dataset
  param_set = params_integration
  # update with info for this evaluation
  param_set$label_name = label_name
  param_set$integration_names = available_integrations
  param_set$integration_res_path = paste0(param_set$integration_folder_path,"integration/")
  param_set$evaluation_file = evaluation_knownLabel_asw_file
  param_set$evaluation_file_grouped = evaluation_knownLabel_asw_grouped_file
  param_set$merged_file_h5ad = paste0(param_set$integration_folder_path,param_set$new_name_suffix,".h5ad")
  param_set$label_file_name = paste0(label_folder,label_file)
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
