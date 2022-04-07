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
params_integration = jsonlite::read_json("data/parameters_integration_v2_1.json")
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
evaluation_folder = paste0(param_set$integration_folder_path,"evaluation/")
system(paste0("mkdir -p ",paste0(evaluation_folder)))

##########
### [6] Run evaluation: rf mixing
##########

files_per_batch = 10

## prepare mixing rf integration
evaluation_mixingrf_folder = paste0(param_set$integration_folder_path,"evaluation/mixing_rf/")
system(paste0("mkdir -p ",paste0(evaluation_mixingrf_folder)))
evaluation_mixingrf_file = paste0(evaluation_folder,"all_mixing_rf.txt")
dummy = data.frame(value="",reduction="")
data.table::fwrite(dummy,evaluation_mixingrf_file)

# A: Find which results to evaluate
# Read all files with integration results
all_integration_files=list.files(paste0(param_set$integration_folder_path,"integration/"),recursive = TRUE,pattern = ".txt")
available_integrations= gsub(".txt","",all_integration_files)
available_integrations = as.character(sapply(available_integrations,function(x){ strsplit(x,split = "/")[[1]][length(strsplit(x,split = "/")[[1]])]}))
# read all files with evaluation results
evaluation_mixingrf_table = data.table::fread(evaluation_mixingrf_file,data.table = FALSE)
available_evaluations= evaluation_mixingrf_table$reductions %>% na.omit() %>% as.character()
available_evaluations = as.character(sapply(available_evaluations,function(x){ strsplit(x,split = "/")[[1]][length(strsplit(x,split = "/")[[1]])]}))
#Determine which results have to evaluated
integrations_to_evaluate = setdiff(available_integrations,available_evaluations)

# determine number of jobs:
if(length(integrations_to_evaluate) > files_per_batch){
  cut_values = as.character(cut(1:length(integrations_to_evaluate),ceiling(length(integrations_to_evaluate)/files_per_batch)))
  cut_levels = unique(cut_values)
}else{
  cut_values = rep("1",length(all_files))
  cut_levels = unique(cut_values)
}

# B: start slurm jobs using result from A

# start one job per latent space set size
slurm_id_6_per_cut = vector()

# start one job per feature set size
for(i in 1:length(cut_levels)){

  # set parameters for current dataset
  param_set = params_integration
  # update with info for this evaluation
  param_set$integration_names = integrations_to_evaluate[which(cut_values==cut_levels[i])]
  param_set$integration_res_path = paste0(param_set$integration_folder_path,"integration/")
  param_set$evaluation_file = evaluation_mixingrf_file
  param_set$seurat_merged_metadata = gsub(".rds","_metadata.txt",param_set$merged_file)

  # make unique id:
  job_id=digest::digest(param_set)
  param_set$job_id = job_id
  # write to JSON as transfer file
  param_file = paste0(param_path,"evaluation_mixing_rf_params_",job_id,".json")
  writeList_to_JSON(list_with_rows = param_set,filename = param_file)
  # execute job
  script_path = "R/run_scripts/evaluate_mixing_rf.R"
  # set sbatch params:
  jobname = paste0("evaluation_mixing_rf","_",job_id)
  outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
  errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
  output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
  slurm_id_6_per_cut[i] = stringr::str_remove(output_message,pattern = "Submitted batch job ")
}
