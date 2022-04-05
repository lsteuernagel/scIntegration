##########
### [0] Load
##########

source("R/functions.R")
singularity_path = "~/Documents/r_scvi_v10.simg"

# direct output and logs to some files on the local filesystem:
# where to store temporary json's with params for jobs:
param_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/hypoMap_v2_params/"
# where to save log files --> use this path in the slurm.sh files!
log_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/hypoMap_v2_slurmlogs/"

# load json file with all other information
params_pre_processing = jsonlite::read_json("data/parameters_integration_v2_1.json")
# if some fields are lists --> unlist
params_pre_processing = lapply(params_pre_processing,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

### try to creat dir if necessary:
system(paste0("mkdir -p ",paste0(param_path)))
system(paste0("mkdir -p ",paste0(log_path)))


##########
### [1] Run feature detection
##########


##########
### [2] Run celltype detection
##########

slurm_id_1_per_dataset = vector()

for(i in 1:nrow(dataset_table)){
  # set parameters for current dataset
  param_set = params_pre_processing
  param_set$dataset_name = dataset_table$Dataset[i]
  param_set$raw_file = dataset_table$seurat_file[i]
  param_set$target_cluster_number = dataset_table$estimated_total_clusters[i] # need to make a good estimate
  param_set$exclude_author = dataset_table$exclude_author[i]
  # make unique id:
  job_id=digest::digest(param_set)
  # write to JSON as transfer file
  param_file = paste0(param_path,"preprocessing_params_",job_id,".json")
  writeList_to_JSON(list_with_rows = param_set,filename = param_file)
  # execute job
  script_path = "R/run_scripts/preprocessing.R"
  # set sbatch params:
  jobname = paste0("preprocessing_",dataset_table$Dataset[i],"_",job_id)
  outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
  errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
  #run job
  output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
  slurm_id_1_per_dataset[dataset_table$Dataset[i]] = stringr::str_remove(output_message,pattern = "Submitted batch job ")
}



##########
### [3] Run downsampling
##########




##########
### [4] Run PCA integration + evaluation
##########

# Integration depends on 1

# Evaluation depends on 2+3

##########
### [5] Run scVI integration + evaluation
##########





