##########
### [0] Load
##########

# source("R/integration_functions.R")
# source("R/evaluation_functions.R")
singularity_path = "~/Documents/r_scvi_v10.simg"

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


##########
### [1] Run feature detection
##########

# set params
param_set = params_integration
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"find_features_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "R/run_scripts/find_features.R"
# set sbatch params:
jobname = paste0("find_features_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = ""
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_1 = stringr::str_remove(output_message,pattern = "Submitted batch job ")


##########
### [2] Run celltype detection
##########

# set params
param_set = params_integration
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"detect_celltypes_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "R/run_scripts/detect_celltypes.R"
# set sbatch params:
jobname = paste0("detect_celltypes_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = ""
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_2 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

##########
### [3] Run downsampling
##########

# set params
param_set = params_integration
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"downsample_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "R/run_scripts/downsample_cells.R"
# set sbatch params:
jobname = paste0("downsample_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = ""
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_3 = stringr::str_remove(output_message,pattern = "Submitted batch job ")


##########
### [4] Run PCA integration + evaluation
##########

# make folder:
integration_pca_folder = paste0(param_set$integration_folder_path,"integration/pca/")
system(paste0("mkdir -p ",paste0(integration_pca_folder)))

# start one job per feature_set_sizes
slurm_id_4_per_size = vector()

for(i in 1:length(params_integration$feature_set_sizes)){

  #### Step 1: Integration:

  # set parameters for current dataset
  param_set = params_integration
  param_set$feature_set_size = params_integration$feature_set_sizes[i]
  param_set$output_folder = integration_pca_folder
  param_set$feature_set_file = paste0(parameter_list$integration_folder_path,"features/feature_sets.json")
  # make unique id:
  job_id=digest::digest(param_set)
  # write to JSON as transfer file
  param_file = paste0(param_path,"integration_pca_params_",job_id,".json")
  writeList_to_JSON(list_with_rows = param_set,filename = param_file)
  # execute job
  script_path = "R/run_scripts/integrate_pca.R"
  # set sbatch params:
  jobname = paste0("integration_pca_",dataset_table$Dataset[i],"_",job_id)
  outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
  errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
  dependency_ids = slurm_id_1
  output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",dependency_ids," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
  slurm_id_integration = stringr::str_remove(output_message,pattern = "Submitted batch job ")

  #### Step 2a: evaluation (R)

  param_set = params_integration
  param_set$ndims = params_integration$latent_space_sizes[i]
  param_set$input_integration = ""

  # make unique id:
  job_id=digest::digest(param_set)
  # write to JSON as transfer file
  param_file = paste0(param_path,"integration_pca_params_",job_id,".json")
  writeList_to_JSON(list_with_rows = param_set,filename = param_file)
  # execute job
  script_path = "R/run_scripts/integrate_pca.R"


  #### Step 2b: evaluation (python)

}

# Integration depends on 1

# Evaluation depends on 2+3 and coresponding integration run

##########
### [5] Run scVI integration + evaluation
##########

# make folder:
integration_scvi_folder = paste0(param_set$integration_folder_path,"integration/scvi/")
system(paste0("mkdir -p ",paste0(integration_scvi_folder)))

# start one job per latent space set size
slurm_id_5_per_size = vector()

# start one job per feature set size (and latent space ?)

# require param df

# Integration depends on 1

# Evaluation depends on 2+3 and coresponding integration run



