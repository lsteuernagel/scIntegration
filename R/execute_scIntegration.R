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

##########
### [0] prepare
##########

# set params
param_set = params_integration
# make unique id:
job_id=digest::digest(param_set)
# write to JSON as transfer file
param_file = paste0(param_path,"prepare_integration_params_",job_id,".json")
writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "R/run_scripts/prepare_integration.R"
# set sbatch params:
jobname = paste0("prepare_integration_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = ""
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_0 = stringr::str_remove(output_message,pattern = "Submitted batch job ")

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
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
slurm_id_1 = stringr::str_remove(output_message,pattern = "Submitted batch job ")


##########
### [2] Run celltype detection
##########

# set params
param_set = params_integration
#param_set$mapped_celltypes_auc_file
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
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
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
output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
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
  # set parameters for current dataset
  param_set = params_integration
  param_set$feature_set_size = params_integration$feature_set_sizes[i]
  param_set$output_folder = integration_pca_folder
  param_set$feature_set_file = paste0(params_integration$integration_folder_path,"features/feature_sets.json")
  # make unique id:
  job_id=digest::digest(param_set)
  # write to JSON as transfer file
  param_file = paste0(param_path,"integration_pca_params_",job_id,".json")
  writeList_to_JSON(list_with_rows = param_set,filename = param_file)
  # execute job
  script_path = "R/run_scripts/integrate_pca.R"
  # set sbatch params:
  jobname = paste0("integration_pca_",params_integration$feature_set_sizes[i],"_",job_id)
  outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
  errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
  dependency_ids = c(slurm_id_0,slurm_id_1)
  output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes R/run_scripts/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
  slurm_id_4_per_size[i] = stringr::str_remove(output_message,pattern = "Submitted batch job ")
}

# Integration depends on 1

##########
### [5] Run scVI integration
##########

library(magrittr)
n_models_to_test = 5

# make folder:
integration_scvi_folder = paste0(param_set$integration_folder_path,"integration/scvi/")
system(paste0("mkdir -p ",paste0(integration_scvi_folder)))

# start one job per latent space set size
slurm_id_5_per_size = vector()

# start one job per feature set size
for(i in 1:length(params_integration$feature_set_sizes)){

  # require param df
  require(magrittr)
  # read scvi params from json
  scVI_fullarg_list = jsonlite::read_json("data/parameters_scvi_1.json")
  # make a param data frame
  scvi_full_param_df = scVI_fullarg_list %>% purrr::cross_df()
  #subsample to random subset
  set.seed(params_integration$global_seed)
  param_df_random = scvi_full_param_df[sample(nrow(scvi_full_param_df),n_models_to_test),]
  # save this to the jobfile and use in job
  param_df_id=digest::digest(param_df_random)
  data.table::fwrite(param_df_random,paste0(param_path,param_df_id,"paramd_df.txt"),sep="\t")

  # set parameters for current dataset
  param_set = params_integration
  param_set$feature_set_size = params_integration$feature_set_sizes[i]
  param_set$output_folder = integration_scvi_folder
  param_set$data_filepath_full = gsub(".rds",".h5ad",param_set$merged_file)
  param_set$feature_set_file = paste0(params_integration$integration_folder_path,"features/feature_sets.json")
  param_set$hvgs_set_name = paste0(params_integration$assay_name,".log.", "vst", ".split_", params_integration$batch_var, ".features.",params_integration$feature_set_sizes[i])
  param_set$categorical_covariates = c("Dataset",param_set$batch_var) #character(0) # c("inferred_sex"), # need to pass batch var here !
  param_set$continuous_covariates =character(0)
  param_set$use_cuda =FALSE
  param_set$param_path = paste0(param_path,param_df_id,"paramd_df.txt")
  # make unique id:
  job_id=digest::digest(param_set)
  param_set$job_id = job_id
  # write to JSON as transfer file
  param_file = paste0(param_path,"integration_scvi_params_",job_id,".json")
  writeList_to_JSON(list_with_rows = param_set,filename = param_file)
  # execute job
  script_path = "python/run_scVI_v015.py"
  # set sbatch params:
  jobname = paste0("integration_scvi_", params_integration$feature_set_sizes[i],"_",job_id)
  outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
  errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
  dependency_ids = c(slurm_id_0,slurm_id_1)
  output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",paste0(dependency_ids,collapse = ":")," --kill-on-invalid-dep=yes python/run_Python_slurm.sh ",singularity_path," ",script_path," ",param_file),intern = TRUE)
  slurm_id_5_per_size[i] = stringr::str_remove(output_message,pattern = "Submitted batch job ")
}







