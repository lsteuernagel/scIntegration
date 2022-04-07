##########
### [6] Run evaluation: rf mixing
##########

# Run other script first

# R: Evaluation depends on 2+3, slurm_id_4_per_size and slurm_id_5_per_size

# A: Check script: Returns json with all files not in evaluation folder but in integration folder and estimates a number of jobs to start

# read all files with integration results
all_integration_files=list.files(parameter_list$integration_res_path,recursive = TRUE,pattern = ".txt")
available_integrations= gsub(".txt","",all_integration_files)
available_integrations = as.character(sapply(available_integrations,function(x){ strsplit(x,split = "/")[[1]][length(strsplit(x,split = "/")[[1]])]}))
# read all files with evaluation results
all_evaluation_files=list.files(parameter_list$evaluation_res_path,recursive = TRUE,pattern = ".txt")
available_evaluations= gsub(".txt","",all_evaluation_files)
available_evaluations = as.character(sapply(available_evaluations,function(x){ strsplit(x,split = "/")[[1]][length(strsplit(x,split = "/")[[1]])]}))
#Determine which results have to evaluated
integrations_to_evaluate = setdiff(available_integrations,available_evaluations)

# B: start slurm jobs using result from A


param_set$evaluation_file = paste0(param_set$evaluation_file_path,"evaluation_mixing_rf.",
                                   param_set$batch_var,".",
                                   param_set$ntrees_mixing,".",
                                   param_set$sampsize_pct,".", ".txt")

# C: Checks evaluation folder and saves all files in there in merged table. depends on finishing of B
