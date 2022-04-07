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

##########
###  Determine which results have to evaluated
##########

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



##########
### Save to file
##########



