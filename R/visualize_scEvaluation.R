
params_integration = jsonlite::read_json("data/parameters_integration_v2_3.json")
evaluation_folder = paste0(params_integration$integration_folder_path,"evaluation/")
evaluation_mixingrf_file = paste0(evaluation_folder,"all_mixing_rf.txt")
evaluation_purityknn_file = paste0(evaluation_folder,"all_purity_knn.txt")
evaluation_mixingknn_file = paste0(evaluation_folder,"all_mixing_knn.txt")
evaluation_purityasw_file = paste0(evaluation_folder,"all_purity_asw.txt")
require(magrittr)

data.table::fwrite(,file = evaluation_mixingrf_file)

## read files:
evaluation_mixingrf = data.table::fread(evaluation_mixingrf_file,data.table = F) %>% dplyr::rename( mixingrf = value)
evaluation_mixingknn = data.table::fread(evaluation_mixingknn_file,data.table = F) %>% dplyr::rename( mixingknn = value)
evaluation_purityknn = data.table::fread(evaluation_purityknn_file,data.table = F) %>% dplyr::rename( purityknn = value)
#evaluation_purityasw = data.table::fread(evaluation_purityasw_file,data.table = F)

## join together for simple overview
evaluation_purityknn_mean = evaluation_purityknn %>% dplyr::group_by(reduction) %>% dplyr::summarise(mean_knn_purity = mean(purityknn))
all_evaluation_results = dplyr::full_join(evaluation_mixingrf,evaluation_mixingknn,by="reduction")
all_evaluation_results = dplyr::full_join(all_evaluation_results,evaluation_purityknn_mean,by="reduction")
# all_evaluation_results = dplyr::full_join(all_evaluation_results,evaluation_purityasw,by="reduction")
all_evaluation_results = all_evaluation_results %>% dplyr::select(reduction, mixingrf, mixingknn , mean_knn_purity)

require(ggplot2)
ggplot2::ggplot(all_evaluation_results,aes(mixingrf,mixingknn))+geom_point()#+geom_abline(slope=1)

## norm scores
min_max = function(x,minVal=NULL,maxVal=NULL){
  if(is.null(minVal)){minVal =min(x,na.rm = TRUE) }
  if(is.null(maxVal)){maxVal =max(x,na.rm = TRUE) }
  return((x-minVal) / (maxVal-minVal))
}
all_evaluation_results$mean_purity_knn_norm = round(min_max(all_evaluation_results$mean_knn_purity,minVal = 0,maxVal = NULL),4)*100
all_evaluation_results$mixingknn_norm = round(min_max(all_evaluation_results$mixingknn,minVal = 0,maxVal = NULL),4)*100
all_evaluation_results$mixingrf_norm = round(min_max(all_evaluation_results$mixingrf,minVal = 0.8,maxVal = NULL),4)*100
all_evaluation_results$mixing_score = all_evaluation_results$mixingrf_norm*0.5+all_evaluation_results$mixingknn_norm*0.5
all_evaluation_results$purity_score =  all_evaluation_results$mean_purity_knn_norm#*0.75 + all_evaluation_results$silhouette_score_norm*0.25

# add information
all_evaluation_results$assay=all_evaluation_results$reduction %>% stringr::str_extract(pattern="\\.\\..[a-zA-Z]+\\.\\.") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "")
all_evaluation_results$ndim=all_evaluation_results$reduction %>% stringr::str_extract(pattern="\\.\\..[0-9]+\\.\\.") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "")
all_evaluation_results$features=all_evaluation_results$reduction %>% stringr::str_extract(pattern="\\.\\..*\\.features") %>% stringr::str_replace(pattern = "\\.\\..*\\.\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.features",replacement = "")
all_evaluation_results$method = stringr::str_extract(all_evaluation_results$reduction,pattern="[a-zA-Z]+_") %>% stringr::str_replace(pattern = "_",replacement = "")
all_evaluation_results$features_ngenes = stringr::str_extract(all_evaluation_results$reduction,pattern="features\\.[0-9]+") %>% stringr::str_replace(pattern = "features\\.",replacement = "") %>% as.numeric()
