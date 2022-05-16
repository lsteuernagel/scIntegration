##########
### function to laod eval results
##########

load_evaluation_results_withSCVI = function(param_file,evaluation_folder = paste0(params_integration$integration_folder_path,"evaluation/"),cell_types_purity=NULL){

  params_integration = jsonlite::read_json(param_file)

  evaluation_mixingrf_file = paste0(evaluation_folder,"all_mixing_rf.txt")
  evaluation_purityknn_file = paste0(evaluation_folder,"all_purity_knn.txt")
  evaluation_mixingknn_file = paste0(evaluation_folder,"all_mixing_knn.txt")
  evaluation_purityasw_file = paste0(evaluation_folder,"all_purity_asw.txt")
  require(magrittr)

  ## read files:
  evaluation_mixingrf = data.table::fread(evaluation_mixingrf_file,data.table = F) %>% dplyr::rename( mixingrf = value)
  evaluation_mixingknn = data.table::fread(evaluation_mixingknn_file,data.table = F) %>% dplyr::rename( mixingknn = value)
  evaluation_purityknn = data.table::fread(evaluation_purityknn_file,data.table = F) %>% dplyr::rename( purityknn = value)
  colnames(evaluation_purityknn) = c("celltype", "reduction", "purityknn")
  evaluation_purityasw = data.table::fread(evaluation_purityasw_file,data.table = F) #%>% dplyr::rename( asw = value)
  if(colnames(evaluation_purityasw)[2] == "reduction"){
    colnames(evaluation_purityasw) =c(colnames(evaluation_purityasw)[2:length(colnames(evaluation_purityasw))],"placeholder")
    evaluation_purityasw$reduction = evaluation_purityasw$placeholder
    evaluation_purityasw = evaluation_purityasw %>% dplyr::select(- placeholder)
  }

  ## join together for simple overview
  if(is.null(cell_types_purity)){cell_types_purity=unique(evaluation_purityknn$celltype)}
  evaluation_purityknn_mean = evaluation_purityknn %>% dplyr::filter(celltype %in% cell_types_purity) %>% dplyr::group_by(reduction) %>% dplyr::summarise(mean_knn_purity = mean(purityknn))
  all_evaluation_results = dplyr::full_join(evaluation_mixingrf,evaluation_mixingknn,by="reduction")
  all_evaluation_results = dplyr::full_join(all_evaluation_results,evaluation_purityknn_mean,by="reduction")
  all_evaluation_results = dplyr::full_join(all_evaluation_results,evaluation_purityasw[,c("reduction","silhouette_score_euclidean")],by="reduction")
  all_evaluation_results = all_evaluation_results %>% dplyr::select(reduction, mixingrf, mixingknn , mean_knn_purity,asw = silhouette_score_euclidean)

  # add information
  all_evaluation_results$method = stringr::str_extract(all_evaluation_results$reduction,pattern="[a-zA-Z]+_") %>% stringr::str_replace(pattern = "_",replacement = "")
  # pca:
  all_evaluation_results$ndim[all_evaluation_results$method=="PCA"]=all_evaluation_results$reduction[all_evaluation_results$method=="PCA"] %>% stringr::str_extract(pattern="\\.[0-9]+\\.") %>% stringr::str_replace(pattern = "\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.",replacement = "")
  all_evaluation_results$features_ngenes[all_evaluation_results$method=="PCA"]=all_evaluation_results$reduction[all_evaluation_results$method=="PCA"] %>% stringr::str_extract(pattern="[0-9]+\\.[0-9]+") %>% stringr::str_replace(pattern = "[0-9]+\\.",replacement = "") #%>% stringr::str_replace(pattern = "\\.",replacement = "")
  # scvi
  all_evaluation_results$ndim[all_evaluation_results$method=="scVI"]=all_evaluation_results$reduction[all_evaluation_results$method=="scVI"] %>% stringr::str_extract(pattern="\\.\\..[0-9]+\\.\\.") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "")
  all_evaluation_results$features_ngenes[all_evaluation_results$method=="scVI"] = stringr::str_extract(all_evaluation_results$reduction[all_evaluation_results$method=="scVI"],pattern="features\\.[0-9]+") %>% stringr::str_replace(pattern = "features\\.",replacement = "") %>% as.numeric()
  all_evaluation_results$scvi_params[all_evaluation_results$method=="scVI"] = all_evaluation_results$reduction[all_evaluation_results$method=="scVI"] %>%
    stringr::str_extract(pattern="scVI_[0-9]+_[0-9]+_0\\.[0-9]+_[0-9]+_[0-9]+") %>% stringr::str_replace(pattern = "scVI_[0-9]+_",replacement = "") #%>% as.numeric()
  all_evaluation_results$scvi_params[all_evaluation_results$method=="PCA"] = "PCA"
  all_evaluation_results$nlayers = all_evaluation_results$scvi_params %>% stringr::str_extract(pattern = "_[0-9]{1,}_")%>% stringr::str_replace(pattern = "_",replacement = "") %>% stringr::str_replace(pattern = "_",replacement = "")
  all_evaluation_results$nlayers = all_evaluation_results$scvi_params %>% stringr::str_extract(pattern = "_[0-9]{1,}_")%>% stringr::str_replace(pattern = "_",replacement = "") %>% stringr::str_replace(pattern = "_",replacement = "")
  all_evaluation_results$max_epochs  = all_evaluation_results$scvi_params %>% stringr::str_extract(pattern = "[0-9]+_0\\.")%>% stringr::str_replace(pattern = "_0\\.",replacement = "") #%>% stringr::str_replace(pattern = "_",replacement = "")
  all_evaluation_results$dropout_rate = all_evaluation_results$scvi_params %>% stringr::str_extract(pattern = "_0\\.[0-9]+")%>% stringr::str_replace(pattern = "_",replacement = "") #%>% stringr::str_replace(pattern = "_",replacement = "")

  return(all_evaluation_results)
}


##########
###  run basic comaprison neurons
##########


# all_K169_asw_v1 = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_integration/evaluation/all_K169_asw.txt",data.table = F)
# grouped_K169_asw_v1 = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_integration/evaluation/grouped_K169_asw.txt",data.table = F)

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

##########
###  run basic comaprison neurons
##########


cell_types_to_include = c('Agrp_Acvr1c','Fst_Fezf2','Ghrh_Mbnl3','Gnrh1_Gng8','Gpr50_Pgr15l','Hcrt_Rfx4','Lef1_Wif1','Npw_Nkx24','Ntrk1_Chat','Omp_Sim2','Oxt_Avp','Pirt_Tbx19','Pmch_Parpbp','Pomc_Anxa2','Pomc_Ttr','Prph_Hdc','Qrfp_Car8','Shox2_Gbx2','Slc6a3_Rxfp2','Sln_Ctxn3','Tac2_Kiss1','Ttn_Foxb1','Vip_Grp','Vip_Nov')

all_evaluation_results_full = load_evaluation_results_withSCVI(param_file = "data/parameters_integration_v2_3.json",cell_types_purity=cell_types_to_include)
all_evaluation_results_neurons = load_evaluation_results_withSCVI(param_file = "data/parameters_integration_v2_neurons_1.json",cell_types_purity=cell_types_to_include)


compare_eval = dplyr::full_join(all_evaluation_results_full,all_evaluation_results_neurons,suffix=c("_full","_neuron"),by=c("ndim" = "ndim", "features_ngenes" = "features_ngenes", "scvi_params"= "scvi_params"))
compare_eval2 = compare_eval[compare_eval$method_full=="scVI" & !is.na(compare_eval$reduction_full) & !is.na(compare_eval$reduction_neuron),]

### check on others

# define file names
evaluation_folder = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_neurons_integration/evaluation/"
label_name = "Affinati10x_AuthorCellTypes.csv"
evaluation_knownLabel_asw_file = paste0(evaluation_folder,label_name,"_AuthorCellType_all_asw.txt")
evaluation_knownLabel_asw_grouped_file = paste0(evaluation_folder,label_name,"_AuthorCellType_grouped_asw.txt")

evaluation_knownLabel_asw = data.table::fread(evaluation_knownLabel_asw_file,data.table = F)


