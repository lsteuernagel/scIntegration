
##########
### Run evaluation: asw
##########

evaluation_folder = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/evaluation/")

all_eval_files = list.files(evaluation_folder,pattern="grouped")
all_eval_files = all_eval_files[grepl("_AuthorCellTypes",all_eval_files)]
table_list = list()
for(i in 1:length(all_eval_files)){
  label = gsub("AuthorCellTypes.csv_AuthorCellType_","",all_eval_files[[i]])
  label = gsub(".txt","",label)
  table_list[[label]] = data.table::fread(paste0(evaluation_folder,all_eval_files[[i]]))
  table_list[[label]]$Dataset = label
}
asw_grouped_author_celltype = do.call(rbind,table_list)
colnames(asw_grouped_author_celltype) = c("reduction","celltype","asw_eulid","asw_cosine","dataset")
asw_grouped_author_celltype$integration_type = "full"

##########
### Run evaluation: asw
##########

evaluation_folder = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_neurons_integration/evaluation/")

all_eval_files = list.files(evaluation_folder,pattern="grouped")
all_eval_files = all_eval_files[grepl("_AuthorCellTypes",all_eval_files)]
table_list = list()
for(i in 1:length(all_eval_files)){
  label = gsub("AuthorCellTypes.csv_AuthorCellType_","",all_eval_files[[i]])
  label = gsub(".txt","",label)
  table_list[[label]] = data.table::fread(paste0(evaluation_folder,all_eval_files[[i]]))
  table_list[[label]]$Dataset = label
}
asw_grouped_author_celltype_neurons = do.call(rbind,table_list)
colnames(asw_grouped_author_celltype_neurons) = c("reduction","celltype","asw_eulid","asw_cosine","dataset")
asw_grouped_author_celltype_neurons$integration_type = "neurons"

##########
### join / rbind
##########

asw_grouped_both = dplyr::bind_rows(asw_grouped_author_celltype %>% dplyr::select(-asw_cosine) %>% dplyr::filter(celltype %in% asw_grouped_author_celltype_neurons$celltype),
                                    asw_grouped_author_celltype_neurons %>% dplyr::select(-asw_cosine))

# add information
asw_grouped_both$method = stringr::str_extract(asw_grouped_both$reduction,pattern="[a-zA-Z]+_") %>% stringr::str_replace(pattern = "_",replacement = "")
# pca:
asw_grouped_both$ndim[asw_grouped_both$method=="PCA"]=asw_grouped_both$reduction[asw_grouped_both$method=="PCA"] %>% stringr::str_extract(pattern="\\.[0-9]+\\.") %>% stringr::str_replace(pattern = "\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.",replacement = "")
asw_grouped_both$features_ngenes[asw_grouped_both$method=="PCA"]=asw_grouped_both$reduction[asw_grouped_both$method=="PCA"] %>% stringr::str_extract(pattern="[0-9]+\\.[0-9]+") %>% stringr::str_replace(pattern = "[0-9]+\\.",replacement = "") #%>% stringr::str_replace(pattern = "\\.",replacement = "")
# scvi
asw_grouped_both$ndim[asw_grouped_both$method=="scVI"]=asw_grouped_both$reduction[asw_grouped_both$method=="scVI"] %>% stringr::str_extract(pattern="\\.\\..[0-9]+\\.\\.") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "")
asw_grouped_both$features_ngenes[asw_grouped_both$method=="scVI"] = stringr::str_extract(asw_grouped_both$reduction[asw_grouped_both$method=="scVI"],pattern="features\\.[0-9]+") %>% stringr::str_replace(pattern = "features\\.",replacement = "") %>% as.numeric()
asw_grouped_both$scvi_params[asw_grouped_both$method=="scVI"] = asw_grouped_both$reduction[asw_grouped_both$method=="scVI"] %>%
  stringr::str_extract(pattern="scVI_[0-9]+_[0-9]+_0\\.[0-9]+_[0-9]+_[0-9]+") %>% stringr::str_replace(pattern = "scVI_[0-9]+_",replacement = "") #%>% as.numeric()
asw_grouped_both$scvi_params[asw_grouped_both$method=="PCA"] = "PCA"

# set to nuemric
asw_grouped_both$asw_eulid = as.numeric(asw_grouped_both$asw_eulid)

##
a1 = asw_grouped_both %>% dplyr::distinct(integration_type,method,ndim,features_ngenes,scvi_params,.keep_all = TRUE) %>% dplyr::group_by()
a1 = a1 %>% dplyr::group_by(method,ndim,features_ngenes,scvi_params)


## a
a2a =  asw_grouped_both %>% dplyr::distinct(reduction,integration_type,method,ndim,features_ngenes,scvi_params)
a2 = asw_grouped_both %>% dplyr::group_by(reduction) %>% dplyr::summarise(mean_across_celltypes = mean(asw_eulid)) %>%
  dplyr::left_join(a2a,by="reduction")



##########
###  run basic comaprison neurons
##########

params_integration = jsonlite::read_json("data/parameters_integration_v2_neurons_1.json")
evaluation_folder = paste0(params_integration$integration_folder_path,"evaluation/")
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
colnames(evaluation_purityasw) =c(colnames(evaluation_purityasw)[2:length(colnames(evaluation_purityasw))],"placeholder")
evaluation_purityasw$reduction = evaluation_purityasw$placeholder

## join together for simple overview
evaluation_purityknn_mean = evaluation_purityknn %>% dplyr::group_by(reduction) %>% dplyr::summarise(mean_knn_purity = mean(purityknn))
all_evaluation_results = dplyr::full_join(evaluation_mixingrf,evaluation_mixingknn,by="reduction")
all_evaluation_results = dplyr::full_join(all_evaluation_results,evaluation_purityknn_mean,by="reduction")
all_evaluation_results = dplyr::full_join(all_evaluation_results,evaluation_purityasw[,c("reduction","silhouette_score_euclidean")],by="reduction")
all_evaluation_results = all_evaluation_results %>% dplyr::select(reduction, mixingrf, mixingknn , mean_knn_purity,asw = silhouette_score_euclidean)

##
asw_grouped_author_celltype_neurons_acrosscelltypes = asw_grouped_author_celltype_neurons %>% dplyr::group_by(reduction) %>%
  dplyr::summarise(asw_celltypes= mean(asw_eulid))
all_evaluation_results_extended =  dplyr::full_join(all_evaluation_results,asw_grouped_author_celltype_neurons_acrosscelltypes,by="reduction")

