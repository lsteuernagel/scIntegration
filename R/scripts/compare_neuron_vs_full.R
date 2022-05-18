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
###  results from asw per dataset function
##########


# all_K169_asw_v1 = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_integration/evaluation/all_K169_asw.txt",data.table = F)
# grouped_K169_asw_v1 = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_integration/evaluation/grouped_K169_asw.txt",data.table = F)

load_evaluation_results_dataset_asw = function(param_file,evaluation_folder = paste0(params_integration$integration_folder_path,"evaluation/"),type="all",clean=TRUE){
  # or type = "grouped"

  params_integration = jsonlite::read_json(param_file)
  all_eval_files = list.files(evaluation_folder,pattern=type)
  if(clean){
    all_eval_files = all_eval_files[grepl("_AuthorCellTypes",all_eval_files)]
  }
  table_list = list()
  for(i in 1:length(all_eval_files)){
    label = gsub("_AuthorCellTypes.csv","",all_eval_files[[i]])
    table_list[[label]] = data.table::fread(paste0(evaluation_folder,all_eval_files[[i]]))
    table_list[[label]]$Dataset = label
  }

  all_asw_results = do.call(rbind,table_list)

  if(type=="grouped"){colnames(all_asw_results)[1:4] = c("reduction","celltype","asw_euclidean","asw_cosine")}else{colnames(all_asw_results)[1:3] = c("reduction","asw_euclidean","asw_cosine")}
  return(all_asw_results)
}

##########
###  run basic comaprison neurons
##########

# use only these cell types for üurity knn:
cell_types_to_include = c('Agrp_Acvr1c','Fst_Fezf2','Ghrh_Mbnl3','Gnrh1_Gng8','Gpr50_Pgr15l','Hcrt_Rfx4','Lef1_Wif1','Npw_Nkx24','Ntrk1_Chat','Omp_Sim2','Oxt_Avp','Pirt_Tbx19','Pmch_Parpbp','Pomc_Anxa2','Pomc_Ttr','Prph_Hdc','Qrfp_Car8','Shox2_Gbx2','Slc6a3_Rxfp2','Sln_Ctxn3','Tac2_Kiss1','Ttn_Foxb1','Vip_Grp','Vip_Nov')

# load
all_evaluation_results_full = load_evaluation_results_withSCVI(param_file = "data/parameters_integration_v2_3.json",cell_types_purity=cell_types_to_include)
all_evaluation_results_neurons = load_evaluation_results_withSCVI(param_file = "data/parameters_integration_v2_neurons_1.json",cell_types_purity=cell_types_to_include)

### check on others
per_dataset_asw_full_all = load_evaluation_results_dataset_asw(param_file = "data/parameters_integration_v2_3.json")
per_dataset_asw_full_grouped = load_evaluation_results_dataset_asw(param_file = "data/parameters_integration_v2_3.json",type = "grouped")
per_dataset_asw_neurons_all = load_evaluation_results_dataset_asw(param_file = "data/parameters_integration_v2_neurons_1.json")
per_dataset_asw_neurons_grouped = load_evaluation_results_dataset_asw(param_file = "data/parameters_integration_v2_neurons_1.json",type = "grouped")

## use grouped data and subset full to same ones as neuron:
per_dataset_asw_full_grouped_subset = per_dataset_asw_full_grouped[per_dataset_asw_full_grouped$celltype %in% unique(per_dataset_asw_neurons_grouped$celltype),]

# make per reduction summaries:
per_dataset_asw_full_grouped_subset_summary = per_dataset_asw_full_grouped_subset %>% dplyr::group_by(reduction) %>% dplyr::summarise(mean_asw_dataset_celltype = mean(asw_euclidean))
per_dataset_asw_neurons_grouped_summary = per_dataset_asw_neurons_grouped %>% dplyr::group_by(reduction) %>% dplyr::summarise(mean_asw_dataset_celltype = mean(asw_euclidean))

# add other data
all_evaluation_results_full = dplyr::left_join(all_evaluation_results_full,per_dataset_asw_full_grouped_subset_summary,by="reduction")
all_evaluation_results_neurons = dplyr::left_join(all_evaluation_results_neurons,per_dataset_asw_neurons_grouped_summary,by="reduction")

# join everything:
compare_eval = dplyr::full_join(all_evaluation_results_full,all_evaluation_results_neurons[!grepl("fromFull",all_evaluation_results_neurons$reduction),],suffix=c("_full","_neuron"),by=c("ndim" = "ndim", "features_ngenes" = "features_ngenes", "scvi_params"= "scvi_params"))
compare_eval2 = compare_eval[compare_eval$method_full=="scVI" & !is.na(compare_eval$reduction_full) & !is.na(compare_eval$reduction_neuron),]

combine_eval = dplyr::bind_rows(all_evaluation_results_full %>% dplyr::mutate(source="full") , all_evaluation_results_neurons[!grepl("fromFull",all_evaluation_results_neurons$reduction),] %>% dplyr::mutate(source="neuron"))
colnames(combine_eval)

### show all together:
require(ggplot2)
ggplot2::ggplot(combine_eval %>% dplyr::filter(method == "scVI") ,aes(mixingknn,mean_knn_purity,size = asw,color=source))+geom_point()#+geom_abline(slope = 1)


ggplot2::ggplot(combine_eval %>% dplyr::filter(method == "scVI") ,aes(mixingknn, mean_knn_purity ,size = mean_asw_dataset_celltype,color=source))+geom_point()#+geom_abline(slope = 1)


## plot features_ngenes
plot_df = combine_eval[combine_eval$method=="scVI" & combine_eval$features_ngenes != 1250,]
plot_df$features_ngenes= factor(plot_df$features_ngenes,levels = as.character(sort(as.numeric(unique(plot_df$features_ngenes)))))
p_asw_nfeatures = ggplot(plot_df, aes(x=features_ngenes, y=asw)) +
  geom_boxplot(aes(fill=source)) +
  geom_point(position=position_dodge(width=0.75),aes(group=source))+theme(text=element_text(size = 20))
p_asw_nfeatures

p_knnpurity_nfeatures = ggplot(plot_df, aes(x=features_ngenes, y=mean_knn_purity)) +
  geom_boxplot(aes(fill=source)) +
  geom_point(position=position_dodge(width=0.75),aes(group=source))+theme(text=element_text(size = 20))
p_knnpurity_nfeatures

p_mixingknn_nfeatures = ggplot(plot_df, aes(x=features_ngenes, y=mixingknn)) +
  geom_boxplot(aes(fill=source)) +
  geom_point(position=position_dodge(width=0.75),aes(group=source))
p_mixingknn_nfeatures

##########
###  run basic comaprison neurons with full evaluation subset to neurons !
##########

#### with fromFull
all_evaluation_results_neurons$source = "neuron"
all_evaluation_results_neurons$source[grepl("fromFull",all_evaluation_results_neurons$reduction)] = "full"

ggplot2::ggplot(all_evaluation_results_neurons %>% dplyr::filter(method == "scVI") ,aes(mixingknn, mean_knn_purity ,size = asw,color=source))+
  geom_point()+theme(text=element_text(size = 20))

#ggplot2::ggplot(all_evaluation_results_neurons %>% dplyr::filter(method == "scVI"), aes(id,y=mixingknn,fill=source,group=source)) + geom_col( position = "dodge") + facet_wrap(~ features_ngenes)

# subset
onlyShared_versions = all_evaluation_results_neurons %>% dplyr::group_by(ndim,scvi_params,features_ngenes) %>% dplyr::add_count(name = "occurs_multiple") %>% dplyr::filter(occurs_multiple > 1)
onlyShared_versions$id = paste0(onlyShared_versions$scvi_params,"_",onlyShared_versions$ndim,"_",onlyShared_versions$features_ngenes)
# plot params
ggplot2::ggplot(onlyShared_versions, aes(id,y=mixingknn,fill=source,group=source)) + geom_col( position = "dodge") + theme(axis.text.x = element_text(angle = 45, hjust=1,size=15),text = element_text(size=15))
ggplot2::ggplot(onlyShared_versions, aes(id,y=mean_knn_purity,fill=source,group=source)) + geom_col( position = "dodge") + theme(axis.text.x = element_text(angle = 45, hjust=1,size=15),text = element_text(size=15))
#ggplot2::ggplot(onlyShared_versions, aes(id,y=mean_asw_dataset_celltype,fill=source,group=source)) + geom_col( position = "dodge") + theme(axis.text.x = element_text(angle = 45, hjust=1,size=15),text = element_text(size=15))
ggplot2::ggplot(onlyShared_versions, aes(id,y=asw,fill=source,group=source)) + geom_col( position = "dodge") + theme(axis.text.x = element_text(angle = 45, hjust=1,size=15),text = element_text(size=15))











