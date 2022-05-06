
##########
### Define input files
##########

params_integration = jsonlite::read_json("data/parameters_integration_v2_3.json")
evaluation_folder = paste0(params_integration$integration_folder_path,"evaluation/")
evaluation_mixingrf_file = paste0(evaluation_folder,"all_mixing_rf.txt")
evaluation_purityknn_file = paste0(evaluation_folder,"all_purity_knn.txt")
evaluation_mixingknn_file = paste0(evaluation_folder,"all_mixing_knn.txt")
evaluation_purityasw_file = paste0(evaluation_folder,"all_purity_asw.txt")
require(magrittr)

#data.table::fwrite(,file = evaluation_mixingrf_file)

##########
### Read eval results
##########

## read files:
evaluation_mixingrf = data.table::fread(evaluation_mixingrf_file,data.table = F) %>% dplyr::rename( mixingrf = value)
evaluation_mixingknn = data.table::fread(evaluation_mixingknn_file,data.table = F) %>% dplyr::rename( mixingknn = value)
evaluation_purityknn = data.table::fread(evaluation_purityknn_file,data.table = F) %>% dplyr::rename( purityknn = value)
evaluation_purityasw = data.table::fread(evaluation_purityasw_file,data.table = F) #%>% dplyr::rename( asw = value)
colnames(evaluation_purityasw) =c(colnames(evaluation_purityasw)[2:length(colnames(evaluation_purityasw))],"placeholder")
evaluation_purityasw$reduction = evaluation_purityasw$placeholder
## additional:
# all_knownLabel_K169_asw = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/evaluation/all_knownLabel_K169_asw.txt",data.table = F)
# colnames(all_knownLabel_K169_asw) =c("reduction","knownasw_mean_euc","knownasw_mean_cos")
# grouped_knownLabel_K169_asw = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/evaluation/grouped_knownLabel_K169_asw.txt",data.table = F)

##########
### Join for overview
##########

## join together for simple overview
evaluation_purityknn_mean = evaluation_purityknn %>% dplyr::group_by(reduction) %>% dplyr::summarise(mean_knn_purity = mean(purityknn))
all_evaluation_results = dplyr::full_join(evaluation_mixingrf,evaluation_mixingknn,by="reduction")
all_evaluation_results = dplyr::full_join(all_evaluation_results,evaluation_purityknn_mean,by="reduction")
all_evaluation_results = dplyr::full_join(all_evaluation_results,evaluation_purityasw[,c("reduction","silhouette_score_euclidean")],by="reduction")
#all_evaluation_results = dplyr::full_join(all_evaluation_results,all_knownLabel_K169_asw[,c("reduction","knownasw_mean_euc")],by="reduction")
all_evaluation_results = all_evaluation_results %>% dplyr::select(reduction, mixingrf, mixingknn , mean_knn_purity,asw = silhouette_score_euclidean,label_asw = knownasw_mean_euc)

require(ggplot2)
ggplot2::ggplot(all_evaluation_results,aes(mixingrf,mixingknn))+geom_point()#+geom_abline(slope=1)

##########
### Normalize and add information from reduction name
##########

## norm scores
min_max = function(x,minVal=NULL,maxVal=NULL){
  if(is.null(minVal)){minVal =min(x,na.rm = TRUE) }
  if(is.null(maxVal)){maxVal =max(x,na.rm = TRUE) }
  return((x-minVal) / (maxVal-minVal))
}
all_evaluation_results$mean_purity_knn_norm = round(min_max(all_evaluation_results$mean_knn_purity,minVal = 0.3,maxVal = NULL),4)*100
all_evaluation_results$mixingknn_norm = round(min_max(all_evaluation_results$mixingknn,minVal = 0,maxVal = NULL),4)*100
all_evaluation_results$mixingrf_norm = round(min_max(all_evaluation_results$mixingrf,minVal = 0.8,maxVal = NULL),4)*100
all_evaluation_results$asw_norm = round(min_max(all_evaluation_results$asw,minVal = NULL,maxVal = NULL),4)*100
all_evaluation_results$mixing_score = all_evaluation_results$mixingrf_norm*0.5+all_evaluation_results$mixingknn_norm*0.5
all_evaluation_results$purity_score =  all_evaluation_results$mean_purity_knn_norm*0.5 + all_evaluation_results$asw_norm*0.5

# add information
# all_evaluation_results$assay=all_evaluation_results$reduction %>% stringr::str_extract(pattern="\\.\\..[a-zA-Z]+\\.\\.") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "")
# all_evaluation_results$ndim=all_evaluation_results$reduction %>% stringr::str_extract(pattern="\\.\\..[0-9]+\\.\\.") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.\\.",replacement = "")
# all_evaluation_results$features=all_evaluation_results$reduction %>% stringr::str_extract(pattern="\\.\\..*\\.features") %>% stringr::str_replace(pattern = "\\.\\..*\\.\\.",replacement = "") %>% stringr::str_replace(pattern = "\\.features",replacement = "")
# all_evaluation_results$method = stringr::str_extract(all_evaluation_results$reduction,pattern="[a-zA-Z]+_") %>% stringr::str_replace(pattern = "_",replacement = "")
# all_evaluation_results$features_ngenes = stringr::str_extract(all_evaluation_results$reduction,pattern="features\\.[0-9]+") %>% stringr::str_replace(pattern = "features\\.",replacement = "") %>% as.numeric()
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


##########
### Visualize
##########

colnames(all_evaluation_results)

require(ggplot2)
ggplot2::ggplot(all_evaluation_results,aes(mixing_score,purity_score,color=method))+geom_point()#+geom_abline(slope=1)

ggplot2::ggplot(all_evaluation_results[all_evaluation_results$method=="scVI",],aes(mean_knn_purity,asw,color=as.numeric(mixing_score)))+geom_point()#+geom_abline(slope=1)
ggplot2::ggplot(all_evaluation_results[all_evaluation_results$method=="scVI",],aes(mixingknn,asw,color=as.numeric(mean_knn_purity)))+geom_point()#+geom_abline(slope=1)


plot_df = all_evaluation_results[all_evaluation_results$method=="scVI" & all_evaluation_results$ndim %in% c("110","65","85"),]
plot_df$ndim= factor(plot_df$ndim,levels = c("65","85","110"))
plot_df$max_epochs= factor(plot_df$max_epochs,levels = c("50","100","150","200","300","400"))
plot_df$dropout_rate= factor(plot_df$dropout_rate,levels = c("0.01","0.05","0.1"))
plot_df$features_ngenes= factor(plot_df$features_ngenes,levels = as.character(sort(as.numeric(unique(plot_df$features_ngenes)))))

#### asw

## plot ndim
p_asw_ndim = ggplot(plot_df, aes(x=ndim, y=asw)) +
  geom_boxplot(aes(fill=nlayers)) +
  geom_point(position=position_dodge(width=0.75),aes(group=nlayers))

## plot epochs
p_asw_epochs = ggplot(plot_df, aes(x=max_epochs, y=asw)) +
  geom_boxplot(aes(fill=nlayers)) +
  geom_point(position=position_dodge(width=0.75),aes(group=nlayers))

## plot dropout_rate
p_asw_dropout = ggplot(plot_df, aes(x=dropout_rate, y=asw)) +
  geom_boxplot(aes(fill=nlayers)) +
  geom_point(position=position_dodge(width=0.75),aes(group=nlayers))

## plot features_ngenes
p_asw_nfeatures = ggplot(plot_df, aes(x=features_ngenes, y=asw)) +
  geom_boxplot(aes(fill=nlayers)) +
  geom_point(position=position_dodge(width=0.75),aes(group=nlayers))

cowplot::plot_grid(p_asw_ndim,p_asw_epochs,p_asw_dropout,p_asw_nfeatures)+ggtitle("ASW")

#### mean_knn_purity

## plot ndim
p_mean_knn_purity_ndim = ggplot(plot_df, aes(x=ndim, y=mean_knn_purity)) +
  geom_boxplot(aes(fill=nlayers)) +
  geom_point(position=position_dodge(width=0.75),aes(group=nlayers))

## plot epochs
p_mean_knn_purity_epochs = ggplot(plot_df, aes(x=max_epochs, y=mean_knn_purity)) +
  geom_boxplot(aes(fill=nlayers)) +
  geom_point(position=position_dodge(width=0.75),aes(group=nlayers))

## plot dropout_rate
p_mean_knn_purity_dropout = ggplot(plot_df, aes(x=dropout_rate, y=mean_knn_purity)) +
  geom_boxplot(aes(fill=nlayers)) +
  geom_point(position=position_dodge(width=0.75),aes(group=nlayers))

## plot features_ngenes
p_mean_knn_purity_nfeatures = ggplot(plot_df, aes(x=features_ngenes, y=mean_knn_purity)) +
  geom_boxplot(aes(fill=nlayers)) +
  geom_point(position=position_dodge(width=0.75),aes(group=nlayers))

#### mixingknn

## plot ndim
p_mixingknn_ndim = ggplot(plot_df, aes(x=ndim, y=mixingknn)) +
  geom_boxplot(aes(fill=nlayers)) +
  geom_point(position=position_dodge(width=0.75),aes(group=nlayers))

## plot epochs
p_mixingknn_epochs = ggplot(plot_df, aes(x=max_epochs, y=mixingknn)) +
  geom_boxplot(aes(fill=nlayers)) +
  geom_point(position=position_dodge(width=0.75),aes(group=nlayers))

## plot dropout_rate
p_mixingknn_dropout = ggplot(plot_df, aes(x=dropout_rate, y=mixingknn)) +
  geom_boxplot(aes(fill=nlayers)) +
  geom_point(position=position_dodge(width=0.75),aes(group=nlayers))

## plot features_ngenes
p_mixingknn_nfeatures = ggplot(plot_df, aes(x=features_ngenes, y=mixingknn)) +
  geom_boxplot(aes(fill=nlayers)) +
  geom_point(position=position_dodge(width=0.75),aes(group=nlayers))

# ## plot features_ngenes
# p_mixingknn_nlayers = ggplot(plot_df, aes(x=nlayers, y=mixingknn)) +
#   geom_boxplot(aes(fill=nlayers)) +
#   geom_point(position=position_dodge(width=0.75),aes(group=nlayers))

#### mixingrf

## plot ndim
p_mixingrf_ndim = ggplot(plot_df, aes(x=ndim, y=mixingrf)) +
  geom_boxplot(aes(fill=nlayers)) +
  geom_point(position=position_dodge(width=0.75),aes(group=nlayers))

## plot epochs
p_mixingrf_epochs = ggplot(plot_df, aes(x=max_epochs, y=mixingrf)) +
  geom_boxplot(aes(fill=nlayers)) +
  geom_point(position=position_dodge(width=0.75),aes(group=nlayers))

## plot dropout_rate
p_mixingrf_dropout = ggplot(plot_df, aes(x=dropout_rate, y=mixingrf)) +
  geom_boxplot(aes(fill=nlayers)) +
  geom_point(position=position_dodge(width=0.75),aes(group=nlayers))

## plot features_ngenes
p_mixingrf_nfeatures = ggplot(plot_df, aes(x=features_ngenes, y=mixingrf)) +
  geom_boxplot(aes(fill=nlayers)) +
  geom_point(position=position_dodge(width=0.75),aes(group=nlayers))


#### arrange plots

cowplot::plot_grid(p_mixingknn_nfeatures,p_mixingrf_nfeatures,p_mean_knn_purity_nfeatures,p_asw_nfeatures)
cowplot::plot_grid(p_mixingknn_dropout,p_mixingrf_dropout,p_mean_knn_purity_dropout,p_asw_dropout)
cowplot::plot_grid(p_mixingknn_epochs,p_mixingrf_epochs,p_mean_knn_purity_epochs,p_asw_epochs)
cowplot::plot_grid(p_mixingknn_ndim,p_mixingrf_ndim,p_mean_knn_purity_ndim,p_asw_ndim)

##########
### per cell type ASW:
##########

### load per cell type ASW:
all_eval_files = list.files(evaluation_folder,pattern="all")
all_eval_files = all_eval_files[grepl("_AuthorCellTypes",all_eval_files)]
table_list = list()
for(i in 1:length(all_eval_files)){
  label = gsub("AuthorCellTypes.csv_AuthorCellType_","",all_eval_files[[i]])
  label = gsub(".txt","",label)
  table_list[[label]] = data.table::fread(paste0(evaluation_folder,all_eval_files[[i]]))
  table_list[[label]]$Dataset = label
}
asw_all_author_celltype = do.call(rbind,table_list)
colnames(asw_all_author_celltype) = c("reduction","asw_eulid","asw_cosine","dataset")

asw_all_author_celltype_wide = asw_all_author_celltype %>% dplyr::select(-asw_cosine) %>% tidyr::spread(key="dataset",value="asw_eulid")
asw_all_author_celltype_wide_norm=as.data.frame(cbind(reduction=asw_all_author_celltype_wide$reduction,apply(asw_all_author_celltype_wide[,2:ncol(asw_all_author_celltype_wide)],2,function(x){ ( x-min(x) ) / (max(x) - min(x)) })))
asw_all_author_celltype_wide_norm[,2:ncol(asw_all_author_celltype_wide_norm)] = apply(asw_all_author_celltype_wide_norm[,2:ncol(asw_all_author_celltype_wide_norm)],2,as.numeric)
all_evaluation_results_extended = dplyr::left_join(all_evaluation_results,asw_all_author_celltype_wide_norm,by="reduction")

### load per cell type ASW grouped:

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

asw_grouped_author_celltype_summarized = asw_grouped_author_celltype %>% dplyr::group_by(reduction) %>%
  dplyr::summarise(mean_across_celltypes = mean(asw_eulid))

all_evaluation_results_extended = dplyr::left_join(all_evaluation_results,asw_grouped_author_celltype_summarized,by="reduction")

# asw_grouped_author_celltype_wide = asw_grouped_author_celltype %>% dplyr::select(-asw_cosine) %>% tidyr::spread(key="reduction",value="asw_eulid")
# colnames(asw_grouped_author_celltype_wide)
# ggplot(asw_grouped_author_celltype_wide,aes(x=scvi_hypoMap_v1,y=scvi_hypoMap_neurons_v1,color=dataset))+geom_point()+geom_abline(slope = 1)
#

ggplot(all_evaluation_results_extended,aes(asw,mean_across_celltypes,color=method))+geom_point()

##############

comapre_purity_knn_directly = evaluation_purityknn %>% dplyr::filter(reduction %in% c("scVI_0_400_0.1_4_256_gene_zinb_cov2..scVI..65..RNA.log.vst.split_Batch_ID.features.3000_e9b11dacd25decfcd89c1623473e5352",
                                                        "scVI_0_400_0.05_2_256_gene_zinb_cov2..scVI..110..RNA.log.vst.split_Batch_ID.features.2000_2d15d1b4808a8f4d44f3e609493a7c3a"))



comapre_purity_knn_directly_wide = comapre_purity_knn_directly %>% tidyr::spread(key="reduction",value="purityknn")
comapre_purity_knn_directly_wide$diff = comapre_purity_knn_directly_wide[,2] - comapre_purity_knn_directly_wide[,3]
comapre_purity_knn_directly_wide$label_highdiff =NA
  comapre_purity_knn_directly_wide$label_highdiff[abs(comapre_purity_knn_directly_wide$diff) > 0.03] = comapre_purity_knn_directly_wide$celltype[abs(comapre_purity_knn_directly_wide$diff) > 0.03]

ggplot(comapre_purity_knn_directly_wide %>% dplyr::rename(high_asw=scVI_0_400_0.1_4_256_gene_zinb_cov2..scVI..65..RNA.log.vst.split_Batch_ID.features.3000_e9b11dacd25decfcd89c1623473e5352,
                                                          high_purity = scVI_0_400_0.05_2_256_gene_zinb_cov2..scVI..110..RNA.log.vst.split_Batch_ID.features.2000_2d15d1b4808a8f4d44f3e609493a7c3a),
       aes(high_asw,high_purity))+geom_point()+geom_abline(slope=1)+geom_text(aes(label=label_highdiff))

###### key todos

#ggplot(comapre_purity_knn_directly,

#### make qrfp subset

# does the subcluster make sense ?

#### highlight gpr50 / pgl

#####

comapre_asw_grouped_author_celltype_directly = asw_grouped_author_celltype %>% dplyr::select(asw_eulid,reduction,celltype) %>% dplyr::filter(reduction %in% c("scVI_0_400_0.1_4_256_gene_zinb_cov2..scVI..65..RNA.log.vst.split_Batch_ID.features.3000_e9b11dacd25decfcd89c1623473e5352",
                                                                                      "scVI_0_400_0.05_2_256_gene_zinb_cov2..scVI..110..RNA.log.vst.split_Batch_ID.features.2000_2d15d1b4808a8f4d44f3e609493a7c3a"))
comapre_asw_grouped_author_celltype_directly = comapre_asw_grouped_author_celltype_directly %>% dplyr::distinct(reduction,celltype,.keep_all = TRUE)
comapre_asw_grouped_author_celltype_directly_wide = comapre_asw_grouped_author_celltype_directly %>% tidyr::spread(key="reduction",value="asw_eulid")
comapre_asw_grouped_author_celltype_directly_wide$diff = comapre_asw_grouped_author_celltype_directly_wide[,2] - comapre_asw_grouped_author_celltype_directly_wide[,3]
comapre_asw_grouped_author_celltype_directly_wide$label_highdiff =NA
comapre_asw_grouped_author_celltype_directly_wide$label_highdiff[abs(comapre_asw_grouped_author_celltype_directly_wide$diff) > 0.15] =
  comapre_asw_grouped_author_celltype_directly_wide$celltype[abs(comapre_asw_grouped_author_celltype_directly_wide$diff) > 0.15]

ggplot(comapre_asw_grouped_author_celltype_directly_wide %>% dplyr::rename(high_asw=scVI_0_400_0.1_4_256_gene_zinb_cov2..scVI..65..RNA.log.vst.split_Batch_ID.features.3000_e9b11dacd25decfcd89c1623473e5352,
                                                          high_purity = scVI_0_400_0.05_2_256_gene_zinb_cov2..scVI..110..RNA.log.vst.split_Batch_ID.features.2000_2d15d1b4808a8f4d44f3e609493a7c3a),
       aes(high_asw,high_purity))+geom_point()+geom_abline(slope=1)+geom_text(aes(label=label_highdiff))







