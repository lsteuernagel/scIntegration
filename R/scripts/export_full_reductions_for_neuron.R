list_of_reductions = unique(c(compare_eval2$reduction_full,all_evaluation_results_full$reduction[all_evaluation_results_full$mixingknn > 0.43 & all_evaluation_results_full$mean_asw_dataset_celltype > 0.06 & all_evaluation_results_full$mean_knn_purity > 0.3]))

# path to reductions:
integration_path =  "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/integration/scvi/"
integration_path_output =  "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_neurons_integration/integration/scvi_full/"
system(paste0("mkdir -p ",integration_path_output))

# source functions:
source("R/evaluation_functions.R")

# read metadata from full map:
full_map_meta = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/hypoMap_integrated_metadata.txt",data.table = F)
rownames(full_map_meta) = full_map_meta$Cell_ID

# go over all reductions:
for(i in 1:length(list_of_reductions)){
  message(i," of ",length(list_of_reductions))
  # get current name
  current_reduction = list_of_reductions[i]
  integration_file_name = paste0(current_reduction,".txt") #"scVI_0_300_0.1_3_256_gene_zinb_cov2..scVI..65..RNA.log.vst.split_Batch_ID.features.2000_faf95328efd432d8d2539ad16c0852db.txt"

  # define full  name and read
  full_name = paste0(integration_path,integration_file_name)
  embedding = read_embedding(filename_withpath = full_name,seurat_object_metadata = full_map_meta)

  # subset
  embedding = embedding[rownames(hypoMap_seurat_neurons@meta.data),]

  # save:
  current_matrix = cbind(Cell_ID=rownames(hypoMap_seurat_neurons@meta.data),embedding)
  filename = paste0(integration_path_output,current_reduction,"_fromFull.txt")
  data.table::fwrite(data.table::as.data.table(current_matrix),file=filename,sep="\t",col.names = TRUE)

}

