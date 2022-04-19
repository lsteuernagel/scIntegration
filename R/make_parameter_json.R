## This script creates the json with general parameters --> make other jsons or edit this script if other params should be used
# requires some manually decisions which are added here.

param_list = list()

# must be loaded from params:
param_list$integration_folder_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/"#"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_integration/"
#param_list$processed_suffix = "_seurat_processed"

# for final merged object:
param_list$merged_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/hypoMap_merged_filtered.rds"#"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/hypoMap_merged_filtered.rds"
param_list$new_name_suffix = "hypoMap_integrated"#"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/hypoMap_merged_filtered.rds"

# signature for evaluation
param_list$celltype_signature_file = "data/hypothalamus_celltype_signatures.json"
param_list$auc_backup_file = "data/hypoMap_celltype_auc_per_cell_result.txt"
param_list$genes_to_exclude_file = "data/features_exclude_list.json"

# general
param_list$n_cores = 50
param_list$id_column = "Cell_ID"
param_list$global_seed = 123456
param_list$sample_column = "Sample_ID"
param_list$batch_var = "Batch_ID"
param_list$feature_set_sizes = c(750,1000,1250,1500,2000,2500,3000)
param_list$assay_name = "RNA"

# aucell
param_list$auc_max_rank=700
param_list$block_size=10000
param_list$alpha=1.5
param_list$thrP=0.01
param_list$smallestPopPercent=0.01
param_list$auc_max_pos_thresh=0.05
param_list$auc_min_pos_thresh=0.15
param_list$detected_cells_filename="detected_celltypes.json"

# downsample
param_list$id_file_name = "downsampled_ids_for_evaluation.json"
param_list$target_sub_sample = 38700
param_list$stepsize = 200

# integration
param_list$scvi_models_to_test = 30 # depends on the number of the models in scvi arg list
param_list$latent_space_sizes = c(50,65,80,95,110,140)

# evaluation:
param_list$ntrees_mixing = 20000
param_list$sampsize_pct = 0.3333
param_list$max_for_norm = 0.01
param_list$k_param = 20
param_list$dist_type = "cosine"
## asw
param_list$subset_cells = TRUE
param_list$target_clusterN = 200
param_list$start_res_asw = 5
param_list$end_res_asw = 15

# save
scUtils::writeList_to_JSON(param_list,filename = "data/parameters_integration_v2_3.json")
