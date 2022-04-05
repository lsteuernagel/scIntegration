## This script creates the json with general parameters --> make other jsons or edit this script if other params should be used
# requires some manually decisions which are added here.

param_list = list()

# must be loaded from params:
param_list$integration_folder_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_integration_test/"#"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_integration/"
#param_list$processed_suffix = "_seurat_processed"

# for final merged object:
param_list$merged_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_integration_test/hypoMap_merged_filtered.rds"#"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/hypoMap_merged_filtered.rds"
param_list$feature_set_sizes = c(750,1000,1250,1500,2000,2500,3000)

# signature for evaluation
param_list$celltype_signature_path = "data/hypothalamus_celltype_signatures.json"
param_list$genes_to_exclude_file = "data/features_exclude_list.json"

# processing
param_list$n_cores = 50
param_list$id_column = "Cell_ID"
param_list$global_seed = 123456

# save
scUtils::writeList_to_JSON(param_list,filename = "data/parameters_integration_v2_1.json")
