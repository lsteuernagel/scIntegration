
## read full file:
hypoMap_merged_filtered = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/hypoMap_merged_filtered.rds")

# downsample to 10000
hypoMap_merged_filtered_downsample = scUtils::downsample_balanced_iterative(hypoMap_merged_filtered@meta.data,target_sample_size = 10000,predictor_var = "Dataset",stepsize = 100,global_seed = 123456)

# subset
hypoMap_merged_filtered_downsampled = subset(hypoMap_merged_filtered,cells = hypoMap_merged_filtered_downsample$Cell_ID)

# save in test folder
test_folder = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_integration_test/"
system(paste0("mkdir -p ",paste0(test_folder)))
saveRDS(hypoMap_merged_filtered_downsampled,paste0(test_folder,"hypoMap_merged_filtered.rds"))
