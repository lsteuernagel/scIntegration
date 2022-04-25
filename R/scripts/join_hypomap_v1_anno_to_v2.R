

# load hypomap v1
#hypoMap_seurat_neurons_v1 = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_largeFiles/hypothalamus_neurons_map.rds")

# load hypomap v2
#hypoMap_seurat = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/hypoMap_merged_filtered.rds")


## associated v1 Cell_IDs to v2
hypoMap_seurat_neurons_v1_metadata = hypoMap_seurat_neurons_v1@meta.data
hypoMap_seurat_v2_metadata = hypoMap_seurat@meta.data
hypoMap_seurat_neurons_v1_metadata$Barcode = stringr::str_extract(hypoMap_seurat_neurons_v1_metadata$Cell_ID,pattern = "[ACGTN]{5,}") # \\-?[0-9]?
hypoMap_seurat_v2_metadata$Barcode = stringr::str_extract(hypoMap_seurat_v2_metadata$Cell_ID,pattern = "[ACGTN]{5,}") # \\-?[0-9]?
hypoMap_seurat_v2_metadata$Barcode[hypoMap_seurat_v2_metadata$Dataset=="Mousebrainorg10x"] = substr(hypoMap_seurat_v2_metadata$Barcode[hypoMap_seurat_v2_metadata$Dataset=="Mousebrainorg10x"],3,14)
hypoMap_seurat_v2_metadata$Dataset_old = gsub("10x|Dropseq","",hypoMap_seurat_v2_metadata$Dataset)
hypoMap_seurat_v2_metadata$Dataset_old[hypoMap_seurat_v2_metadata$Dataset_old == "KimDev"] = "kimDev"
hypoMap_seurat_v2_metadata$Dataset_old[hypoMap_seurat_v2_metadata$Dataset_old == "Lee"] = "Lee_Idol"
hypoMap_seurat_v2_metadata$Dataset_old[hypoMap_seurat_v2_metadata$Dataset_old == "wen"] = "wenDropSeq"
hypoMap_seurat_v2_metadata$Dataset_old[hypoMap_seurat_v2_metadata$Dataset_old == "Wen"] = "wen10x"
# make barcode+dataset combination
hypoMap_seurat_neurons_v1_metadata$pasted_barcode = paste0(hypoMap_seurat_neurons_v1_metadata$Dataset,"_",hypoMap_seurat_neurons_v1_metadata$Barcode)
hypoMap_seurat_v2_metadata$pasted_barcode = paste0(hypoMap_seurat_v2_metadata$Dataset_old,"_",hypoMap_seurat_v2_metadata$Barcode)
# join
joined_metadata = dplyr::left_join(hypoMap_seurat_neurons_v1_metadata %>% dplyr::select(pasted_barcode,Cell_ID,K169_pruned,Dataset) ,
                                   hypoMap_seurat_v2_metadata %>% dplyr::select(pasted_barcode,Cell_ID) ,by=c("pasted_barcode"="pasted_barcode"),suffix = c("_v1", "_v2"))
# there will be duplicates: remove
duplicated_barcodes = joined_metadata$pasted_barcode[duplicated(joined_metadata$pasted_barcode)]
joined_metadata = joined_metadata %>% dplyr::filter(! pasted_barcode %in% duplicated_barcodes) %>% dplyr::filter(!is.na(Cell_ID_v2))

table(joined_metadata$Dataset)

## make dataframe with v2 Cell IDs and v1 cluster ids
data_to_export = joined_metadata %>% dplyr::select(Cell_ID = Cell_ID_v2,cell_labels = K169_pruned) %>% as.data.frame()
table(data_to_export$cell_labels)
data.table::fwrite(data_to_export,"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/cell_labels_K169_v1.csv")


