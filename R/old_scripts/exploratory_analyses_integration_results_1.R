
integration_file_name = "scVI_0_200_0.01_4_256_gene_zinb_cov2..scVI..50..RNA.log.vst.split_Batch_ID.features.1500_e06093f25a8fc603ec7f9eb5186b44f8.txt"

integration_file_name = "scVI_0_400_0.05_2_256_gene_zinb_cov2..scVI..110..RNA.log.vst.split_Batch_ID.features.2000_2d15d1b4808a8f4d44f3e609493a7c3a.txt"

integration_file_name = "scVI_0_300_0.05_3_256_gene_zinb_cov2..scVI..110..RNA.log.vst.split_Batch_ID.features.1250_fb23407b78b330da59bc5a1e8f539053.txt"

integration_file_name = "scVI_0_400_0.1_4_256_gene_zinb_cov2..scVI..65..RNA.log.vst.split_Batch_ID.features.3000_e9b11dacd25decfcd89c1623473e5352.txt"

integration_file_name = "scVI_0_300_0.01_3_256_gene_zinb_cov2..scVI..85..RNA.log.vst.split_Batch_ID.features.1500_99244d8dce5a5588784f504eab99f0eb.txt"


new_name = gsub(".txt","",integration_file_name )#"scvi"

# load seurat
#hypoMap_seurat = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/hypoMap_merged_filtered.rds")

# get inetgration
integration_path =  "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/integration/scvi/"
full_name = paste0(integration_path,integration_file_name)

source("R/evaluation_functions.R")
embedding = read_embedding(filename_withpath = full_name,seurat_object = hypoMap_seurat)

#add
new_name = gsub(".txt","",integration_file_name )#"scvi"
dimred <- Seurat::CreateDimReducObject(
  embeddings = as.matrix(embedding),
  stdev = as.numeric(apply(embedding, 2, stats::sd)),
  assay = "RNA",
  key = new_name
)
# add
hypoMap_seurat@reductions[[new_name]] = dimred

# Run UMAP
hypoMap_seurat = Seurat::RunUMAP(hypoMap_seurat,reduction = new_name,dims=1:ncol(embedding), seed.use = 123456,reduction.name =paste0("umap_",new_name),reduction.key =paste0("umap_",new_name))

# plot
p1= Seurat::DimPlot(hypoMap_seurat,group.by =  "Dataset",raster = F,order = F,shuffle = TRUE,reduction=paste0("umap_",new_name))#reduction=paste0("umap_",new_name)
p1r = scUtils::rasterize_ggplot(p1,pixel_raster = 2048)
p1r


# plot
p2= Seurat::FeaturePlot(hypoMap_seurat,features = "Fezf1",raster = F,order = TRUE,reduction=paste0("umap_",new_name))#reduction=paste0("umap_",new_name)
p2r = scUtils::rasterize_ggplot(p2,pixel_raster = 2048)
p2r

#### evelaution cell types

celltype_id_list = jsonlite::read_json("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/detected_celltypes.json")
celltype_id_list = lapply(celltype_id_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

library(purrr)
celltype_df <- purrr::map_df(celltype_id_list, ~as.data.frame(.x), .id="celltype")
colnames(celltype_df) = c("detected_celltype","Cell_ID")

temp_meta = dplyr::left_join(hypoMap_seurat@meta.data,celltype_df,by=c("Cell_ID")) %>% as.data.frame()
rownames(temp_meta) = temp_meta$Cell_ID
#hypoMap_seurat@meta.data = temp_meta

# plot
p3= Seurat::DimPlot(hypoMap_seurat,group.by =  "detected_celltype",raster = F,order = TRUE,shuffle = TRUE,reduction=paste0("umap_",new_name))#reduction=paste0("umap_",new_name)
p3r = scUtils::rasterize_ggplot(p3,pixel_raster = 2048)
p3r

####check on agrp neurons in high asw result

# use cell selector to get relevant cells by hand
selected_cells = Seurat::CellSelector(plot = p3)

# make a subset object
agrp_subset_high_asw = subset(hypoMap_seurat,cells = selected_cells)

## re run
DimPlot(object = agrp_subset_high_asw,group.by = "Dataset",reduction=paste0("umap_",new_name))

## re run
agrp_subset_high_asw_processed = RunUMAP(agrp_subset_high_asw,reduction = new_name,dims = 1:ncol(agrp_subset_high_asw_processed@reductions[[new_name]]@cell.embeddings), reduction.name = paste0("umap_","small"),reduction.key = paste0("umap_","small"))
agrp_subset_high_asw_processed = FindNeighbors(agrp_subset_high_asw_processed,reduction = new_name,dims = 1:ncol(agrp_subset_high_asw_processed@reductions[[new_name]]@cell.embeddings),graph.name = "newsnn")
agrp_subset_high_asw_processed = FindClusters(agrp_subset_high_asw_processed,resolution = 1,random.seed = 123456,graph.name =  "newsnn")

# plot
DimPlot(object = agrp_subset_high_asw_processed,group.by = "seurat_clusters", reduction = paste0("umap_","small"),label=TRUE,label.size = 5)
DimPlot(object = agrp_subset_high_asw_processed,group.by = "Dataset", reduction = paste0("umap_","small"))

# markers
Idents(agrp_subset_high_asw_processed) = "seurat_clusters"
#agrp_subset_all_markers = FindAllMarkers(agrp_subset_high_asw_processed,logfc.threshold = 0.3,min.pct = 0.1,min.diff.pct = 0.05,only.pos = TRUE)
agrp_subset_all_markers$specificity = (agrp_subset_all_markers$pct.1 / agrp_subset_all_markers$pct.2) * agrp_subset_all_markers$avg_log2FC

# FetaurePlo
FeaturePlot( agrp_subset_high_asw_processed,features = "Sst", reduction = paste0("umap_","small"))
FeaturePlot( agrp_subset_high_asw_processed,features = "Gpr101", reduction = paste0("umap_","small"))
FeaturePlot( agrp_subset_high_asw_processed,features = "Kcnmb2", reduction = paste0("umap_","small"))

###
feature_list = jsonlite::read_json("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/features/feature_sets.json")
feature_list = lapply(feature_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

a1=data.frame(feature_list$RNA.log.vst.split_Batch_ID.features.3000)






