
feature_exclude_list = unlist(jsonlite::read_json("data/features_exclude_list.json"))
# high pruity:
integration_file_name = "scVI_0_400_0.1_3_256_gene_zinb_cov2..scVI..85..RNA.log.vst.split_Batch_ID.features.1000_dee51dad39ec51d5b20cec93d2df348f.txt"
# this one has a very high asw but detailed analysis showed clear over correctioN:
integration_file_name = "scVI_0_400_0.15_4_256_gene_zinb_cov2..scVI..85..RNA.log.vst.split_Batch_ID.features.750_27466e1704ffc1d3f82787db0b6e7b8d.txt"

# this one has higher purity:
integration_file_name = "scVI_0_400_0.1_2_256_gene_zinb_cov2..scVI..85..RNA.log.vst.split_Batch_ID.features.1000_6a908a78b57b96b7d1c35e4300f8c194.txt"

# a good balanced neuron version:
integration_file_name = "scVI_0_300_0.1_3_256_gene_zinb_cov2..scVI..65..RNA.log.vst.split_Batch_ID.features.1000_562216427a8b35b2658273fecaf712a2.txt"

# the best neuron version when just looking at the metrics:
integration_file_name = "scVI_0_400_0.15_3_256_gene_zinb_cov2..scVI..85..RNA.log.vst.split_Batch_ID.features.750_64203fb2c726221727cd01767d952850.txt"


###load
new_name = gsub(".txt","",integration_file_name )#"scvi"
# load seurat
#hypoMap_seurat_neurons = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_neurons_integration/hypoMap_neurons.rds")

# get inetgration
integration_path =  "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_neurons_integration/integration/scvi/"
full_name = paste0(integration_path,integration_file_name)

source("R/evaluation_functions.R")
embedding = read_embedding(filename_withpath = full_name,seurat_object = hypoMap_seurat_neurons)
#add
new_name = gsub(".txt","",integration_file_name )#"scvi"
dimred <- Seurat::CreateDimReducObject(
  embeddings = as.matrix(embedding),
  stdev = as.numeric(apply(embedding, 2, stats::sd)),
  assay = "RNA",
  key = new_name
)
# add
hypoMap_seurat_neurons@reductions[[new_name]] = dimred

# Run UMAP
hypoMap_seurat_neurons = Seurat::RunUMAP(hypoMap_seurat_neurons,reduction = new_name,dims=1:ncol(embedding), seed.use = 123456,reduction.name =paste0("umap_",new_name),reduction.key =paste0("umap_",new_name))

# plot
p1= Seurat::DimPlot(hypoMap_seurat_neurons,group.by =  "Dataset",raster = F,order = F,shuffle = TRUE,reduction=paste0("umap_",new_name))#reduction=paste0("umap_",new_name)
p1r = scUtils::rasterize_ggplot(p1,pixel_raster = 2048)
p1r

# plot
p2= Seurat::FeaturePlot(hypoMap_seurat_neurons,features = "Grp",raster = F,order = TRUE,reduction=paste0("umap_",new_name))#reduction=paste0("umap_",new_name)
p2r = scUtils::rasterize_ggplot(p2,pixel_raster = 2048)
p2r

p1= Seurat::DimPlot(hypoMap_seurat_neurons,group.by =  "seurat_clusters",raster = F,order = F,shuffle = TRUE,reduction=paste0("umap_",new_name),label = TRUE)+NoLegend()#reduction=paste0("umap_",new_name)
p1r = scUtils::rasterize_ggplot(p1,pixel_raster = 2048)
p1r


#### evelaution cell types

celltype_id_list = jsonlite::read_json("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_neurons_integration/detected_celltypes.json")
celltype_id_list = lapply(celltype_id_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

library(purrr)
celltype_df <- purrr::map_df(celltype_id_list, ~as.data.frame(.x), .id="celltype")
colnames(celltype_df) = c("detected_celltype","Cell_ID")

temp_meta = dplyr::left_join(hypoMap_seurat_neurons@meta.data,celltype_df,by=c("Cell_ID")) %>% as.data.frame()
rownames(temp_meta) = temp_meta$Cell_ID
hypoMap_seurat_neurons@meta.data = temp_meta

# plot
p3= Seurat::DimPlot(hypoMap_seurat_neurons,group.by =  "detected_celltype.y",raster = F,order = TRUE,shuffle = TRUE,reduction=paste0("umap_",new_name))#reduction=paste0("umap_",new_name)
p3r = scUtils::rasterize_ggplot(p3,pixel_raster = 2048)
p3r


cellsh = hypoMap_seurat_neurons@meta.data$Cell_ID[hypoMap_seurat_neurons@meta.data$detected_celltype.y == "Vip_Nov"]
p3= Seurat::DimPlot(hypoMap_seurat_neurons,group.by =  "detected_celltype.y",raster = F,cells.highlight = cellsh,sizes.highlight = 0.1,reduction=paste0("umap_",new_name))#reduction=paste0("umap_",new_name)
p3r = scUtils::rasterize_ggplot(p3,pixel_raster = 2048)
p3r

##### load full map reduction and add:

integration_file_name = "scVI_0_300_0.1_3_256_gene_zinb_cov2..scVI..65..RNA.log.vst.split_Batch_ID.features.2000_faf95328efd432d8d2539ad16c0852db.txt"
source("R/evaluation_functions.R")
embedding = read_embedding(filename_withpath = full_name,seurat_object = hypoMap_seurat_neurons)


####
dimred_sub = hypoMap_seurat@reductions$scVI_1_400_0.1_2_256_gene_zinb_cov2..scVI..85..RNA.log.vst.split_Batch_ID.features.2000_5e31aa323c46e0fb15d6f94afc210cce
dimred_sub = subset(dimred_sub,cells = hypoMap_seurat_neurons@meta.data$Cell_ID)

# add
hypoMap_seurat_neurons@reductions[["from_full_map_good"]] = dimred_sub
# Run UMAP
hypoMap_seurat_neurons = Seurat::RunUMAP(hypoMap_seurat_neurons,reduction = "from_full_map_good",dims=1:ncol(dimred_sub), seed.use = 123456,
                                         reduction.name =paste0("umap_","from_full_map_good"),reduction.key =paste0("umap_","from_full_map_good"))
# plot
p1= Seurat::DimPlot(hypoMap_seurat_neurons,group.by =  "Dataset",raster = F,order = F,shuffle = TRUE,reduction=paste0("umap_","from_full_map_good"))#reduction=paste0("umap_",new_name)
p1r = scUtils::rasterize_ggplot(p1,pixel_raster = 2048)
p1r
# plot
p2= Seurat::FeaturePlot(hypoMap_seurat_neurons,features = "Tmem215",raster = F,order = TRUE,reduction=paste0("umap_","from_full_map_good"))#reduction=paste0("umap_",new_name)
p2r = scUtils::rasterize_ggplot(p2,pixel_raster = 2048)
p2r

p2= Seurat::FeaturePlot(hypoMap_seurat_neurons,features = c("Pnoc","Lepr"),blend = TRUE,blend.threshold = 0.1,raster = F,order = TRUE,reduction=paste0("umap_","from_full_map_good"))#reduction=paste0("umap_",new_name)


########## POMC

#selected_cells = Seurat::CellSelector(p2)
pomc_subset = subset(hypoMap_seurat_neurons,cells=selected_cells)
Seurat::DimPlot(pomc_subset,group.by =  "Dataset",raster = F,order = F,shuffle = TRUE,reduction=paste0("umap_","from_full_map"))
Seurat::FeaturePlot(pomc_subset,features = "Pomc",raster = F,order = TRUE,reduction=paste0("umap_","from_full_map"))

Seurat::DimPlot(pomc_subset,group.by =  "Dataset",raster = F,order = F,shuffle = TRUE,reduction=paste0("umap_",new_name))
Seurat::FeaturePlot(pomc_subset,features = "Pomc",raster = F,order = TRUE,reduction=paste0("umap_",new_name))

pomc_subset = FindNeighbors(pomc_subset,reduction =paste0("umap_","from_full_map"),dims = 1:ncol(pomc_subset@reductions[[paste0("umap_","from_full_map")]]@cell.embeddings))
pomc_subset = FindClusters(pomc_subset,resolution = 0.15)

Seurat::DimPlot(pomc_subset,group.by =  "seurat_clusters",raster = F,order = F,shuffle = TRUE,label=TRUE,reduction=paste0("umap_","from_full_map"))
pomc_sub_markers = FindMarkers(pomc_subset,ident.1 = c("0","3"),min.pct = 0.1,min.diff.pct = 0.05,only.pos = TRUE)
pomc_sub_markers$gene = rownames(pomc_sub_markers)
pomc_sub_markers$specificity = (pomc_sub_markers$pct.1 / pomc_sub_markers$pct.2) * pomc_sub_markers$avg_log2FC

#Seurat::FeaturePlot(pomc_subset,features = "Pcsk2",raster = F,order = TRUE,reduction=paste0("umap_",new_name))
Seurat::FeaturePlot(pomc_subset,features = "Marchf1",raster = F,order = TRUE,reduction=paste0("umap_","from_full_map"))


########## PMCH

#selected_cells = Seurat::CellSelector(p2)
pmch_subset = subset(hypoMap_seurat_neurons,cells=selected_cells)
Seurat::DimPlot(pmch_subset,group.by =  "Dataset",raster = F,order = F,shuffle = TRUE,reduction=paste0("umap_","from_full_map_good"))
Seurat::FeaturePlot(pmch_subset,features = "Pmch",raster = F,order = TRUE,reduction=paste0("umap_","from_full_map_good"))

Seurat::DimPlot(pmch_subset,group.by =  "Dataset",raster = F,order = F,shuffle = TRUE,reduction=paste0("umap_",new_name))
Seurat::FeaturePlot(pmch_subset,features = "Pmch",raster = F,order = TRUE,reduction=paste0("umap_",new_name))

pmch_subset = FindNeighbors(pmch_subset,reduction =paste0("umap_","from_full_map_good"),dims = 1:ncol(pmch_subset@reductions[[paste0("umap_","from_full_map_good")]]@cell.embeddings))
pmch_subset = FindClusters(pmch_subset,resolution = 0.15)

Seurat::DimPlot(pmch_subset,group.by =  "seurat_clusters",raster = F,order = F,shuffle = TRUE,label=TRUE,reduction=paste0("umap_","from_full_map_good"))
pmch_sub_markers = FindMarkers(pmch_subset,ident.1 = c("2"),min.pct = 0.1,min.diff.pct = 0.05,only.pos = TRUE)
pmch_sub_markers$gene = rownames(pmch_sub_markers)
pmch_sub_markers$specificity = (pmch_sub_markers$pct.1 / pmch_sub_markers$pct.2) * pmch_sub_markers$avg_log2FC
pmch_sub_markers = pmch_sub_markers[ ! pmch_sub_markers$gene %in% feature_exclude_list,]

#Seurat::FeaturePlot(pmch_subset,features = "Pcsk2",raster = F,order = TRUE,reduction=paste0("umap_",new_name))
Seurat::FeaturePlot(pmch_subset,features = "Parpbp",raster = F,order = TRUE,reduction=paste0("umap_","from_full_map_good"))


####
markers_111 = FindMarkers(hypoMap_seurat_neurons,ident.1 = c("111"),min.pct = 0.1,min.diff.pct = 0.05,only.pos = TRUE,max.cells.per.ident = 5000)
markers_111$gene = rownames(markers_111)
markers_111$specificity = (markers_111$pct.1 / markers_111$pct.2) * markers_111$avg_log2FC
markers_111_sub = markers_111[ ! markers_111$gene %in% feature_exclude_list,]


############
evaluation_purityknn_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_neurons_integration/evaluation/all_purity_knn.txt"
evaluation_purityknn = data.table::fread(evaluation_purityknn_file,data.table = F) %>% dplyr::rename( purityknn = value)
evaluation_purityknn_filtered = evaluation_purityknn %>% dplyr::filter(reduction %in% c("scVI_0_300_0.1_3_256_gene_zinb_cov2..scVI..65..RNA.log.vst.split_Batch_ID.features.1000_562216427a8b35b2658273fecaf712a2",
                                                                                        "scVI_0_400_0.15_4_256_gene_zinb_cov2..scVI..85..RNA.log.vst.split_Batch_ID.features.750_27466e1704ffc1d3f82787db0b6e7b8d",
                                                                                        "scVI_0_400_0.15_3_256_gene_zinb_cov2..scVI..85..RNA.log.vst.split_Batch_ID.features.750_64203fb2c726221727cd01767d952850"))
evaluation_purityknn_filtered_wide = evaluation_purityknn_filtered %>% tidyr::spread(key = reduction, value = purityknn)





