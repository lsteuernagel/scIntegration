
integration_file_name = "scVI_0_200_0.01_4_256_gene_zinb_cov2..scVI..50..RNA.log.vst.split_Batch_ID.features.1500_e06093f25a8fc603ec7f9eb5186b44f8.txt"

integration_file_name = "scVI_0_400_0.05_2_256_gene_zinb_cov2..scVI..110..RNA.log.vst.split_Batch_ID.features.2000_2d15d1b4808a8f4d44f3e609493a7c3a.txt"

integration_file_name = "scVI_0_300_0.05_3_256_gene_zinb_cov2..scVI..110..RNA.log.vst.split_Batch_ID.features.1250_fb23407b78b330da59bc5a1e8f539053.txt"

integration_file_name = "scVI_0_400_0.1_4_256_gene_zinb_cov2..scVI..65..RNA.log.vst.split_Batch_ID.features.3000_e9b11dacd25decfcd89c1623473e5352.txt"

integration_file_name = "scVI_0_300_0.01_3_256_gene_zinb_cov2..scVI..85..RNA.log.vst.split_Batch_ID.features.1500_99244d8dce5a5588784f504eab99f0eb.txt"

# high label asw + author label asw:
integration_file_name = "scVI_0_400_0.05_4_256_gene_zinb_cov2..scVI..110..RNA.log.vst.split_Batch_ID.features.2000_7e37d6e7f4cd4d3603bfc40f122cdf2b.txt"
# relatively high author label asw, also high in chen and moffit. also higher than other asw in purity asw
integration_file_name = "scVI_1_400_0.1_2_256_gene_zinb_cov2..scVI..85..RNA.log.vst.split_Batch_ID.features.1250_03c10730a4b1685057ad0314b7f86cf4.txt"

## good knn mixing and purity withd decent asw:
integration_file_name = "scVI_1_400_0.1_2_256_gene_zinb_cov2..scVI..85..RNA.log.vst.split_Batch_ID.features.2000_5e31aa323c46e0fb15d6f94afc210cce.txt"

#
integration_file_name = "scVI_1_400_0.1_2_256_gene_zinb_cov2..scVI..85..RNA.log.vst.split_Batch_ID.features.2000_5e31aa323c46e0fb15d6f94afc210cce.txt"

integration_file_name = "scVI_1_400_0.1_2_256_gene_zinb_cov2..scVI..85..RNA.log.vst.split_Batch_ID.features.2000_5e31aa323c46e0fb15d6f94afc210cce.txt"

integration_file_name = "scVI_0_300_0.1_3_256_gene_zinb_cov2..scVI..65..RNA.log.vst.split_Batch_ID.features.3000_5a67c718bb47fffe6e6c8b87f9d1b722.txt"

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
p1r = scUtils::rasterize_ggplot(p1,pixel_raster = 4096)
p1r


# plot
p2= Seurat::FeaturePlot(hypoMap_seurat,features = "P2ry13",raster = F,order = TRUE,reduction=paste0("umap_",new_name))#reduction=paste0("umap_",new_name)
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


###################

integration_file_name = "scVI_0_400_0.1_4_256_gene_zinb_cov2..scVI..65..RNA.log.vst.split_Batch_ID.features.3000_e9b11dacd25decfcd89c1623473e5352.txt"
new_name = gsub(".txt","",integration_file_name )#"scvi"
clusterRes = 10
# hypoMap_seurat <- Seurat::FindNeighbors(hypoMap_seurat,reduction = new_name, dims = 1:ncol(hypoMap_seurat@reductions[[new_name]]@cell.embeddings),k.param = 30,verbose = F,graph.name="SNN_scvi")
# hypoMap_seurat <- Seurat::FindClusters(hypoMap_seurat,graph.name="SNN_scvi", resolution = clusterRes,verbose = F,random.seed = 123456)

# plot
p3= Seurat::DimPlot(hypoMap_seurat,group.by =  "seurat_clusters",raster = F,order = TRUE,shuffle = TRUE,reduction=paste0("umap_",new_name),label=TRUE)+NoLegend()#reduction=paste0("umap_",new_name)
p3r = scUtils::rasterize_ggplot(p3,pixel_raster = 2048)
p3r

#data.table::fwrite(hypoMap_seurat_meta[,c("Cell_ID","SNN_scvi_res.10")],"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/hypoMap_v2_perliminary_clusters_290422.tsv",sep="\t")

hypoMap_seurat_meta = hypoMap_seurat@meta.data

hypoMap_seurat_meta = hypoMap_seurat@meta.data[,1:33]
hypoMap_seurat_meta = dplyr::left_join(hypoMap_seurat_meta, hypoMap_v2_perliminary_clusters[,c("Cell_ID","SNN_scvi_res.10")],by=c("Cell_ID"="Cell_ID"))
hypoMap_seurat_meta$seurat_clusters = as.character(hypoMap_seurat_meta$SNN_scvi_res.10)
hypoMap_seurat_meta$Author_Class[hypoMap_seurat_meta$Author_Class=="Vascular"] ="Endothelial"
hypoMap_seurat_meta$Author_Class[hypoMap_seurat_meta$Author_Class %in% c("Unassigned","Mixed")] = NA
class_per_cluster = hypoMap_seurat_meta %>% dplyr::group_by(seurat_clusters) %>% dplyr::add_count(name="cluster_total") %>%
  dplyr::group_by(seurat_clusters,Author_Class,cluster_total)%>%
  dplyr::count(name="cluster_class_count") %>% dplyr::mutate(pct = round(cluster_class_count / cluster_total *100,4)) %>%
  dplyr::filter(!is.na(Author_Class)) %>% dplyr::group_by(seurat_clusters) %>% dplyr::mutate(annotated_pct = sum(pct)) %>%
  dplyr::ungroup() %>%  dplyr::filter(pct > 1)

class_per_cluster_top = class_per_cluster %>% dplyr::group_by(seurat_clusters)  %>% dplyr::filter(pct == max(pct)) %>%
  dplyr::mutate(fraction_of_annotated = pct / annotated_pct) %>% dplyr::distinct(seurat_clusters,.keep_all = TRUE)

class_per_cluster_top$Author_Class_Curated = class_per_cluster_top$Author_Class
class_per_cluster_top$Author_Class_Curated[class_per_cluster_top$fraction_of_annotated < 0.67] = "Unknown"
class_per_cluster_top$Author_Class_Curated[class_per_cluster_top$fraction_of_annotated < 0.8 & class_per_cluster_top$annotated_pct < 10] = "Unknown"
#class_per_cluster_top$seurat_clusters= as.factor(class_per_cluster_top$seurat_clusters)

temp_meta = dplyr::left_join(hypoMap_seurat_meta,class_per_cluster_top[,c("seurat_clusters","Author_Class_Curated")],by=c("seurat_clusters"="seurat_clusters"))
rownames(temp_meta) = temp_meta$Cell_ID
temp_meta$Author_Class_Curated[is.na(temp_meta$Author_Class_Curated)] ="Unknown"
hypoMap_seurat@meta.data = temp_meta
#

# plot
p3= Seurat::DimPlot(hypoMap_seurat,group.by =  "Author_Class_Curated",raster = F,order = TRUE,shuffle = TRUE,reduction=paste0("umap_",new_name))#+NoLegend()#reduction=paste0("umap_",new_name)
p3r = scUtils::rasterize_ggplot(p3,pixel_raster = 2048)
p3r

# 25

# 83 mural

# 126 fibroblasts

# what are 187: Ppbp+ and 195: weird lots of gm and rpl genes

# 117 looks problematic

clusterid="142"
table(hypoMap_seurat_meta$Dataset[hypoMap_seurat_meta$SNN_scvi_res.10==clusterid])

table(hypoMap_seurat_meta$Author_CellType[hypoMap_seurat_meta$SNN_scvi_res.10==clusterid])[table(hypoMap_seurat_meta$Author_CellType[hypoMap_seurat_meta$SNN_scvi_res.10==clusterid])>10]


# plot
cellsh=hypoMap_seurat_meta$Cell_ID[hypoMap_seurat_meta$SNN_scvi_res.10==clusterid]
p3= Seurat::DimPlot(hypoMap_seurat,group.by =  "Author_Class_Curated",cells.highlight = cellsh,sizes.highlight = 0.1,raster = F,order = TRUE,shuffle = TRUE,reduction=paste0("umap_",new_name))#+NoLegend()#reduction=paste0("umap_",new_name)
p3r = scUtils::rasterize_ggplot(p3,pixel_raster = 2048)
p3r


# plot
p3= Seurat::DimPlot(hypoMap_seurat,group.by =  "seurat_clusters",raster = F,order = TRUE,shuffle = TRUE,reduction=paste0("umap_",new_name),label=TRUE)+NoLegend()#reduction=paste0("umap_",new_name)
p3r = scUtils::rasterize_ggplot(p3,pixel_raster = 2048)
p3r

#################

hypoMap_celltype_auc_per_cell_result = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/hypoMap_celltype_auc_per_cell_result.txt",data.table = F)
hypoMap_celltype_auc_per_cell_result = cbind(hypoMap_seurat@meta.data[,c("Cell_ID","Dataset","seurat_clusters")],hypoMap_celltype_auc_per_cell_result) %>% as.data.frame()
hypoMap_celltype_auc_per_cell_result_wide = hypoMap_celltype_auc_per_cell_result %>% tidyr::gather(key = "celltype",value = "auc",-Cell_ID,-Dataset,-seurat_clusters)#dplyr::group_by()
hypoMap_celltype_auc_per_cell_result_wide_stat = hypoMap_celltype_auc_per_cell_result_wide %>% dplyr::group_by(seurat_clusters,celltype) %>%
  dplyr::summarise(mean_auc= mean(auc))
hypoMap_celltype_auc_per_cell_result_wide_stat$seurat_clusters = as.factor(hypoMap_celltype_auc_per_cell_result_wide_stat$seurat_clusters)



# plot
p3= Seurat::DimPlot(hypoMap_seurat,group.by =  "seurat_clusters",raster = F,order = TRUE,shuffle = TRUE,reduction="umap",label=TRUE)+NoLegend()#reduction=paste0("umap_",new_name)
p3r = scUtils::rasterize_ggplot(p3,pixel_raster = 2048)
p3r

### subset for marker detection:
options(scipen=5) # prevent scientific notation
hypoMap_seurat_downsample_meta = scUtils::downsample_balanced_iterative(hypoMap_seurat@meta.data,target_sample_size = 40000,predictor_var = "seurat_clusters",stepsize = 500,global_seed = 123456)

hypoMap_seurat_subset = subset(hypoMap_seurat,cells = hypoMap_seurat_downsample_meta$Cell_ID)

marker_genes = FindMarkers(hypoMap_seurat_subset,ident.1 = "142",logfc.threshold = 0.3,min.pct = 0.1,min.diff.pct = 0.05,max.cells.per.ident = 5000)
marker_genes$gene = rownames(marker_genes)
marker_genes$specificity = (marker_genes$pct.1 / marker_genes$pct.2) * marker_genes$avg_log2FC

# Mobp
Seurat::FeaturePlot(hypoMap_seurat_subset,features = "Mobp",raster = F,order = TRUE,reduction=paste0("umap_",new_name))#reduction=paste0("umap_",new_name)

# find Gnrh1 neurons
aa1= cbind(seurat_clusters=hypoMap_seurat_subset@meta.data$seurat_clusters,FetchData(hypoMap_seurat_subset,vars = "Gnrh1"))


# plot
p2= Seurat::FeaturePlot(hypoMap_seurat,features = "Top2a",raster = F,order = TRUE,reduction=paste0("umap_",new_name))#reduction=paste0("umap_",new_name)
p2r = scUtils::rasterize_ggplot(p2,pixel_raster = 2048)
p2r

####### manually add annoation
class_per_cluster_top[class_per_cluster_top$Author_Class_Curated=="Unknown",]
hypoMap_seurat@meta.data$Author_Class_Curated[hypoMap_seurat@meta.data$seurat_clusters=="14"] = "Astrocytes"
hypoMap_seurat@meta.data$Author_Class_Curated[hypoMap_seurat@meta.data$seurat_clusters=="82"] = "Tanycytes"
hypoMap_seurat@meta.data$Author_Class_Curated[hypoMap_seurat@meta.data$seurat_clusters=="83"] = "Mural"
hypoMap_seurat@meta.data$Author_Class_Curated[hypoMap_seurat@meta.data$seurat_clusters=="126"] = "Fibroblast"
hypoMap_seurat@meta.data$Author_Class_Curated[hypoMap_seurat@meta.data$seurat_clusters=="127"] = "Astrocytes"
hypoMap_seurat@meta.data$Author_Class_Curated[hypoMap_seurat@meta.data$seurat_clusters=="142"] = "Differentiating"
hypoMap_seurat@meta.data$Author_Class_Curated[hypoMap_seurat@meta.data$seurat_clusters=="145"] = "Doublet"
hypoMap_seurat@meta.data$Author_Class_Curated[hypoMap_seurat@meta.data$seurat_clusters=="167"] = "Doublet"
hypoMap_seurat@meta.data$Author_Class_Curated[hypoMap_seurat@meta.data$seurat_clusters=="168"] = "Doublet"
hypoMap_seurat@meta.data$Author_Class_Curated[hypoMap_seurat@meta.data$seurat_clusters=="177"] = "Doublet"
hypoMap_seurat@meta.data$Author_Class_Curated[hypoMap_seurat@meta.data$seurat_clusters=="184"] = "Neurons"
hypoMap_seurat@meta.data$Author_Class_Curated[hypoMap_seurat@meta.data$seurat_clusters=="185"] = "Doublet"
hypoMap_seurat@meta.data$Author_Class_Curated[hypoMap_seurat@meta.data$seurat_clusters=="191"] = "Cd52+ immune"
hypoMap_seurat@meta.data$Author_Class_Curated[hypoMap_seurat@meta.data$seurat_clusters=="195"] = "Unknown"
hypoMap_seurat@meta.data$Author_Class_Curated[hypoMap_seurat@meta.data$seurat_clusters=="187"] = "Krt19+ cells"

# plot
p3= Seurat::DimPlot(hypoMap_seurat,group.by =  "Author_Class_Curated",raster = F,order = TRUE,shuffle = TRUE,reduction=paste0("umap_",new_name),label=TRUE)#+NoLegend()#reduction=paste0("umap_",new_name)
p3r = scUtils::rasterize_ggplot(p3,pixel_raster = 2048)
p3r

#data.table::fwrite(hypoMap_seurat_meta[,c("Cell_ID","SNN_scvi_res.10","Author_Class_Curated")],"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/hypoMap_v2_perliminary_clusters_290422.tsv",sep="\t")
#data.table::fwrite(hypoMap_seurat@meta.data %>% dplyr::select(Cell_ID,SNN_scvi_res.10=seurat_clusters,Author_Class_Curated),"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/hypoMap_v2_preliminary_clusters_090522.tsv",sep="\t")

########## export
v2_neurons_folder ="/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_neurons_integration/"

hypoMap_seurat_neurons = subset(hypoMap_seurat,subset=Author_Class_Curated == "Neurons")
hypoMap_seurat_neurons@reductions = list()
hypoMap_seurat_neurons
saveRDS(hypoMap_seurat_neurons,paste0(v2_neurons_folder,"hypoMap_neurons.rds"))



p2= Seurat::FeaturePlot(hypoMap_seurat,features = "Pnoc",raster = F,order = TRUE,reduction=paste0("umap_",new_name))#reduction=paste0("umap_",new_name)
p2r = scUtils::rasterize_ggplot(p2,pixel_raster = 2048)
p2r
