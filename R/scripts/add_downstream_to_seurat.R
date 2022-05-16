
# save the annoations:
#data.table::fwrite(hypoMap_seurat@meta.data %>% dplyr::select(Cell_ID,SNN_scvi_res.10=seurat_clusters,Author_Class_Curated),"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/hypoMap_v2_preliminary_clusters_090522.tsv",sep="\t")

# save umap
# embed = hypoMap_seurat@reductions$umap_scVI_0_300_0.1_3_256_gene_zinb_cov2..scVI..65..RNA.log.vst.split_Batch_ID.features.3000_5a67c718bb47fffe6e6c8b87f9d1b722@cell.embeddings
# embed = cbind(Cell_ID=rownames(hypoMap_seurat@meta.data),embed)
# data.table::fwrite(as.data.frame(embed),"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/hypoMap_v2_preliminary_scvi300653000.tsv",sep="\t")
#


# load hypoma
hypoMap_seurat = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/hypoMap_merged_filtered.rds")

# load others and add
hypoMap_v2_preliminary_clusters = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/hypoMap_v2_preliminary_clusters_090522.tsv",data.table = F)

# add
temp_meta = dplyr::left_join(hypoMap_seurat@meta.data[,1:33],hypoMap_v2_preliminary_clusters],by=c("Cell_ID"="Cell_ID"))
rownames(temp_meta) = temp_meta$Cell_ID
hypoMap_seurat@meta.data = temp_meta

# load
source("R/evaluation_functions.R")
embedding = read_embedding(filename_withpath = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/hypoMap_v2_preliminary_scvi300653000.tsv",seurat_object = hypoMap_seurat)

#add
new_name = "umap_scvi"#"scvi"
dimred <- Seurat::CreateDimReducObject(
  embeddings = as.matrix(embedding),
  stdev = as.numeric(apply(embedding, 2, stats::sd)),
  assay = "RNA",
  key = new_name
)
# add
hypoMap_seurat@reductions[[new_name]] = dimred

# ready for plots:
# plot
p3= Seurat::DimPlot(hypoMap_seurat,group.by =  "Author_Class_Curated",raster = F,order = TRUE,shuffle = TRUE,reduction=paste0("umap_",new_name))#+NoLegend()#reduction=paste0("umap_",new_name)
p3r = scUtils::rasterize_ggplot(p3,pixel_raster = 2048)
p3r



