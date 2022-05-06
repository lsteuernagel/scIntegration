
##########
### Author Anno v2 full map
##########

## export all author annotations

author_celltype_df = hypoMap_seurat@meta.data[,c("Cell_ID","Author_CellType","Dataset")]

author_celltype_df_stats = author_celltype_df %>% dplyr::group_by(Author_CellType) %>% dplyr::count()
author_celltype_df_stats$flagged = "no"
author_celltype_df_stats$flagged[grepl("unassign|other|unknown|Other|Doublets|Mixed|mixed|NA|Unassigned",author_celltype_df_stats$Author_CellType)] = "yes"
author_celltype_df_stats$remove = FALSE
author_celltype_df_stats$remove[author_celltype_df_stats$flagged=="yes"|author_celltype_df_stats$n<10] = TRUE
remove_celltypes = author_celltype_df_stats$Author_CellType[author_celltype_df_stats$remove]

author_celltype_df_filtered = author_celltype_df %>% dplyr::filter(!is.na(Author_CellType) & ! Author_CellType %in% remove_celltypes)

# need to split by dataset
author_celltype_df_filtered_list = split(author_celltype_df_filtered[,1:2],f = author_celltype_df_filtered$Dataset)

# export two column file
label_folder ="/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/cell_labels/"
system(paste0("mkdir -p ",label_folder))
for( i in 1:length(author_celltype_df_filtered_list)){
  data.table::fwrite(author_celltype_df_filtered_list[[i]],file = paste0(label_folder,names(author_celltype_df_filtered_list)[i],"_AuthorCellTypes.csv"))
}

## run asw per Author cellpye (from each dataset!)

##########
### Author Anno v1 full map
##########

### using v1:

author_celltype_df_v1 = hypoMap_seurat_neurons_v1@meta.data[,c("Cell_ID","Author_CellType","Dataset")]

author_celltype_df_v1_stats = author_celltype_df_v1 %>% dplyr::group_by(Author_CellType) %>% dplyr::count()
author_celltype_df_v1_stats$flagged = "no"
author_celltype_df_v1_stats$flagged[grepl("unassign|other|unknown|Other|Doublets|Mixed|mixed|NA|Unassigned",author_celltype_df_v1_stats$Author_CellType)] = "yes"
author_celltype_df_v1_stats$remove = FALSE
author_celltype_df_v1_stats$remove[author_celltype_df_v1_stats$flagged=="yes"|author_celltype_df_v1_stats$n<20] = TRUE
remove_celltypes = author_celltype_df_v1_stats$Author_CellType[author_celltype_df_v1_stats$remove]

##
author_celltype_df_full_v1 = hypoMap_seurat_v1@meta.data[,c("Cell_ID","Author_CellType","Dataset")]
author_celltype_df_full_v1_stats = author_celltype_df_full_v1 %>% dplyr::group_by(Author_CellType) %>% dplyr::count(name="n_all")

author_celltype_df_v1_stats = dplyr::left_join(author_celltype_df_v1_stats,author_celltype_df_full_v1_stats,by=c("Author_CellType"="Author_CellType"))

author_celltype_df_v1_filtered = author_celltype_df_v1 %>% dplyr::filter(!is.na(Author_CellType) & ! Author_CellType %in% remove_celltypes)
author_celltype_df_v1_filtered =  author_celltype_df_v1_filtered %>% dplyr::filter(Cell_ID %in% author_celltype_df_full_v1$Cell_ID)

# need to split by dataset
author_celltype_df_v1_filtered_list = split(author_celltype_df_v1_filtered[,1:2],f = author_celltype_df_v1_filtered$Dataset)

# export two column file
label_folder ="/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_integration/cell_labels/"
system(paste0("mkdir -p ",label_folder))
for( i in 1:length(author_celltype_df_v1_filtered_list)){
  data.table::fwrite(author_celltype_df_v1_filtered_list[[i]],file = paste0(label_folder,names(author_celltype_df_v1_filtered_list)[i],"_AuthorCellTypes.csv"))
}

##########
### Author Anno v2 neuron map
##########

hypoMap_seurat_neurons = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_neurons_integration/hypoMap_neurons.rds")

## export all author annotations

author_celltype_df = hypoMap_seurat_neurons@meta.data[,c("Cell_ID","Author_CellType","Dataset")]

author_celltype_df_stats = author_celltype_df %>% dplyr::group_by(Author_CellType) %>% dplyr::count()
author_celltype_df_stats$flagged = "no"
author_celltype_df_stats$flagged[grepl("unassign|other|unknown|Other|Doublets|Mixed|mixed|NA|Unassigned|Astro|Oligo|Tany|Ependy|Endoth|Mural",author_celltype_df_stats$Author_CellType)] = "yes"
author_celltype_df_stats$remove = FALSE
author_celltype_df_stats$remove[author_celltype_df_stats$flagged=="yes"|author_celltype_df_stats$n<10] = TRUE
remove_celltypes = author_celltype_df_stats$Author_CellType[author_celltype_df_stats$remove]

author_celltype_df_filtered = author_celltype_df %>% dplyr::filter(!is.na(Author_CellType) & ! Author_CellType %in% remove_celltypes)

# need to split by dataset
author_celltype_df_filtered_list = split(author_celltype_df_filtered[,1:2],f = author_celltype_df_filtered$Dataset)

# export two column file
label_folder ="/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_neurons_integration/cell_labels/"
system(paste0("mkdir -p ",label_folder))
for( i in 1:length(author_celltype_df_filtered_list)){
  data.table::fwrite(author_celltype_df_filtered_list[[i]],file = paste0(label_folder,names(author_celltype_df_filtered_list)[i],"_AuthorCellTypes.csv"))
}


