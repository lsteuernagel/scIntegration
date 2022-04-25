


author_celltype_df = hypoMap_seurat@meta.data[,c("Cell_ID","Author_CellType","Dataset")]

author_celltype_df_stats = author_celltype_df %>% dplyr::group_by(Author_CellType) %>% dplyr::count()
author_celltype_df_stats$flagged = "no"
author_celltype_df_stats$flagged[grepl("unassign|other|unknown|Other|Doublets|Mixed|mixed|NA|Unassigned",author_celltype_df_stats$Author_CellType)] = "yes"
author_celltype_df_stats$remove = FALSE
author_celltype_df_stats$remove[author_celltype_df_stats$flagged=="yes"|author_celltype_df_stats$n<10] = TRUE
remove_celltypes = author_celltype_df_stats$Author_CellType[author_celltype_df_stats$remove]

author_celltype_df_filtered = author_celltype_df %>% dplyr::filter(!is.na(Author_CellType) & ! Author_CellType %in% remove_celltypes)

# need to split by dataset

# export two column file

## run asw per uthor cellpye (from each dataset!)
