##########
### map_celltype_signatures
##########

#' Map gene signatures characterizing celltypes to cells in the merged dataset using AUC. Runs AUCell batch wise to avoid large dense matrices ans returns results per cell.
#' @param seurat_object Merged seurat object after computation of reductions (PCA)
#' @param block_size block size for batch wise sequential processing
#' @param aucMaxRank_n how many genes are used for AUC ranking
#' @param gene_set_list set of gene signatures characterizing celltypes
#' @param min_rowSum rowsums per gene when building AUC ranking
#' @param global_seed
#' @return a list of cells per signature.

# set.seed(123)
# exprMatrix = seurat_object@assays$RNA@counts[,sample(1:ncol(seurat_object@assays$RNA@counts),10240)]

map_celltype_signatures2 = function(exprMatrix,block_size=10000,aucMaxRank_n,gene_set_list,min_rowSum=10,global_seed =123){
  #subset to drop some zero genes
  exprMatrix <- exprMatrix[Matrix::rowSums(exprMatrix)>min_rowSum,]

  message(Sys.time(),": Running map_celltype_signatures2 on matrix with ",length(gene_set_list)," signatures" )
  require(tidyverse)
  require(Seurat)
  require(AUCell)
  # # quick fix to this https://github.com/aertslab/AUCell/issues/10
  # library("AUCell",lib.loc = "/beegfs/scratch/bruening_scratch/lsteuernagel/R/user_lib/x86_64-pc-linux-gnu-library/4.0/")
  require(Matrix)

  # get levels to split matrix
  if(ncol(exprMatrix)/block_size > 1){
    cut_levels = as.character(cut(1:ncol(exprMatrix),ncol(exprMatrix)/block_size))
    cut_levels_unique = unique(as.character(cut_levels))
  }else{
    cut_levels = rep("1",ncol(exprMatrix))
    cut_levels_unique = unique(as.character(cut_levels))
  }
  all_signature_auc_mats = list()
  for(i in 1:length(cut_levels_unique)){
    idx = which(cut_levels==cut_levels_unique[i])
    ### Build AUC ranking
    exprMatrix_subset <- exprMatrix[,idx]

    # Build gene-expression rankings for each cell
    message("Building AUC ranking ",i)
    set.seed(global_seed)
    cells_rankings <- AUCell::AUCell_buildRankings(exprMatrix_subset, nCores=1, plotStats=F,verbose = F)

    ### Calculate AUC for signatures and identify cells
    message("Calculate AUC for signatures and extract mapped cells ",i)
    signature_AUC <- AUCell::AUCell_calcAUC(gene_set_list, cells_rankings,nCores=1,verbose=F,aucMaxRank=aucMaxRank_n)
    signature_AUC_mat = t(signature_AUC@assays@data@listData$AUC)
    all_signature_auc_mats[[i]] = signature_AUC_mat
  }
  signature_AUC_mat_full = do.call(rbind,all_signature_auc_mats)
  #gc()
  return(signature_AUC_mat_full)

}
