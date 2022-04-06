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

##########
### read_embedding
##########

#' Load an emebedding with cells x lowDims from flatfile, ensuring consistency with a Seurat object (or metadata only for faster usage)
#' @param filename_withpath filepath
#' @param seurat_object seuratobject associated with current embedding. If specified metadata does not have to be set explicitly.
#' @param seurat_object_metadata metadata only of seuratobject associated with current embedding
#' @return

read_embedding = function(filename_withpath,seurat_object=NULL,seurat_object_metadata=NULL){

  #get metadata
  if(!is.null(seurat_object_metadata)){
    metadata = seurat_object_metadata
  }else{
    if(!is.null(seurat_object)){
      metadata = seurat_object@meta.data
    }else{
      stop("Please provide either a dataframe with metadata or a seurat object with metadata that can be exctracted!")
    }
  }
  # load
  current_embedding = data.table::fread(filename_withpath,data.table = F)
  # use first col as rownames
  if(is.character(current_embedding[,1])){
    rnames = current_embedding[,1]
    current_embedding = current_embedding[,2:ncol(current_embedding)]
    rownames(current_embedding)=rnames
    # reorder to align with rest of object
    if(any(is.na(match(rownames(metadata),rownames(current_embedding))))){
      stop("Cell names from loaded reduction and new object are not matching exactly. Stopping import.")
    }
    current_embedding = current_embedding[match(rownames(metadata),rownames(current_embedding)),]
  }else{
    warning("First column of loaded file is not of type character, using rownames of metadata as rownames of added reduction. This can induce bugs if the order changed due to split/merge of the Seurat object!")
    rownames(current_embedding) = rownames(metadata)
  }
  return(current_embedding)
}


##########
### evaluate_purity_knn
##########

#' Wrapper that runs kNN based evaluation of embeddings for celltype purity.
#' @param seurat_object_metadata metadata associated with the mapped celltypes and emebdding that will be evaluated
#' @param cells_sets list of vectors (per mapped celltypes) with celllabels matching rownames of metadata
#' @param k_param cell ids of subset that is used for evaluation
#' @param control_size_factor_max
#' @param integration_files names of files that contain embeddings
#' @param integration_path directory where to find integration_files
#' @param evaluation_file filename where to store results
#' @param ncores how many cores to use by doParallel
#' @param dist_type all batches with a probability > max_for_norm are considered when normalizing
#' @param global_seed seed
#' @return Nothing. Writes files

evaluate_purity_knn = function(seurat_object_metadata,cells_sets,k_param=20,integration_files,integration_path,evaluation_file,ncores=1,dist_type="cosine",global_seed=123){
  # load
  require(randomForest)
  require(foreach)
  require(doParallel)
  require(doRNG)
  # check that files can be found at integration path
  # make names
  integration_to_run=c()
  for( j in 1:length(integration_files)){
    current_file = as.character(integration_files[j])
    split_filepath = as.character(strsplit(current_file,split = "/")[[1]])
    if(length(split_filepath)>1){
      current_file = gsub(".txt","",split_filepath[length(split_filepath)])
    }else{
      current_file = gsub(".txt","",current_file)
    }
    integration_to_run = c(integration_to_run,current_file)
  }
  #init result lists
  celltype_purity_knn_list = list()
  # parallel
  registerDoParallel(cores=ncores)
  doRNG::registerDoRNG(seed = global_seed)
  res_par_1 <- foreach(embed_idx = 1:length(integration_to_run), .combine='cbind') %dopar% {
    # get current result
    current_name = gsub(".txt","",integration_to_run[embed_idx])
    message(current_name)
    current_file = integration_files[which(grepl(paste0(integration_to_run[embed_idx],".txt"),integration_files))][1]
    current_embedding = read_embedding(paste0(integration_path,current_file),seurat_object_metadata = seurat_object_metadata)
    current_embedding = as.matrix(current_embedding)
    ## knn based:
    ## get all neighbors
    message("Start NN")
    n_dim=ncol(current_embedding)
    neighbor_dists = run_seurat_annoy_manually(current_embedding,dim=n_dim,k=k_param+1,metric=dist_type)
    #neighbor_dists = run_seurat_annoy_manually(current_embedding,dim=20,k=20+1,metric="cosine")
    nn_idx=neighbor_dists$nn.idx
    message("Built NN")
    # init
    all_occ_scores = c()
    # for all cell types --> ids in cells_sets must correspond to rownames
    for(j in 1:length(cells_sets)){
      # which idx are in given cell set
      pos_idx = which(seurat_object_metadata$Cell_ID %in% cells_sets[[j]])
      ## subset of cells belonging to cell types
      pos_nn_idx = nn_idx[pos_idx,2:ncol(nn_idx)]
      # set to 1 if neighbor is also in cell types, 0 otherwise
      pos_nn_scores = apply(pos_nn_idx,2,function(indices,pos_indices){y=rep(0,length(indices));y[indices %in% pos_indices]=1;return(y)},pos_indices=pos_idx)
      # calculate the fraction of neighbors belonging to the same celltype
      pos_nn_occ = apply(pos_nn_scores,1,sum) / ncol(pos_nn_scores)
      all_occ_scores = c(all_occ_scores,mean(pos_nn_occ)) # add the mean per celltype ! (not all cell values!)
    }
    message("Calculated knn purity")
    #store per cell type and return
    names(all_occ_scores) = names(cells_sets)
    all_occ_scores
  }
  message(paste0(dim(res_par_1)))
  colnames(res_par_1) = integration_to_run
  # save updated files
  message("Saving results to files.")
  data.table::fwrite(as.data.frame(res_par_1),file = evaluation_file,sep =  "\t")
}

##########
### evaluate_knn_entropy_all --> mixing evaluation using kNN
##########

#' Wrapper that runs kNN based evaluation of embeddings for mixing.
#' @param seurat_object_metadata metadata associated with the mapped celltypes and emebdding that will be evaluated
#' @param cells_sets list of vectors (per mapped celltypes) with celllabels matching rownames of metadata
#' @param k_param cell ids of subset that is used for evaluation
#' @param control_size_factor_max
#' @param integration_files names of files that contain embeddings
#' @param integration_path directory where to find integration_files
#' @param evaluation_file filename where to store results
#' @param ncores how many cores to use by doParallel
#' @param dist_type all batches with a probability > max_for_norm are considered when normalizing
#' @param global_seed seed
#' @return Nothing. Writes files

evaluate_mixing_knn = function(seurat_object_metadata,batch_var="Batch_ID",k_param=20,integration_files,integration_path,evaluation_file,ncores=1,dist_type="cosine",global_seed=123){

  require(randomForest)
  require(foreach)
  require(doParallel)
  require(doRNG)
  message("evaluate_knn_entropy")
  # check that files can be found at integration path

  # make names
  #   split_filepath = sapply(integration_files,function(x){as.character(strsplit(x,split = "/")[[1]])})
  integration_to_run=c()
  for( j in 1:length(integration_files)){
    current_file = as.character(integration_files[j])
    split_filepath = as.character(strsplit(current_file,split = "/")[[1]])
    if(length(split_filepath)>1){
      current_file = gsub(".txt","",split_filepath[length(split_filepath)])
    }else{
      current_file = gsub(".txt","",current_file)
    }
    integration_to_run = c(integration_to_run,current_file)
  }
  #init result lists
  celltype_purity_knn_list = list()
  #celltypeProb_purity_norm_list  = list()
  # run classProb on all missing results
  #for(i in 1:length(integration_to_run)){
  # parallel
  registerDoParallel(cores=ncores)
  doRNG::registerDoRNG(seed = global_seed)
  #for(embed_idx in 1:2){
  res_par_1 <- foreach(embed_idx = 1:length(integration_to_run), .combine='cbind') %dopar% {
    # get current result
    current_name = gsub(".txt","",integration_to_run[embed_idx])
    message(current_name)
    current_file = integration_files[which(grepl(paste0(integration_to_run[embed_idx],".txt"),integration_files))][1]
    # message(">>>>>",paste0(integration_path,current_file))
    current_embedding = read_embedding(paste0(integration_path,current_file),seurat_object_metadata = seurat_object_metadata)
    current_embedding = as.matrix(current_embedding)
    # set cell names --> assumes that the order of cells in metadata is the same as in embedding!
    #rownames(current_embedding) = rownames(seurat_object_metadata)

    ## knn based:
    ## get all neighbors
    message("Start NN")
    n_dim=ncol(current_embedding)
    neighbor_dists = run_seurat_annoy_manually(current_embedding,dim=n_dim,k=k_param+1,metric=dist_type)
    #neighbor_dists = run_seurat_annoy_manually(current_embedding,dim=20,k=20+1,metric="cosine")
    nn_idx=neighbor_dists$nn.idx
    message("Built NN")

    cell_labels = as.character(seurat_object_metadata[,batch_var])

    pos_nn_idx = nn_idx[,2:ncol(nn_idx)] # get idx for current cells
    entropy = apply(pos_nn_idx,1,function(row,cell_labels){ # run entropy on batch distribution for celltype cells
      freq_batch = table(cell_labels[row])/length(cell_labels[row])
      freq_batch = freq_batch[freq_batch > 0]
      entropy = entropy_fun(freq_batch,logfun ="log2")
      entropy
    }, cell_labels = cell_labels)
    all_entropy_scores = entropy / max(1,log2(length(unique(cell_labels)))) # normalize by this number so that max entropy =1
    all_entropy_scores
  }
  message(paste0(dim(res_par_1)))
  colnames(res_par_1) = integration_to_run
  # save updated files
  message("Saving results to files.")
  data.table::fwrite(as.data.frame(res_par_1),file =  evaluation_file,sep =  "\t")

}

##########
### evaluate_mixing_rf
##########

#' Wrapper that runs random forest based evaluation of embeddings for mixing.
#' @param seurat_object_metadata metadata associated with the mapped celltypes and emebdding that will be evaluated
#' @param integration_files names of files that contain embeddings
#' @param integration_path directory where to find integration_files
#' @param evaluation_file filename where to store results
#' @param batch_var name of column in metadata
#' @param subset_cell_ids cell ids of subset that is used for evaluation
#' @param ntrees Number of trees, should be rather high to obtain accurate oob probabilities
#' @param sampsize_pct pct of samples in bag
#' @param ncores how many cores to use by doParallel
#' @param returnNormalized whether to write normalized or non-normalized
#' @param returnCollapsed whether to take the mean before writing results
#' @param scale_to_max scale non-normalized entropy to maximum number of batches
#' @param max_for_norm all batches with a probability > max_for_norm are considered when normalizing
#' @param global_seed seed
#' @return Nothing. Writes files

evaluate_mixing_rf = function(seurat_object_metadata,integration_files,integration_path,evaluation_file,batch_var,subset_cell_ids=NULL,ntrees,sampsize_pct,ncores,scale_to_max,returnNormalized=TRUE,returnCollapsed =TRUE,max_for_norm,global_seed){

  # check that files can be found at integration path
  message("meta_data ",dim(seurat_object_metadata)[1])
  # make names
  #   split_filepath = sapply(integration_files,function(x){as.character(strsplit(x,split = "/")[[1]])})
  integration_to_run=c()
  for( j in 1:length(integration_files)){
    current_file = as.character(integration_files[j])
    #message(current_file)
    split_filepath = as.character(base::strsplit(current_file,split = "\\/")[[1]])
    if(length(split_filepath)>1){
      current_file = gsub(".txt","",split_filepath[length(split_filepath)])
    }else{
      current_file = gsub(".txt","",current_file)
    }
    integration_to_run = c(integration_to_run,current_file)
  }
  # prepare other info
  batch_labels = as.character(seurat_object_metadata[,batch_var])
  message("check nas ",any(is.na(batch_labels)))

  #init result lists
  classProb_entropy_list = list()
  classProb_entropy_norm_list = list()

  # run classProb on all missing results
  if(length(integration_to_run)>0){
    for(i in 1:length(integration_to_run)){
      # get current result
      current_name = gsub(".txt","",integration_to_run[i])
      current_file = integration_files[which(grepl(paste0(integration_to_run[i],".txt"),integration_files))]
      current_embedding = read_embedding(paste0(integration_path,current_file),seurat_object_metadata=seurat_object_metadata)
      message("Rownames" ,length(rownames(current_embedding))," ",rownames(current_embedding)[1]," ")
      # set cell names --> assumes that the order of cells in metadata is the same as in embedding!
      if(!is.null(subset_cell_ids)){
        idx = which(rownames(current_embedding) %in% subset_cell_ids)
        current_embedding = current_embedding[idx,]
        batch_labels_subset = batch_labels[idx]
      }else{
        batch_labels_subset = batch_labels
      }
      # entropy on class probabilities
      message(Sys.time(),"Evaluating ",current_name)
      classProb_temp= classProb(train_predictors=current_embedding,
                                train_response=batch_labels_subset,
                                trees=ntrees,
                                sampsize_pct = sampsize_pct,
                                ncores=ncores,
                                scale_to_max=scale_to_max,
                                max_for_norm = max_for_norm,
                                global_seed =global_seed)

      #store results
      classProb_entropy_list[[current_name]] = classProb_temp$entropy
      classProb_entropy_norm_list[[current_name]] = classProb_temp$entropy_norm
    }
  }

  # save updated files
  # only attempt saving if results are there
  if(length(classProb_entropy_list)>0){
    # combine results
    classProb_entropy = do.call(cbind,classProb_entropy_list)
    classProb_entropy_norm = do.call(cbind,classProb_entropy_norm_list)
    message("Saving results.")
    # only write normalized OR raw
    if(returnNormalized){
      # optionally collapse to median
      if(returnCollapsed){
        classProb_entropy_norm = data.frame(mixing = apply(classProb_entropy_norm,2,median),reduction=as.character(names(apply(classProb_entropy_norm,2,median))))
      }
      data.table::fwrite(classProb_entropy_norm, file = evaluation_file, sep ="\t")
    }else{
      if(returnCollapsed){
        classProb_entropy = data.frame(mixing = apply(classProb_entropy,2,median),reduction=as.character(names(apply(classProb_entropy,2,median))))
      }
      data.table::fwrite(classProb_entropy, file = evaluation_file, sep ="\t")
    }
  }
}

##########
### classProb
##########

#' Calculate entropy of a vector of probabilities, setting log(0) = 0
#' @param x Merged seurat object after computation of reductions (PCA)
#' @param logfun which function for log ?
#' @return shannon entropy

entropy_fun = function(x,logfun ="log2"){
  log_vec = do.call(logfun,list(x))
  log_vec[is.infinite(log_vec)] = 0
  log_vec[is.nan(log_vec)] = 0
  return(-sum(x * log_vec))
}

#' Train a random forest and obtain oob probabilities per batch (class). Then report the entropy per cell as a measure of batch mixing.
#' Additionally reports a normalized entropy which accounts for batches where the celltype is not expected (via cutoff)
#' @param train_predictors cell x features (genes,PCs,etc)
#' @param train_response batch labels as classes
#' @param trees Number of trees, should be rather high to obtain accurate oob probabilities
#' @param sampsize_pct percentage of samples in bag
#' @param ncores how many cores to use by doParallel
#' @param scale_to_max scale non-normalized entropy to maximum number of batches
#' @param max_for_norm all batches with a probability > max_for_norm are considered when normalizing
#' @param global_seed seed
#' @return list of two vectors: entropy and entropy normalized with log2(expected number of batches)

classProb = function(train_predictors,train_response,trees=500,sampsize_pct=0.632,ncores=4,scale_to_max = FALSE,max_for_norm = 0.01,global_seed = 123){
  require(randomForest)
  require(foreach)
  require(doParallel)
  require(doRNG)
  # length of current_labels
  n_batches = length(unique(train_response))
  # run random forest
  registerDoParallel(cores=ncores)
  registerDoRNG(seed = global_seed)
  rf_res <- foreach(ntree=rep(trees/ncores, ncores), .combine=randomForest::combine,.multicombine=TRUE, .packages='randomForest') %dopar% {
    randomForest(x=train_predictors, y=as.factor(train_response), ntree=ntree,sampsize=ceiling(sampsize_pct*nrow(train_predictors)))
  }
  # obtain class probabilities
  votes = rf_res$votes / ncores
  # use entropy to quantify batch effects --> still descirbes situation with n_batches, NOT the current clustering
  entropy = apply(votes,1,entropy_fun)
  # scale to number of expected batches per cell
  norm_factor = log2(apply(votes,1,function(x,max_n){return(length(x[x>max_n]))},max_n=max_for_norm))
  norm_factor[norm_factor<1] = 1 # 1 as minimum
  entropy_norm = entropy / norm_factor
  # scale if told to
  if(scale_to_max){
    entropy = entropy / log2(ncol(votes))
  }
  return(list(entropy = entropy, entropy_norm = entropy_norm))
}

