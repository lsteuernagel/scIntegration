# These functions are used to evaluate the results from integrations of sc-seq data in terms of batch mixing and biological purity
# using low dimensional embeddings as input.

##########
### map_celltype_signatures
##########

#' Map gene signatures characterizing celltypes to cells in the merged dataset using AUC. These mapped celltypes can be used to asses biological purity.
#' Requires to define a set of signatures beforehand! Run reduce_celltype_sets afterwards to make a more reliable set of celltypes for evaluation.
#' @param seurat_object Merged seurat object after computation of reductions (PCA)
#' @param subset_cell_ids subset matrix to these cell ids. NULL for all
#' @param aucMaxRank_n how many genes are used for AUC ranking
#' @param gene_set_list set of gene signatures characterizing celltypes
#' @param min_rowSum rowsums per gene when building AUC ranking
#' @param summary_quantile_cutoff quantile from random gene set AUC used as cutoff 
#' @param n_perm how many random gene sets are used to find cutoffs
#' @param alpha linearly scale the threshold obtained by summary_quantile_cutoff with 1+alpha*threshold. Defaults to 0
#' @param global_seed 
#' @return a list of cells per signature.

map_celltype_signatures = function(seurat_object,subset_cell_ids,aucMaxRank_n,gene_set_list,min_rowSum=10,summary_quantile_cutoff = 0.99,
                                   n_perm = 100,alpha=1 ,global_seed,thrP=0.01, smallestPopPercent=0.01){
  
  message(Sys.time(),": Running map_celltype_signatures on seurat object: ",seurat_object@project.name," with ",length(gene_set_list)," signatures" )
  require(tidyverse)
  require(Seurat)
  #require(AUCell)
  # quick fix to this https://github.com/aertslab/AUCell/issues/10
  library("AUCell",lib.loc = "/beegfs/scratch/bruening_scratch/lsteuernagel/R/user_lib/x86_64-pc-linux-gnu-library/4.0/")
  require(Matrix)
  
  ### Build AUC ranking 
  exprMatrix <- seurat_object@assays$RNA@counts
  if(!is.null(subset_cell_ids)){
    idx = which(seurat_object@meta.data[,"Cell_ID"] %in% subset_cell_ids)
    message("Subsetting matrix to ",length(idx)," indices.")
    exprMatrix <- exprMatrix[,idx]
  }
  exprMatrix <- exprMatrix[Matrix::rowSums(exprMatrix)>min_rowSum,]
  # Build gene-expression rankings for each cell
  message("Building AUC ranking")
  set.seed(123)
  cells_rankings <- AUCell::AUCell_buildRankings(exprMatrix, nCores=1, plotStats=F,verbose = F)
  rm(exprMatrix)
  gc()
  
  ### Build permutation based reference for cutoffs
  # message("Find permutation based thresholds using quantile: ",summary_quantile_cutoff," and alpha: ",alpha)
  # setsizes_for_permutation = lapply(gene_set_list,length)
  # result_list_orig=list()
  # cc=0
  # for(i in setsizes_for_permutation){
  #   cc=cc+1
  #   message("Find permutation based thresholds for gene set of length ",i)
  #   genes = rownames(cells_rankings@assays@data@listData$ranking)
  #   gene_set = list()
  #   for( j in 1:n_perm){
  #     set.seed(global_seed+j)
  #     gene_set[[as.character(j)]] = genes[sample(1:length(genes),i)]
  #   }
  #   #Calculate enrichment for gene signatures (AUC)
  #   permutation_AUC <- AUCell_calcAUC(gene_set, cells_rankings,nCores=1,verbose=TRUE,aucMaxRank=aucMaxRank_n)
  #   permutation_AUC_mat = permutation_AUC@assays@data@listData$AUC
  #   permutation_AUC_orig = apply(permutation_AUC_mat,2,quantile,probs=summary_quantile_cutoff)
  #   result_list_orig[[names(gene_set_list)[cc]]] = permutation_AUC_orig
  # }
  # permutation_cutoffs = as.data.frame(do.call(rbind,result_list_orig),stringsAsFactors=F)
  
  ### Calculate AUC for signatures and identify cells
  message("Calculate AUC for signatures and extract mapped cells")
  signature_AUC <- AUCell_calcAUC(gene_set_list, cells_rankings,nCores=1,verbose=F,aucMaxRank=aucMaxRank_n)
  signature_AUC_mat = signature_AUC@assays@data@listData$AUC
  #get cells
  cells_per_gene_set = list()
  for(i in 1:length(gene_set_list)){
    auc_vector = signature_AUC_mat[names(gene_set_list)[i],]
    glProb <- 1-(thrP/length(auc_vector) + smallestPopPercent)  
    Global_k1 = qnorm(glProb,mean=mean(auc_vector),sd=sd(auc_vector))
    cutoff = Global_k1*alpha
    cells_per_gene_set[[names(gene_set_list)[i]]] = names(auc_vector)[auc_vector > cutoff]
    #cells_per_gene_set[[names(gene_set_list)[i]]] = names(auc_vector)[auc_vector > (cutoffs*(1+alpha))]
    message("For gene set: ",names(gene_set_list)[i]," , ",length(cells_per_gene_set[[names(gene_set_list)[i]]])," cells were detected above the treshold, with a threshold: ",cutoff)
  }
  
  return(cells_per_gene_set)
  
} 

##########
### reduce_celltype_sets
##########

#' Reduce mapped celltypes. use after running map_celltype_signatures
#' @param seurat_object Merged seurat object after computation of reductions (PCA)
#' @param cells_per_gene_set cells from seurat_object mapped to signatures
#' @param percentage_unique_min after mapping cells to a signature, how many cells must uniquely map to a signature to keep it.
#' @param minimum_cells after mapping cells to a signature, how many cells must at least be mapped to a signature
#' @param global_seed 
#' @return a list of cells per signature after reduction to relevant ones.

reduce_celltype_sets = function(seurat_object,cells_per_gene_set,percentage_unique_min = 50, minimum_cells = 50,global_seed){
  
  ### Reduce celltypes
  message(Sys.time(),"Reducing to relevant mapped celltypes")
  # filter for minimum cells
  cells_per_gene_set_filtered = cells_per_gene_set[lapply(cells_per_gene_set,length)>minimum_cells]
  
  # go top down through all celltype and throw out sets with too few unique cells
  sorted_names = names(sort(unlist(lapply(cells_per_gene_set_filtered,length)),decreasing = TRUE))
  # keep vector changes dynamically!
  keep_vector= sorted_names
  for(i in 1:length(sorted_names)){
    current_celltype=sorted_names[i]
    current_cells = cells_per_gene_set_filtered[[current_celltype]]
    ## only with not previouly excluded celltypes:
    cells_per_gene_set_excluded = cells_per_gene_set_filtered[keep_vector]
    ## make cell summary
    df_list = list()
    for(i in 1:length(cells_per_gene_set_excluded)){
      df_list[[i]]=data.frame(cell = cells_per_gene_set_excluded[[i]],celltype=names(cells_per_gene_set_excluded)[i],stringsAsFactors = F)
    }
    cells_per_gene_set_long = do.call(rbind,df_list)
    per_cell_summary = cells_per_gene_set_long %>% dplyr::group_by(cell) %>% dplyr::count()
    
    # # count how many cells are unique to celltpe and how many other part of other celltypes as well
    unique_cells = nrow(per_cell_summary[per_cell_summary$cell %in% current_cells & per_cell_summary$n == 1,])
    nonunique_cells = nrow(per_cell_summary[per_cell_summary$cell %in% current_cells & per_cell_summary$n > 1,])
    percent_unique = round(100*unique_cells/(unique_cells+nonunique_cells),3)
    #print(paste0(current_celltype," : ",percent_unique,"% " ))
    # if more than x%  are non-unique --> exclude celltype
    if(percent_unique<percentage_unique_min){
      message("Exclude: ",current_celltype)
      keep_vector = keep_vector[keep_vector != current_celltype]
    }
    
  }
  # subset to celltypes that are kept and return
  cells_per_gene_set_final = cells_per_gene_set_filtered[keep_vector]
  return(cells_per_gene_set_final)
  
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
#' @param ncores how many cores to use by doParallel
#' @param scale_to_max scale non-normalized entropy to maximum number of batches
#' @param max_for_norm all batches with a probability > max_for_norm are considered when normalizing 
#' @param global_seed seed
#' @return list of two vectors: entropy and entropy normalized with log2(expected number of batches)

classProb = function(train_predictors,train_response,trees=500,ncores=4,scale_to_max = FALSE,max_for_norm = 0.01,global_seed = 123){
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
    randomForest(x=train_predictors, y=as.factor(train_response), ntree=ntree)
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

##########
### celltypeProb
##########

#' Train multiple random forests per celltype and obtain oob probabilities of binary classification (yes or no). Then report the entropy per cell as a measure of batch mixing. 
#' Celltypes with unbalanced numbers of mapped cells will have 
#' @param train_predictors cell x features (genes,PCs,etc)
#' @param cells_sets sets of cells from mapped celltypes
#' @param control_size_factor_max how much bigger can the control sample size be (downsampled if exceeds)
#' @param trees Number of trees, should be rather high to obtain accurate oob probabilities
#' @param ncores how many cores to use by doParallel
#' @param global_seed seed
#' @return list of two vectors: entropy and entropy normalized with log2(expected number of batches)

celltypeProb = function(train_predictors,cells_sets,control_size_factor_max = 5,trees=1000,ncores=4,global_seed = 123){
  require(randomForest)
  require(foreach)
  require(doParallel)
  require(doRNG)
  if(any(!unique(as.character(unlist(cells_sets))) %in% rownames(train_predictors))){
    warning("cell_sets should be rownames of train predictors.")
  }
  
  # if(subsample>0){
  #   target_cells = rownames(train_predictors)[rownames(train_predictors) %in% unique(as.character(unlist(cells_sets)))]
  #   other_cells =  rownames(train_predictors)[!rownames(train_predictors) %in% unique(as.character(unlist(cells_sets)))]
  #   set.seed(global_seed)
  #   sampled = base::sample(other_cells,size=max(subsample,length(target_cells)),replace = FALSE)
  #   train_predictors = train_predictors[c(target_cells,sampled),]
  # }
  
  # loop through celltypes
  all_celltypes_celltypeProb_list = list()
  for(i in 1:length(cells_sets)){
    current_celltype = names(cells_sets)[i]
    message(" ",current_celltype)
    # subset
    target_cells = rownames(train_predictors)[rownames(train_predictors) %in% unique(as.character(unlist(cells_sets[[current_celltype]])))]
    other_cells =  rownames(train_predictors)[!rownames(train_predictors) %in% target_cells]
    if(length(other_cells)>(length(target_cells)*control_size_factor_max)){
      message("Downsampling since: ",length(other_cells), " > ",length(target_cells)*control_size_factor_max)
      set.seed(global_seed)
      sampled = base::sample(other_cells,size=length(target_cells)*control_size_factor_max,replace = FALSE)
      train_predictors_subset = train_predictors[c(target_cells,sampled),]
    }else{
      train_predictors_subset=train_predictors
    }
    # make labels 
    labels = rep("no",nrow(train_predictors_subset))
    labels[which(rownames(train_predictors_subset) %in% cells_sets[[current_celltype]])] = "yes"
    message("Celltype samples (no/yes): ",paste0(table(labels),collapse = " ; "))
    # run random forest
    registerDoParallel(cores=ncores)
    doRNG::registerDoRNG(seed = global_seed)
    rf_res <- foreach(ntree=rep(trees/ncores, ncores), .combine=randomForest::combine,.multicombine=TRUE, .packages='randomForest') %dopar% {
      randomForest(x=train_predictors_subset, y=as.factor(labels), ntree=ntree)
    }
    # obtain class probabilities
    votes = rf_res$votes / ncores
    probs_cells = votes[rownames(votes) %in% cells_sets[[current_celltype]],"yes"]
    # add
    all_celltypes_celltypeProb_list[[current_celltype]] = probs_cells
  }
  return(all_celltypes_celltypeProb_list)
}


##########
### evaluate_mixing
##########

#' Wrapper that runs random forest based evaluation of embeddings for mixing. 
#' @param seurat_object_metadata metadata associated with the mapped celltypes and emebdding that will be evaluated
#' @param integration_files names of files that contain embeddings
#' @param integration_path directory where to find integration_files
#' @param evaluation_file filename where to store results
#' @param batch_var name of column in metadata
#' @param subset_cell_ids cell ids of subset that is used for evaluation
#' @param ntrees Number of trees, should be rather high to obtain accurate oob probabilities
#' @param ncores how many cores to use by doParallel
#' @param max_dim max features in random forest
#' @param scale_to_max scale non-normalized entropy to maximum number of batches
#' @param max_for_norm all batches with a probability > max_for_norm are considered when normalizing 
#' @param global_seed seed
#' @return Nothing. Writes files

evaluate_mixing = function(seurat_object_metadata,integration_files,integration_path,evaluation_file,batch_var,subset_cell_ids=NULL,ntrees,ncores,max_dim = 100,scale_to_max,max_for_norm,global_seed){
  
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
      # reduce to max dim -->    # ignore max dim!
      #current_embedding=current_embedding[,1:min(ncol(current_embedding),max_dim)]

      # set cell names --> assumes that the order of cells in metadata is the same as in embedding!
      #rownames(current_embedding) = rownames(seurat_object_metadata)
      message("length batch labels pre",length(batch_labels))
      if(!is.null(subset_cell_ids)){
        idx = which(rownames(current_embedding) %in% subset_cell_ids)
        current_embedding = current_embedding[idx,]
        batch_labels_subset = batch_labels[idx]
      }else{
        batch_labels_subset = batch_labels
      }
      message("check nas ",any(is.na(batch_labels_subset)))
      message("length batch labels post",length(batch_labels_subset))
      message("dim",dim(current_embedding))
      # message("headfeh",head(idx))
      # entropy on class probabilities
      message(Sys.time(),"Evaluating ",current_name)
      classProb_temp= classProb(train_predictors=current_embedding,train_response=batch_labels_subset,trees=ntrees,ncores=ncores,scale_to_max=scale_to_max,max_for_norm = max_for_norm,global_seed =global_seed)
      
      #store results
      classProb_entropy_list[[current_name]] = classProb_temp$entropy
      classProb_entropy_norm_list[[current_name]] = classProb_temp$entropy_nor
    }
  }
  
  # save updated files
  # only attempt saving if results are there
  if(length(classProb_entropy_list)>0){
    # combine results
    classProb_entropy = do.call(cbind,classProb_entropy_list)
    classProb_entropy_norm = do.call(cbind,classProb_entropy_norm_list)
    message("Saving results.")
    utils::write.table(as.data.frame(classProb_entropy),file =  evaluation_file,sep =  "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
    utils::write.table(as.data.frame(classProb_entropy_norm),file =  gsub("evaluation_mixing_prob","evaluation_mixing_probNorm",evaluation_file), sep =  "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
    # readr::write_delim(as.data.frame(classProb_entropy),path = evaluation_file,delim = "\t",col_names = TRUE)
    # readr::write_delim(as.data.frame(classProb_entropy_norm),path = gsub(".txt","_Norm.txt",evaluation_file),delim = "\t",col_names = TRUE)
    # saveRDS(classProb_entropy, evaluation_file)
    # saveRDS(classProb_entropy_norm, gsub(".txt","_Norm.txt",evaluation_file))
  }
}

##########
### evaluate_purity
##########

#' Wrapper that runs random forest based evaluation of embeddings for celltype purity.
#' @param seurat_object_metadata metadata associated with the mapped celltypes and emebdding that will be evaluated
#' @param cells_sets list of vectors (per mapped celltypes) with celllabels matching rownames of metadata
#' @param subset_cell_ids cell ids of subset that is used for evaluation
#' @param control_size_factor_max
#' @param integration_files names of files that contain embeddings
#' @param integration_path directory where to find integration_files
#' @param evaluation_file filename where to store results
#' @param ntrees Number of trees, should be rather high to obtain accurate oob probabilities
#' @param ncores how many cores to use by doParallel
#' @param max_dim maximum number of features in random forest
#' @param scale_to_max scale non-normalized entropy to maximum number of batches
#' @param max_for_norm all batches with a probability > max_for_norm are considered when normalizing 
#' @param global_seed seed
#' @return Nothing. Writes files

evaluate_purity = function(seurat_object_metadata,cells_sets,subset_cell_ids=NULL,control_size_factor_max=5,integration_files,integration_path,evaluation_file,ntrees,ncores,max_dim=50,global_seed){
  
  # check that files can be found at integration path
  
  # make names
  #   split_filepath = sapply(integration_files,function(x){as.character(strsplit(x,split = "/")[[1]])})
  integration_to_run=c()
  for( j in 1:length(integration_files)){
    current_file = integration_files[j]
    split_filepath = as.character(strsplit(current_file,split = "/")[[1]])
    if(length(split_filepath)>1){
      current_file = gsub(".txt","",split_filepath[length(split_filepath)])
    }else{
      current_file = gsub(".txt","",current_file)
    }
    integration_to_run = c(integration_to_run,current_file)
  }
  #init result lists
  celltypeProb_purity_list = list()
  
  # run classProb on all missing results
  for(i in 1:length(integration_to_run)){
    
    # get current result
    current_name = gsub(".txt","",integration_to_run[i])
    current_file = integration_files[which(grepl(paste0(integration_to_run[i],".txt"),integration_files))]
    # message(">>>>>",paste0(integration_path,current_file))
    current_embedding = read_embedding(paste0(integration_path,current_file),seurat_object_metadata = seurat_object_metadata)
    # reduce to max dim
    current_embedding=current_embedding[,1:min(ncol(current_embedding),max_dim)]
    #subsample
    if(!is.null(subset_cell_ids)){
      idx = which(rownames(current_embedding) %in% subset_cell_ids)
      current_embedding = current_embedding[idx,]
    }
    # set cell names --> assumes that the order of cells in metadata is the same as in embedding!
    #rownames(current_embedding) = rownames(seurat_object_metadata)
    
    # prob on class probabilities
    message(Sys.time(),"Evaluating ",current_name)
    # entropy on class probabilities
    celltypeProb_temp= celltypeProb(train_predictors=current_embedding,cells_sets=cells_sets,control_size_factor_max=control_size_factor_max,trees=ntrees,ncores=ncores,global_seed = global_seed)
    #store results
    celltypeProb_purity_list[[current_name]] = celltypeProb_temp
    
  }
  
  # save updated files
  if(length(celltypeProb_purity_list)>0){
    # summarize results to mean_per_cettype x reduction matrix
    names(celltypeProb_purity_list) = integration_to_run
    celltypeProb_purity_summary = lapply(celltypeProb_purity_list,function(x){return(unlist(lapply(x,mean)))})
    celltypeProb_purity_summary = do.call(cbind,celltypeProb_purity_summary)
    
    message("Saving results to files.")
    #saveRDS(purity_prob, paste0(evaluation_path,"evaluation_purity_prob.",length(cells_sets),".",ntrees,".",max_dim,".",global_seed,".txt"))
    #  readr::write_delim(as.data.frame(celltypeProb_purity_summary),path = evaluation_file,delim = "\t",col_names = TRUE)
    utils::write.table(as.data.frame(celltypeProb_purity_summary),file =  evaluation_file,sep =  "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
    # saveRDS(celltypeProb_purity_summary, evaluation_file)
  }
}

##########
### remove_deleted_results
##########

#' Remove integration results that have been deleted from the evaluation file.
#' @param evaluation_file_path
#' @param integration_res_path

remove_deleted_results = function(evaluation_file_path,integration_res_path){
  
  message("Running remove_deleted_results.")
  #mixing_entropyTR = readRDS
  evaluation_file = readRDS(evaluation_file_path)
  all_evaluations = colnames(evaluation_file)
  available_integrations= gsub(".txt","",list.files(integration_res_path,recursive = TRUE,pattern = ".txt"))
  available_integrations = as.character(sapply(available_integrations,function(x){ strsplit(x,split = "/")[[1]][length(strsplit(x,split = "/")[[1]])]}))
  keep = all_evaluations[all_evaluations %in% available_integrations]
  
  message("Removing ",length(all_evaluations)-length(keep)," of ",length(all_evaluations)," evaluation results..")
  
  out = evaluation_file[,keep]
  
  saveRDS(out, evaluation_file_path)
  
  
}

##########
### check_evaluation_results
##########

#' Which integration results are already evaluated
#' @param evaluation_file_path where to look for full evaluation file
#' @param integration_res_path where to look for integration files

check_evaluation_results = function(evaluation_file_path,integration_res_path){
  
  message("Running check_evaluation_results.")
  # read all files with integration results
  all_files=list.files(integration_res_path,recursive = TRUE,pattern = ".txt")
  
  # if a global evaluation file exists:
  if(file.exists(evaluation_file_path)){
    message("Found file, removing integration results that are already present...")
    # read all evaluation results and extract colnames
    evaluation_file = data.table::fread(evaluation_file_path,data.table = F)#readRDS(evaluation_file_path)
    if(grepl("asw",evaluation_file_path)){
      message("Detected asw keyword")
      all_evaluations = as.character(evaluation_file[,1])
    }else if(grepl("separation",evaluation_file_path)){
      message("Detected separation keyword")
      all_evaluations = as.character(evaluation_file[,"embedding"])
    }else{
      message("Using colnames")
      all_evaluations = colnames(evaluation_file)
    }
    # reformat all_files to be comparable
    available_integrations= gsub(".txt","",all_files)
    available_integrations = as.character(sapply(available_integrations,function(x){ strsplit(x,split = "/")[[1]][length(strsplit(x,split = "/")[[1]])]}))
    # which files are NOT in global evaluation file
    not_yet_evaluated = all_files[!available_integrations %in% all_evaluations]
  }else{
    # evaluate all files
    not_yet_evaluated=all_files
  }
  return(not_yet_evaluated)
  
}

##########
### merge_evaluation_results
##########

#' Merge individual evaluation results into one file. With the current pipeline reset_file should be TRUE.
#' @param evaluation_file_path where to look for files that should be merged
#' @param type_string what type of evaluation result: 'evaluation_mixing_prob', 'evaluation_mixing_probNorm', evaluation_purity_prob
#' @param reset_file overwrite existing file
#' @return the filename with all results

merge_evaluation_results = function(evaluation_file_path,type_string,expected_ntrees="",expected_seed="",expected_clusters="",expected_n_celltypes="",reset_file=TRUE){
  
  message("Running merge_evaluation_results")
  message("Warning: This function expect a defined order in the filenames from the evaluation pipeline, which makes it unflexible when changing file names in the mains scripts!")
  
  ### if this function is not working check that expected_ntrees and expected_seed are included in the call!!!! ####
  
  # find all relevant files
  all_eval_files = list.files(evaluation_file_path,pattern=paste0(type_string,"\\."))
  # exclude summary files
  all_eval_files = all_eval_files[!grepl("_all.txt",all_eval_files)]
  # print(length(all_eval_files))
  # add together
  message("Found ",length(all_eval_files)," files. Load and add together")
  tmp_list=list()
  cc = 1
  for(i in 1:length(all_eval_files)){
    name_details = strsplit(all_eval_files[i],split = "\\.")[[1]]
    if(type_string %in% c("evaluation_mixing_prob","evaluation_purity_prob","evaluation_mixing_probNorm","evaluation_purity_rmse","evaluation_purity_rmse_norm") & name_details[4] == expected_ntrees & name_details[6] == expected_seed){
      name_with_file = all_eval_files[i]
      tmp_list[[cc]] = data.table::fread(paste0(evaluation_file_path,all_eval_files[i]),data.table = F)
      cc = cc+1
    }else if(type_string %in% c("evaluation_asw","evaluation_entropy_knn_all","evaluation_separation")){
      name_with_file = all_eval_files[i]
      tmp_list[[cc]] = data.table::fread(paste0(evaluation_file_path,all_eval_files[i]),data.table = F)
      cc = cc+1
    }else if(type_string %in% c("evaluation_purity_knn","evaluation_entropy_knn") & name_details[3] == expected_n_celltypes ){
      name_with_file = all_eval_files[i]
      tmp_list[[cc]] = data.table::fread(paste0(evaluation_file_path,all_eval_files[i]),data.table = F)
      cc = cc+1
    }else if(type_string %in% c("purity_perBatch_cluster") & name_details[3] == expected_clusters ){
      name_with_file = all_eval_files[i]
      tmp_list[[cc]] = data.table::fread(paste0(evaluation_file_path,all_eval_files[i]),data.table = F)
      cc = cc+1
    }
    # tmp_list[[i]] = utils::read.table(paste0(evaluation_file_path,all_eval_files[i]),col.names = TRUE,row.names = FALSE)
  }
  # print(length(tmp_list))
  
  message("Define file name")
  name_details = strsplit(name_with_file,split = "\\.")[[1]]
  if(grepl("asw|separation",type_string)){
    # when using silhouette score files
    tmp_mat = do.call(rbind,tmp_list)
    full_file_name = paste0(evaluation_file_path,type_string,".",paste0(name_details[2:4],collapse = "."),"_all.txt")
  }else if(grepl("purity_knn",type_string)){
    tmp_mat = do.call(cbind,tmp_list)
    full_file_name = paste0(evaluation_file_path,type_string,".",paste0(name_details[3:5],collapse = "."),"_all.txt")
  }else if(grepl("entropy_knn",type_string)){
    tmp_mat = do.call(cbind,tmp_list)
    full_file_name = paste0(evaluation_file_path,type_string,".",paste0(name_details[3:5],collapse = "."),"_all.txt")
  }else if(grepl("evaluation_entropy_knn_all",type_string)){
    tmp_mat = do.call(cbind,tmp_list)
    full_file_name = paste0(evaluation_file_path,type_string,".",paste0(name_details[3:5],collapse = "."),"_all.txt")
  }else if(grepl("purity_perBatch_cluster",type_string)){
    tmp_mat = do.call(cbind,tmp_list)
    full_file_name = paste0(evaluation_file_path,type_string,".",paste0(name_details[3:4],collapse = "."),"_all.txt")
  }else{
    # when using other file
    tmp_mat = do.call(cbind,tmp_list)
    full_file_name = paste0(evaluation_file_path,type_string,".",paste0(name_details[3:6],collapse = "."),"_all.txt")
  }
  #  print(dim(tmp_mat))
  # check if file exists, filename contains length(cells_sets) ntrees and the seed used, so if this changes a new file will be created
  message("Saving merged file")
  if(!file.exists(full_file_name)|reset_file){
    message("Saving all results in new file: ",full_file_name)
    utils::write.table(as.data.frame(tmp_mat),file = full_file_name,sep =  "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
  }else{
    previous_results = data.table::fread(full_file_name,data.table = F)
    if(grepl("asw",type_string)){
      # when using silhouette score files
      if(ncol(previous_results)!=ncol(tmp_mat)){
        # throw error here?
        warning("Different number of cols. Not adding new results")
        updated_results=previous_results
      }else{
        # else cbind
        print("rg")
        updated_results = rbind(previous_results,tmp_mat)
      }
    }else{
      # when using other files
      if(nrow(previous_results)!=nrow(tmp_mat)){
        # throw error here?
        warning("Different number of rows. Not adding new results")
        updated_results=previous_results
      }else{
        # else cbind
        updated_results = cbind(previous_results,tmp_mat)
      }
    }
    message("Adding results to new file: ",full_file_name)
    utils::write.table(as.data.frame(updated_results),file = full_file_name,sep =  "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
  }
  return(full_file_name)
  
}


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
  #require(AUCell)
  # quick fix to this https://github.com/aertslab/AUCell/issues/10
  library("AUCell",lib.loc = "/beegfs/scratch/bruening_scratch/lsteuernagel/R/user_lib/x86_64-pc-linux-gnu-library/4.0/")
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
    signature_AUC <- AUCell_calcAUC(gene_set_list, cells_rankings,nCores=1,verbose=F,aucMaxRank=aucMaxRank_n)
    signature_AUC_mat = t(signature_AUC@assays@data@listData$AUC)
    all_signature_auc_mats[[i]] = signature_AUC_mat
  }
  signature_AUC_mat_full = do.call(rbind,all_signature_auc_mats)
  #gc()
  return(signature_AUC_mat_full)
  
} 

##########
### evaluate_purity_regression
##########

#' Regression based random forest purity evaluation. Needs output from map_celltype_signatures2
#' @param reduction 
#' @param signature_AUC_mat rows(cells) x cols(celltypes from AUC)
#' @param min_thresh 
#' @param ntrees 
#' @param ncores 
#' @param global_seed 
#' @return 

# set.seed(123)
# exprMatrix = seurat_object@assays$RNA@counts[,sample(1:ncol(seurat_object@assays$RNA@counts),10240)]
# reduction=seurat_object_scVI@reductions$scVI_7_200_0.1_2_128_30..scVI..features_scVI_RNA.log.vst.split_Batch_ID.features.2500_87cb9293d14869b650fabf1c07c03b5d@cell.embeddings[sample(1:ncol(seurat_object@assays$RNA@counts),10240),]

evaluate_purity_regression = function(reduction,signature_AUC_mat,min_thresh=NULL,ntrees=1000,min_total_samples=1000,alpha=1,thrP=0.01,smallestPopPercent=0.01,ncores=1,increase_trees = TRUE,global_seed =123){
  
  require(randomForest)
  require(foreach)
  require(doParallel)
  require(doRNG)
  
  all_celltypes = colnames(signature_AUC_mat)
  all_celltypes_results = c()
  all_celltypes_results_norm =c()
  ### For each celltype:
  for(i in 1:length(all_celltypes)){
    celltype_values = signature_AUC_mat[,i]
    
    # get threshold as
    min_thresh=NULL
    if(is.null(min_thresh)){
      glProb <- 1-(thrP/length(celltype_values) + smallestPopPercent)
      min_thresh = qnorm(glProb,mean=mean(celltype_values),sd=sd(celltype_values)) * alpha
    }
    # get all values above min_thresh
    keep_idx = which(celltype_values>min_thresh)
    # if min_total_samples < positives * 2
    if(min_total_samples > length(keep_idx)*2){
      negatives_size = min_total_samples - length(keep_idx)
    }else{ # negatives = same amount with < Global_k1
      negatives_size = length(keep_idx)
    }
    # sample negatives
    not_keep_idx = which(celltype_values <= min_thresh)
    set.seed(global_seed)
    message("Sampling: ",negatives_size)
    additional_idx = sample(not_keep_idx,negatives_size)
    all_idx = c(keep_idx,additional_idx)
    # subset to train data and
    train_var = reduction[all_idx,]
    predictor_var = celltype_values[all_idx]
    # determine trees:
    ntrees_use=ntrees
    if(increase_trees & length(predictor_var)>ntrees){
      ntrees_use = length(predictor_var)
    }
    if(increase_trees & ntrees > (length(predictor_var)*5)){
      ntrees_use = length(predictor_var)*5
    }
    message("Running regression random forest for: ",all_celltypes[i]," using ",length(predictor_var)," observations and ",ntrees_use," trees. Based on min_tresh: ",min_thresh)
    # run regression random forest
    registerDoParallel(cores=ncores)
    doRNG::registerDoRNG(seed = global_seed)
    rf_res <- foreach(ntree=rep(ntrees_use/ncores, ncores), .combine=randomForest::combine,.multicombine=TRUE, .packages='randomForest') %dopar% {
      randomForest(x=train_var, y=predictor_var, ntree=ntree)
    }
    # save training RMSE
    all_celltypes_results[all_celltypes[i]]=sqrt(mean((rf_res$y-rf_res$predicted)^2,na.rm=TRUE))
    message("RMSE of: ",all_celltypes[i]," : ",sqrt(mean((rf_res$y-rf_res$predicted)^2,na.rm=TRUE)))
    all_celltypes_results_norm[all_celltypes[i]]= sqrt(mean((rf_res$y)^2,na.rm=TRUE)) / sqrt(mean((rf_res$y-rf_res$predicted)^2,na.rm=TRUE))
    message("RMSE of: ",all_celltypes[i]," : ",all_celltypes_results_norm[all_celltypes[i]])
  }
  return(list(all_celltypes_results = all_celltypes_results, all_celltypes_results_norm = all_celltypes_results_norm))
} 


##########
### evaluate_purity2
##########

#' Wrapper that runs random forest based evaluation of embeddings for celltype purity.
#' @param seurat_object_metadata metadata associated with the mapped celltypes and emebdding that will be evaluated
#' @param cells_sets list of vectors (per mapped celltypes) with celllabels matching rownames of metadata
#' @param subset_cell_ids cell ids of subset that is used for evaluation
#' @param control_size_factor_max
#' @param integration_files names of files that contain embeddings
#' @param integration_path directory where to find integration_files
#' @param evaluation_file filename where to store results
#' @param ntrees Number of trees, should be rather high to obtain accurate oob probabilities
#' @param ncores how many cores to use by doParallel
#' @param max_dim maximum number of features in random forest
#' @param scale_to_max scale non-normalized entropy to maximum number of batches
#' @param max_for_norm all batches with a probability > max_for_norm are considered when normalizing 
#' @param global_seed seed
#' @return Nothing. Writes files

evaluate_purity2 = function(seurat_object_metadata,signature_AUC_mat,control_size_factor_max=5,min_thresh=0.05,min_total_samples=1000,alpha=1,integration_files,integration_path,evaluation_file,ntrees,ncores,max_dim=50,global_seed){
  
  # check that files can be found at integration path
  
  # make names
  #   split_filepath = sapply(integration_files,function(x){as.character(strsplit(x,split = "/")[[1]])})
  integration_to_run=c()
  for( j in 1:length(integration_files)){
    current_file = integration_files[j]
    split_filepath = as.character(strsplit(current_file,split = "/")[[1]])
    if(length(split_filepath)>1){
      current_file = gsub(".txt","",split_filepath[length(split_filepath)])
    }else{
      current_file = gsub(".txt","",current_file)
    }
    integration_to_run = c(integration_to_run,current_file)
  }
  #init result lists
  celltypeProb_purity_list = list()
  #celltypeProb_purity_norm_list  = list()
  # run classProb on all missing results
  for(i in 1:length(integration_to_run)){
    
    # get current result
    current_name = gsub(".txt","",integration_to_run[i])
    current_file = integration_files[which(grepl(paste0(integration_to_run[i],".txt"),integration_files))]
    # message(">>>>>",paste0(integration_path,current_file))
    current_embedding = read_embedding(paste0(integration_path,current_file),seurat_object_metadata = seurat_object_metadata)
    # reduce to max dim
    current_embedding=current_embedding[,1:min(ncol(current_embedding),max_dim)]
    # set cell names --> assumes that the order of cells in metadata is the same as in embedding!
    #rownames(current_embedding) = rownames(seurat_object_metadata)
    
    # prob on class probabilities
    message(Sys.time(),"Evaluating ",current_name)
    # entropy on class probabilities
    # celltypeProb_temp= celltypeProb(train_predictors=current_embedding,cells_sets=cells_sets,control_size_factor_max=control_size_factor_max,trees=ntrees,ncores=ncores,global_seed = global_seed)
    
    #regression:
    # celltypeRMSE_list = evaluate_purity_regression(reduction=current_embedding,signature_AUC_mat,min_thresh=min_thresh,ntrees=ntrees,min_total_samples=min_total_samples,alpha=alpha,ncores=ncores,global_seed =global_seed)
    # #store results
    # celltypeProb_purity_list[[current_name]] = celltypeRMSE_list$all_celltypes_results
    # celltypeProb_purity_norm_list[[current_name]] = celltypeRMSE_list$all_celltypes_results_norm
    
    #
    #celltypeProb2 = function(train_predictors,signature_AUC_mat,trees=1000,min_total_samples=2000,alpha=1,min_pos_thresh=0.05,thrP=0.01,smallestPopPercent=0.01,ncores=4,global_seed = 123){
    
    celltypeProb_temp = celltypeProb2(train_predictors=current_embedding,signature_AUC_mat,min_total_samples=min_total_samples,alpha=alpha,min_pos_thresh=min_thresh,
                                      thrP=0.01,smallestPopPercent=0.01,trees=ntrees,ncores=ncores,global_seed =global_seed)
    #store results
    celltypeProb_purity_list[[current_name]] = celltypeProb_temp
    
  }
  
  # save updated files
  if(length(celltypeProb_purity_list)>0){
    # summarize results to mean_per_cettype x reduction matrix
    names(celltypeProb_purity_list) = integration_to_run
    celltypeProb_purity_summary = lapply(celltypeProb_purity_list,function(x){return(unlist(lapply(x,mean)))})
    celltypeProb_purity_summary = do.call(cbind,celltypeProb_purity_summary)
    message("Saving results to files.")
    data.table::fwrite(as.data.frame(celltypeProb_purity_summary),file =  evaluation_file,sep =  "\t")
    #evaluation_file_norm =  gsub("rmse","rmse_norm",evaluation_file) # HARDCODED!
    #data.table::fwrite(as.data.frame(celltypeProb_purity_norm_summary),file =  evaluation_file_norm,sep =  "\t")
  }
}

##########
### celltypeProb2
##########

#' Wrapper that runs random forest based evaluation of embeddings for celltype purity.
#' @param train_predictors reduction
#' @param signature_AUC_mat
#' 
#' 
#' @param global_seed seed
#' @return all_celltypes_celltypeProb_list
#' 
celltypeProb2 = function(train_predictors,signature_AUC_mat,trees=1000,min_total_samples=2000,alpha=1,min_pos_thresh=0.05,thrP=0.01,smallestPopPercent=0.01,ncores=4,global_seed = 123){
  require(randomForest)
  require(foreach)
  require(doParallel)
  require(doRNG)
  
  all_celltypes = colnames(signature_AUC_mat)
  # loop through celltypes
  all_celltypes_celltypeProb_list = list()
  for(i in 1:length(all_celltypes)){
    celltype_values = signature_AUC_mat[,i]
    #print(hist(celltype_values))
    current_celltype = all_celltypes[i]
    message(" ",current_celltype)
    # subset
    # get threshold 
    if(alpha<1){alpha=1}
    glProb <- 1-(thrP/length(celltype_values) + smallestPopPercent)
    pos_thresh = qnorm(glProb,mean=mean(celltype_values),sd=sd(celltype_values)) * alpha
    if(min_pos_thresh>pos_thresh){pos_thresh = min_pos_thresh}
    neg_thresh = pos_thresh / alpha
    message("  Neg_tresh: ",neg_thresh," pos_thresh: ",pos_thresh)
    # get all values above min_thresh
    pos_idx = which(celltype_values>pos_thresh)
    neg_idx = which(celltype_values<neg_thresh)
    # downsample negatives ?
    if(min_total_samples > length(pos_idx)*2){
      negatives_size = min_total_samples - length(pos_idx)
    }else{ # negatives = same amount with < Global_k1
      negatives_size = length(pos_idx)
    }
    if(length(neg_idx)>negatives_size){
      set.seed(global_seed)
      neg_idx = sample(neg_idx,negatives_size)
    }
    # make labels 
    labels = c(rep("yes",length(pos_idx)),rep("no",length(neg_idx)))
    message("  Celltype samples (no/yes): ",paste0(table(labels),collapse = " ; "))
    all_idx = c(pos_idx,neg_idx)
    # subset to train data and
    train_var = train_predictors[all_idx,]
    # run random forest
    registerDoParallel(cores=ncores)
    doRNG::registerDoRNG(seed = global_seed)
    rf_res <- foreach(ntree=rep(trees/ncores, ncores), .combine=randomForest::combine,.multicombine=TRUE, .packages='randomForest') %dopar% {
      randomForest(x=train_var, y=as.factor(labels), ntree=ntree)
    }
    # obtain class probabilities
    votes = rf_res$votes / ncores
    probs_cells = votes[1:length(pos_idx),"yes"]
    message("  Mean prob of positives: ",mean(probs_cells))
    # add
    all_celltypes_celltypeProb_list[[current_celltype]] = probs_cells
  }
  return(all_celltypes_celltypeProb_list)
}

##########
### evaluate_purity_knn
##########

#' Wrapper that runs random forest based evaluation of embeddings for celltype purity.
#' @param seurat_object_metadata metadata associated with the mapped celltypes and emebdding that will be evaluated
#' @param cells_sets list of vectors (per mapped celltypes) with celllabels matching rownames of metadata
#' @param k_param cell ids of subset that is used for evaluation
#' @param control_size_factor_max
#' @param integration_files names of files that contain embeddings
#' @param integration_path directory where to find integration_files
#' @param evaluation_file filename where to store results
#' @param ncores how many cores to use by doParallel
#' @param max_dim maximum number of features in random forest
#' @param dist_type all batches with a probability > max_for_norm are considered when normalizing 
#' @param global_seed seed
#' @return Nothing. Writes files

evaluate_purity_knn = function(seurat_object_metadata,cells_sets,k_param=20,integration_files,integration_path,evaluation_file,ncores=1,dist_type="cosine",max_dim=50,global_seed=123){
  
  require(randomForest)
  require(foreach)
  require(doParallel)
  require(doRNG)
  
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
    # reduce to max dim --> no!
   # current_embedding=current_embedding[,1:min(ncol(current_embedding),max_dim)]
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
    
    all_occ_scores = c()
    for(j in 1:length(cells_sets)){
      # prob on class probabilities
      #message(Sys.time(),"Evaluating ",names(cells_sets)[j])
      #
      pos_idx = which(seurat_object_metadata$Cell_ID %in% cells_sets[[j]])
      ## calculate median over neighbors for all cells
      pos_nn_idx = nn_idx[pos_idx,2:ncol(nn_idx)]
      #pos_nn_scores = apply(pos_nn_idx,2,function(indices,values,thresh){y=values[indices];y[y>=thresh]=1;y[y<1]=0;return(y)},values=celltype_values,thresh=pos_thresh)
      pos_nn_scores = apply(pos_nn_idx,2,function(indices,pos_indices){y=rep(0,length(indices));y[indices %in% pos_indices]=1;return(y)},pos_indices=pos_idx)
      pos_nn_occ = apply(pos_nn_scores,1,sum) / ncol(pos_nn_scores)
      all_occ_scores = c(all_occ_scores,mean(pos_nn_occ))
    }
    message("Calculated knn purity")
    #store
    names(all_occ_scores) = names(cells_sets)
    # celltypeProb_purity_list[[current_name]] = celltypeProb_temp
    # celltypeProb_purity_list
    all_occ_scores
  }
  message(paste0(dim(res_par_1)))
  colnames(res_par_1) = integration_to_run
  # save updated files
  message("Saving results to files.")
  data.table::fwrite(as.data.frame(res_par_1),file =  evaluation_file,sep =  "\t")
  # if(length(celltypeProb_purity_list)>0){
  #   # summarize results to mean_per_cettype x reduction matrix
  #   #names(celltypeProb_purity_list) = integration_to_run
  #   #celltypeProb_purity_summary = lapply(celltypeProb_purity_list,function(x){return(unlist(lapply(x,mean)))})
  #   celltypeProb_purity_summary = do.call(cbind,celltypeProb_purity_list)
  #   message("Saving results to files.")
  #   data.table::fwrite(as.data.frame(celltypeProb_purity_summary),file =  evaluation_file,sep =  "\t")
  #   #evaluation_file_norm =  gsub("rmse","rmse_norm",evaluation_file) # HARDCODED!
  #   #data.table::fwrite(as.data.frame(celltypeProb_purity_norm_summary),file =  evaluation_file_norm,sep =  "\t")
  # }
}

##########
### evaluate_purity_perBatch
##########

#' Wrapper that runs random forest based evaluation of embeddings for celltype purity.
#' @param seurat_object_metadata metadata associated with the mapped celltypes and emebdding that will be evaluated
#' @param cells_sets list of vectors (per mapped celltypes) with celllabels matching rownames of metadata
#' @param k_param cell ids of subset that is used for evaluation
#' @param control_size_factor_max
#' @param integration_files names of files that contain embeddings
#' @param integration_path directory where to find integration_files
#' @param evaluation_file filename where to store results
#' @param ncores how many cores to use by doParallel
#' @param max_dim maximum number of features in random forest
#' @param dist_type all batches with a probability > max_for_norm are considered when normalizing 
#' @param global_seed seed
#' @return Nothing. Writes files

evaluate_purity_perBatch = function(seurat_object_metadata,cell_sets,integration_files,integration_path,evaluation_file,ncores=1,max_dim=100,global_seed=123){
  
  require(randomForest)
  require(foreach)
  require(doParallel)
  require(doRNG)
  
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
  registerDoParallel(cores=ncores)
  doRNG::registerDoRNG(seed = global_seed)
  res_par_1 <- foreach(embed_idx = 1:length(integration_to_run), .combine='cbind') %dopar% {  
    # get current result
    current_name = gsub(".txt","",integration_to_run[embed_idx])
    message(current_name)
    current_file = integration_files[which(grepl(paste0(integration_to_run[embed_idx],".txt"),integration_files))][1]
    # message(">>>>>",paste0(integration_path,current_file))
    current_embedding = read_embedding(paste0(integration_path,current_file),seurat_object_metadata = seurat_object_metadata)
    # reduce to max dim
    current_embedding=current_embedding[,1:min(ncol(current_embedding),max_dim)]
    current_embedding = as.matrix(current_embedding)
    # parallel
    cell_sets_results =list()
    for(i in 1:length(cell_sets)){
    #for(embed_idx in 1:2){
    #celltype_batch_sets_all <- foreach(i = 1:length(cell_sets), .combine='rbind') %dopar% {  
      if(i%%100==0)message(i)
      current_cells = cell_sets[[i]]
      current_name = names(cell_sets)[i]
      mat_cells = as.matrix(current_embedding[current_cells,])
      dist_mat = fast_cosine(mat_cells)
      dist_mat_lower = sort(dist_mat[lower.tri(dist_mat)])
      median_dist = median(dist_mat_lower)
      ratio_dist =median(dist_mat_lower[(length(dist_mat_lower)-floor(length(dist_mat_lower)*0.1)):length(dist_mat_lower)]) / median(dist_mat_lower[1:floor(length(dist_mat_lower)*0.1)])
      all_celltypes_dist_summary = c(cell_set = current_name,median_dist=median_dist,ratio_dist=ratio_dist,ncells=length(current_cells))
      cell_sets_results[[current_name]] = all_celltypes_dist_summary
    }
    cell_sets_results_df = do.call(rbind,cell_sets_results)
    message(paste0(head(cell_sets_results_df),collapse = " | "))
    res_vec = cell_sets_results_df[,2]
    names(res_vec) = cell_sets_results_df[,1]
    res_vec
  }
  # save updated files
  colnames(res_par_1) = integration_to_run
  message("Saving results to files.")
  data.table::fwrite(as.data.frame(res_par_1),file =  evaluation_file,sep =  "\t")
}

## cosine similarity
#Convert to cosine dissimilarity matrix (distance matrix).
#https://stats.stackexchange.com/questions/31565/compute-a-cosine-dissimilarity-matrix-in-r
fast_cosine = function(mat){
  sim <- mat / sqrt(rowSums(mat * mat)) # normalize by sum of each row (vector)
  sim <- sim %*% t(sim) # 
  #dissim <- 1 - sim
  return(sim)
}

##########
### evaluate_knn_entropy
##########

#' Wrapper that runs random forest based evaluation of embeddings for celltype purity + mixing.
#' @param seurat_object_metadata metadata associated with the mapped celltypes and emebdding that will be evaluated
#' @param cells_sets list of vectors (per mapped celltypes) with celllabels matching rownames of metadata
#' @param k_param cell ids of subset that is used for evaluation
#' @param control_size_factor_max
#' @param integration_files names of files that contain embeddings
#' @param integration_path directory where to find integration_files
#' @param evaluation_file filename where to store results
#' @param ncores how many cores to use by doParallel
#' @param max_dim maximum number of features in random forest
#' @param dist_type all batches with a probability > max_for_norm are considered when normalizing 
#' @param global_seed seed
#' @return Nothing. Writes files

evaluate_knn_entropy = function(seurat_object_metadata,batch_var="Batch_ID",cells_sets,k_param=20,integration_files,integration_path,evaluation_file,ncores=1,dist_type="cosine",max_dim=50,global_seed=123){
  
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
    # reduce to max dim --> no!
    # current_embedding=current_embedding[,1:min(ncol(current_embedding),max_dim)]
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
    
    all_entropy_scores = c()
    for(j in 1:length(cells_sets)){
      pos_idx = which(seurat_object_metadata$Cell_ID %in% cells_sets[[j]]) # get current cells
      pos_nn_idx = nn_idx[pos_idx,2:ncol(nn_idx)] # get idx for current cells
      entropy = apply(pos_nn_idx,1,function(row,cell_labels){ # run entropy on batch distribution for celltype cells
        freq_batch = table(cell_labels[row])/length(cell_labels[row])
        freq_batch = freq_batch[freq_batch > 0]
        entropy = entropy_fun(freq_batch,logfun ="log2")
        entropy
      }, cell_labels = cell_labels)
      batch_ident=as.character(seurat_object_metadata[pos_idx,batch_var]) # get number of expected abtches for celltype
      entropy = entropy / max(1,log2(length(unique(batch_ident)))) # normalize by this number so that max entropy =1
      entropy[entropy>1]=1 # if there are more batches than expected in neighborhood entropy could get larger than 1 --> set to one because this is overcorrecting!
      all_entropy_scores = c(all_entropy_scores,median(entropy)) # return mean for celltype
    }
    message("Calculated knn mixing entropy")
    #store
    names(all_entropy_scores) = names(cells_sets)
    # celltypeProb_purity_list[[current_name]] = celltypeProb_temp
    # celltypeProb_purity_list
    all_entropy_scores
  }
  message(paste0(dim(res_par_1)))
  colnames(res_par_1) = integration_to_run
  # save updated files
  message("Saving results to files.")
  data.table::fwrite(as.data.frame(res_par_1),file =  evaluation_file,sep =  "\t")

}

##########
### evaluate_knn_entropy_all
##########

#' Wrapper that runs random forest based evaluation of embeddings for celltype purity + mixing.
#' @param seurat_object_metadata metadata associated with the mapped celltypes and emebdding that will be evaluated
#' @param cells_sets list of vectors (per mapped celltypes) with celllabels matching rownames of metadata
#' @param k_param cell ids of subset that is used for evaluation
#' @param control_size_factor_max
#' @param integration_files names of files that contain embeddings
#' @param integration_path directory where to find integration_files
#' @param evaluation_file filename where to store results
#' @param ncores how many cores to use by doParallel
#' @param max_dim maximum number of features in random forest
#' @param dist_type all batches with a probability > max_for_norm are considered when normalizing 
#' @param global_seed seed
#' @return Nothing. Writes files

evaluate_knn_entropy_all = function(seurat_object_metadata,batch_var="Batch_ID",k_param=20,integration_files,integration_path,evaluation_file,ncores=1,dist_type="cosine",max_dim=50,global_seed=123){
  
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
    # reduce to max dim --> no!
    # current_embedding=current_embedding[,1:min(ncol(current_embedding),max_dim)]
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



# from: https://github.com/immunogenomics/LISI/blob/master/R/utils.R 

# or : https://github.com/almutlue/CellMixS/blob/master/R/otherMetrics.R

# lisi_df <- data.frame(matrix(NA, N, length(label_colnames)))
# lisi_df <- Reduce(cbind, lapply(label_colnames, function(label_colname) {
#   labels <- data.frame(meta_data)[, label_colname, drop = TRUE]
#   if (any(is.na(labels))) {
#     message('Cannot compute LISI on missing values')      
#     return(rep(NA, N))
#   } else {
#     ## don't count yourself in your neighborhood
#     dknn$nn.idx <- dknn$nn.idx[, 2:ncol(dknn$nn.idx)]
#     dknn$nn.dists <- dknn$nn.dists[, 2:ncol(dknn$nn.dists)]
#     labels <- as.integer(factor(labels)) - 1
#     n_batches <- length(unique(labels))
#     simpson <- compute_simpson_index(
#       t(dknn$nn.dists),
#       t(dknn$nn.idx) - 1, 
#       labels,
#       n_batches,
#       perplexity
#     )
#     return(1 / simpson)
#   }
# }))
# lisi_df <- as.data.frame(lisi_df)  
# colnames(lisi_df) <- label_colnames
# row.names(lisi_df) <- row.names(meta_data)

# compute_simpson_index
