# this file contains functions that are called to run Harmony with different inputs and parameters

##########
### run_harmony
##########

# run harmony
# call per feature set and PCA to use
# make kmeans init with global seed
# run harmony parameters
# save individual runs

#' Run harmony
#' @param seurat_object
#' @param reduction_name
#' @param param_df
#' @param dims_to_use
#' @param kmeans_nstart
#' @param nclust_K
#' @param kmeans_celltypeMin
#' @param filepath
#' @param global_seed seed
#' @param random_id
#' @return .

run_harmony = function(seurat_object,reduction_name,param_df,dims_to_use=50,kmeans_nstart,nclust_K=NULL,kmeans_celltypeMin=20,filepath,global_seed=123,random_id=NULL){
  
  require(Seurat)
  require(harmony)
  require(Matrix)
  require(data.table)
  require(readr)
 # require(future)
  
  ## set path
  message(Sys.time()," : Running run_harmony on seurat object: ",seurat_object@project.name," and ",reduction_name," reductions." )
  message("Saving export to: ",filepath)
  system(paste0("mkdir -p ",paste0(filepath)))

  #
  if(is.null(random_id)){
    random_id = round(runif(1,min = 0,max = 10000000))
  }
  ### run kmeans as init
  ## Heuristic based number of clusters: init based on counting author clusters
  # nclust_K=100
  if(is.null(nclust_K)){
    celltype_summary = seurat_object@meta.data %>% dplyr::group_by(Dataset,Tissue,Author_CellType) %>% dplyr::summarise(count=n()) %>% 
      dplyr::filter(count>=kmeans_celltypeMin) %>% group_by(Dataset,Tissue) %>% dplyr::summarise(count_celltypes=n()) %>% dplyr::group_by(Tissue) %>%
      dplyr::filter(count_celltypes==max(count_celltypes))
    nclust_K= sum(celltype_summary$count_celltypes)
  }
  ## run kmeans
  message(Sys.time(),": Running kmeans to find common initialization for all Harmony runs (k clusters: ",nclust_K,", nstarts: ",kmeans_nstart,")")
  # input matrix
  max_dim = min(dims_to_use,ncol(seurat_object@reductions[[reduction_name]]@cell.embeddings))
  input_matrix = seurat_object@reductions[[reduction_name]]@cell.embeddings[,1:max_dim]
  #run kmeans
  set.seed(global_seed)
  kmeans_res= kmeans(x=input_matrix,centers = nclust_K, iter.max = 25, nstart = kmeans_nstart) # iter max as in harmony
  ## transform to input for Harmony
  #cluster prior: nrow: clusters, ncol: cells --> try binary matrix with one entry per column?
  vec = as.character(kmeans_res$cluster)
  names(vec) = names(kmeans_res$cluster)
  make_matrix <- function(vec) {
    #https://stackoverflow.com/questions/25271353/cleaner-way-of-constructing-binary-matrix-from-vector  
    U <- sort(unique(vec))
    M <- matrix(0, nrow = length(vec), 
                ncol = length(U), 
                dimnames = list(names(vec), U))
    M[cbind(seq_len(length(vec)), match(vec, U))] <- 1L
    M
  }
  cluster_prior_matrix = make_matrix(vec)
  cluster_prior_matrix =t(cluster_prior_matrix)
  message(Sys.time()," : kmeans initialization finished")
  
  ### call harmony on params
  for(i in 1:nrow(param_df)){
    
    message(Sys.time()," : Running parameter combination ",i," of ",nrow(param_df))
    # get param set into list
    arguments_as_list = as.list(param_df[i, ! colnames(param_df) %in% c("full_id")]) # drop id for function call
    arguments_as_list$lambda = unlist(as.data.frame(arguments_as_list$lambda))
    arguments_as_list$theta = unlist(as.data.frame(arguments_as_list$theta))
    arguments_as_list$vars_use = as.character(unlist(as.data.frame(arguments_as_list$vars_use)))
    message(paste0(arguments_as_list,collapse = " | "))
  
    ## set some additional input values in param list:
    arguments_as_list[["data_mat"]] = seurat_object@reductions[[reduction_name]]@cell.embeddings
    arguments_as_list[["meta_data"]] = seurat_object@meta.data
    # no verbose
    arguments_as_list[["verbose"]] = FALSE
    # don't return full object
    arguments_as_list[["return_object"]] =FALSE
    # no pca
    arguments_as_list[["do_pca"]] =FALSE
    # add init labels
    arguments_as_list[["cluster_prior"]] =cluster_prior_matrix
    arguments_as_list[["nclust"]] = nclust_K

    message("Before call to harmony")
    ## call harmony and get embedding
    set.seed(global_seed) # better safe than sorry
    harmonyEmbed = do.call(harmony::HarmonyMatrix,args = arguments_as_list)
    
    # save result in file
    message(Sys.time()," Saving result as ",param_df$full_id[i])
    filename = paste0(filepath,param_df$full_id[i],".txt")
    rownames(harmonyEmbed)=rownames(seurat_object@meta.data) # cannot write rownames reliable to file with non base write functions..
    harmonyEmbed = cbind(Cell_ID=rownames(seurat_object@meta.data),harmonyEmbed) # so add as extra col
    #utils::write.table(as.data.frame(harmonyEmbed),file =  filename,sep =  "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
  #  readr::write_delim(as.data.frame(harmonyEmbed),path = filename,delim = "\t",col_names = TRUE)
    data.table::fwrite(data.table::as.data.table(harmonyEmbed),file=filename,sep="\t",row.names = F,col.names = TRUE)
    
    
    
  }
  
  
}

##########
### make_harmony_param_df
##########

#' Set harmony params
#' @param seurat_object
#' @param reduction_name
#' @return a dataframe with one row corresponding to a harmony parameter set for one run
#' 

make_harmony_param_df = function(arg_list = NULL,batch_vars="Batch_ID",theta_list=list(c(2)),lambda_list=list(c(1))){
  
  require(dplyr)
  require(magrittr)
  require(purrr)
  ## check that length fit
  if(length(batch_vars)!=length(theta_list) | length(batch_vars)!=length(lambda_list)){
    stop("Length of batchvars must be the same as theta_list and lambda_list")
  }
  # if ok, set names accordingly
  names(theta_list) = batch_vars
  names(lambda_list) = batch_vars
  
  if(is.null(arg_list)){
    message("Using default arg_list")
    arg_list<- list(npcs = c(20),
                    sigma = c(0.1),
                    tau = c(0),
                    block.size = c(0.05),
                    max.iter.harmony = c(10),
                    max.iter.cluster = c(200),
                    epsilon.cluster= c(1e-05),
                    epsilon.harmony= c(1e-04)) 
  }
  #make base df
  param_df = arg_list %>% purrr::cross_df()
  # add batchvars as nested var
  fml=as.list(batch_vars)
  names(fml) = batch_vars
  tmp1=fml  %>% purrr::cross_df() %>% dplyr::mutate(id=1) %>% dplyr::mutate(id=1:length(id)) %>% group_by(id) %>% tidyr::nest(theta = !!batch_vars) %>% ungroup() %>% dplyr::select(-id)
  param_df$vars_use = rep(tmp1,nrow(param_df))
  
  # add theta and lambda_list as nested variables
  a1=theta_list  %>% purrr::cross_df() %>% dplyr::mutate(id=1) %>% dplyr::mutate(id=1:length(id)) %>% group_by(id) %>% tidyr::nest(theta = !!batch_vars) %>% ungroup() %>% dplyr::select(-id)
  a2=lambda_list  %>% purrr::cross_df() %>% dplyr::mutate(id=1) %>% dplyr::mutate(id=1:length(id)) %>% group_by(id) %>% tidyr::nest(lambda = !!batch_vars) %>% ungroup() %>% dplyr::select(-id)
  
  # put lambda and theta together
  a3_list = list()
  for( i in 1:nrow(a2)){
    suppressWarnings(expr = {a3_list[[i]] = cbind(a1,a2[i,])})
  }
  a3 = do.call(rbind,a3_list)
  
  # add nested vars to rest
  a4_list = list()
  for( i in 1:nrow(a3)){
    suppressWarnings(expr = {a4_list[[i]] = cbind(param_df,a3[i,])})
  }
  all_param_df = do.call(rbind,a4_list)
  
  return(all_param_df)
  
}

# Example:
# batch_vars1=c("Batch_ID","Sex")
# arg_list1 = list(npcs = c(20,30),
#                  lambda = c(1),
#                  tau = c(0),
#                  block.size = c(0.05),
#                  max.iter.harmony = c(10),
#                  max.iter.cluster = c(200),
#                  epsilon.cluster= c(1e-05),
#                  epsilon.harmony= c(1e-04))
# theta_list1 = list(c(2,3),c(3,5))
# lambda_list1 = list(c(0.1,0.2),c(0.1))
# 
# unlist(as.data.frame(a1))
# 
# test_param_df = make_harmony_param_df(arg_list1,batch_vars1,theta_list1,sigma_list1)
# test_param_df = make_harmony_param_df()
# 
# # Then in run function
# arg_list = as.list(test_param_df[1,])
# arg_list$sigma = unlist(as.data.frame(arg_list$sigma))
# arg_list$theta = unlist(as.data.frame(arg_list$theta))
# arg_list$vars_use = as.character(unlist(as.data.frame(arg_list$vars_use)))
# arg_list


