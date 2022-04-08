##########
###  Seurat annoy clustering functions to get distances manually
##########

# from: https://github.com/satijalab/seurat/blob/master/R/clustering.R

# define:
"%||%" <- function(a, b) if (!is.null(a)) a else b

# load:

#library(RANN)
library(Seurat)
library(future)
library(future.apply)


# Internal helper function to dispatch to various neighbor finding methods
#
# @param data Input data
# @param query Data to query against data
# @param k Number of nearest neighbors to compute
# @param method Nearest neighbor method to use: "rann", "annoy"
# @param cache.index Store algorithm index with results for reuse
# @param ... additional parameters to specific neighbor finding method
#
#' @importFrom methods new
#' @importClassesFrom SeuratObject Neighbor
#

NNHelper <- function(data, query = data, k, method, cache.index = FALSE, ...) {
  args <- as.list(x = sys.frame(which = sys.nframe()))
  args <- c(args, list(...))
  results <- (
    switch(
      EXPR = method,
      "rann" = {
        args <- args[intersect(x = names(x = args), y = names(x = formals(fun = nn2)))]
        do.call(what = 'nn2', args = args)
      },
      "annoy" = {
        args <- args[intersect(x = names(x = args), y = names(x = formals(fun = AnnoyNN)))]
        do.call(what = 'AnnoyNN', args = args)
      },
      stop("Invalid method. Please choose one of 'rann', 'annoy'")
    )
  )
  n.ob <- new(
    Class = 'Neighbor',
    nn.idx = results$nn.idx,
    nn.dist = results$nn.dists,
    alg.info = results$alg.info %||% list(),
    cell.names = rownames(x = query)
  )
  if (isTRUE(x = cache.index) && !is.null(x = results$idx)) {
    slot(object = n.ob, name = "alg.idx") <- results$idx
  }
  return(n.ob)
}

# Run annoy
#
# @param data Data to build the index with
# @param query A set of data to be queried against data
# @param metric Distance metric; can be one of "euclidean", "cosine", "manhattan",
# "hamming"
# @param n.trees More trees gives higher precision when querying
# @param k Number of neighbors
# @param search.k During the query it will inspect up to search_k nodes which
# gives you a run-time tradeoff between better accuracy and speed.
# @param include.distance Include the corresponding distances
# @param index optional index object, will be recomputed if not provided
#
AnnoyNN <- function(data,
                    query = data,
                    metric = "euclidean",
                    n.trees = 50,
                    k,
                    search.k = -1,
                    include.distance = TRUE,
                    index = NULL
) {
  idx <- index %||% AnnoyBuildIndex(
    data = data,
    metric = metric,
    n.trees = n.trees)
  nn <- AnnoySearch(
    index = idx,
    query = query,
    k = k,
    search.k = search.k,
    include.distance = include.distance)
  nn$idx <- idx
  nn$alg.info <- list(metric = metric, ndim = ncol(x = data))
  return(nn)
}
# Build the annoy index
#
# @param data Data to build the index with
# @param metric Distance metric; can be one of "euclidean", "cosine", "manhattan",
# "hamming"
# @param n.trees More trees gives higher precision when querying
#
#' @importFrom RcppAnnoy AnnoyEuclidean AnnoyAngular AnnoyManhattan AnnoyHamming
#
AnnoyBuildIndex <- function(data, metric = "euclidean", n.trees = 50) {
  f <- ncol(x = data)
  a <- switch(
    EXPR = metric,
    "euclidean" =  new(Class = RcppAnnoy::AnnoyEuclidean, f),
    "cosine" = new(Class = RcppAnnoy::AnnoyAngular, f),
    "manhattan" = new(Class = RcppAnnoy::AnnoyManhattan, f),
    "hamming" = new(Class = RcppAnnoy::AnnoyHamming, f),
    stop ("Invalid metric")
  )
  for (ii in seq(nrow(x = data))) {
    a$addItem(ii - 1, data[ii, ])
  }
  a$build(n.trees)
  return(a)
}

# Search an Annoy approximate nearest neighbor index
#
# @param Annoy index, built with AnnoyBuildIndex
# @param query A set of data to be queried against the index
# @param k Number of neighbors
# @param search.k During the query it will inspect up to search_k nodes which
# gives you a run-time tradeoff between better accuracy and speed.
# @param include.distance Include the corresponding distances in the result
#
# @return A list with 'nn.idx' (for each element in 'query', the index of the
# nearest k elements in the index) and 'nn.dists' (the distances of the nearest
# k elements)
#
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#
AnnoySearch <- function(index, query, k, search.k = -1, include.distance = TRUE) {
  n <- nrow(x = query)
  idx <- matrix(nrow = n,  ncol = k)
  dist <- matrix(nrow = n, ncol = k)
  convert <- methods::is(index, "Rcpp_AnnoyAngular")
  if (!inherits(x = plan(), what = "multicore")) {
    oplan <- plan(strategy = "sequential")
    on.exit(plan(oplan), add = TRUE)
  }
  res <- future_lapply(X = 1:n, FUN = function(x) {
    res <- index$getNNsByVectorList(query[x, ], k, search.k, include.distance)
    # Convert from Angular to Cosine distance
    if (convert) {
      res$dist <- 0.5 * (res$dist * res$dist)
    }
    list(res$item + 1, res$distance)
  })
  for (i in 1:n) {
    idx[i, ] <- res[[i]][[1]]
    if (include.distance) {
      dist[i, ] <- res[[i]][[2]]
    }
  }
  return(list(nn.idx = idx, nn.dists = dist))
}

##########
### Run clustering
##########

run_seurat_annoy_manually = function(input_data,dim=200,k=20,metric= "cosine"){

  input_data= input_data[,1:max(ncol(input_data),dim)]

  args = list(
    data = input_data,
    k = k,
    metric =metric
  )

  nn.ranked  = do.call(what = 'AnnoyNN', args = args)
  return(nn.ranked)
}
