# 
# R functions for Laplacian Eigenmaps
# See e.g.
#   'A Tutorial on Spectral Clustering', von Luxburg Statistics and Computing, 17 (4), 2007.
# 
# kieran.campbell@sjc.ox.ac.uk
#
#

#' Construct a weighted graph adjacency matrix
#' 
#' @param x A k by n matrix for n samples with k features (probably transpose of what you would expect)
#' @param kernel The choice of kernel. 'nn' will give nearest neighbours, 'dist' gives minimum distance and
#' 'heat' gives a heat kernel. Discussed in detail in 'Laplacian Eigenmaps and Spectral Techniques for Embedding and Clustering',
#' Belkin & Niyogi
#' @param nn Number of nearest neighbours if kernel == 'nn'
#' @param eps Maximum distance parameter if kernel == 'dist'
#' @param t 'time' for heat kernel if kernel == 'heat'
#' 
#' @return An n by n adjacency matrix
weighted_graph <- function(x, kernel=c('nn','dist','heat'), 
                           nn = 20, eps = NULL, t = NULL) {
  ## sanity checking
  kernel <- match.arg(kernel)
  
  ## compute distance matrix
  dm <- as.matrix(dist(t(x)))
  
  W <- NULL
  if(kernel == 'nn' || kernel == 'heat') {
    W <- apply(dm, 2, function(r) {
      ind <- order(r)[1:nn]
      rep_vec <- rep(0, length(r))
      rep_vec[ind] <- 1
      names(rep_vec) <- names(r)
      return( rep_vec )
    })
    ## symmetrise as A -> B => B -> A
    W <- 0.5 * (W + t(W))
  }
  
  if(kernel == 'heat') {
    We <- exp(-dm * dm / t) 
    W <- W * We
  }
  
  if(kernel == 'dist') {
    W <- apply(dm, 1, function(r) as.numeric(r < eps))
  }
  return( W )  
}

#' Laplacian eigenmaps
#' 
#' Construct a laplacian eigenmap embedding
#' @param W Weight matrix
#' @param type Type of laplacian eigenmap (norm for normalised, unorm otherwise)
#' @param p Dimension of the embedded space, default is 2
#' 
#' @return The p-dimensional embedding
laplacian_eigenmap <- function(W, type=c('norm','unorm'), p=2) {
  type <- match.arg(type)
  if(nrow(W) != ncol(W)) {
    print('Input weight matrix W must be symmetric')
    return( NULL )
  }
  
  L  <- diag(rowSums(W)) - W
  l <- nrow(L)
  M <- NULL # object to be returned
  
  if(type == 'norm') {
    Ds <- diag(1/sqrt(rowSums(W)))
    L_sym <- Ds %*% L %*% Ds
    eig <- eigen(L_sym, symmetric=T) # values ordered in decreasing order
    M <- diag(Ds) * eig$vectors[,(l-1):(l-p)]
  }
  if(type == 'unorm') {
    L  <- diag(rowSums(W)) - W
    eig <- eigen(L, symmetric=T) # values ordered in decreasing order
    l <- nrow(L)
    M <- eig$vectors[,(l-1):(l-p)]  
  }
  colnames(M) <- paste('component_', 1:p, sep='')
  rownames(M) <- colnames(W)
  return( data.frame(M) )
}

#' Plot the degree distribution of the weight matrix
#' 
#' @param W Weight matrix
#' @param ignore_weights If TRUE weights are discretised to 0 - 1
#' @return A ggplot histogram of weights
plot_degree_dist <- function(W, ignore_weights = FALSE) {
  x <- NULL
  if(ignore_weights) {
    x <- rowSums(ceil(W))
  } else {
    x <- rowSums(W)
  }

  qplot(x) + theme_bw() +
    xlab('Degree') + ylab('Distribution') +
    ggtitle('Degree distribution of weight graph')
}

cluster_embedding <- function(M, k = 3) {
  km <- kmeans(M, 3)
  M$cluster <- km$cluster
  return( M )
}

fit_curve <- function(M, clusters=NULL) {
  library(dplyr)
  library(princurve)
  Mp <- filter(M, cluster %in% clusters)
  pc <- principal.curve(as.matrix(select(Mp, component_1, component_2)))
  Mp$pseudotime <- pc$lambda
  Mp$trajectory_1 <- pc$s[,1]
  Mp$trajectory_2 <- pc$s[,2]
  Mp <- arrange(Mp, pseudotime)
  return( Mp )
}

plot_embedding <- function(M, color_by = 'cluster') {
  library(ggplot2)
  plt <- ggplot(data=M) + theme_bw(base_size=14)
  if(color_by %in% names(M)) {
    if(color_by == 'pseudotime') {
      mapping_str <- color_by
    } else {
      mapping_str <- paste("as.factor(", color_by, ")")
    }
    plt <- plt + geom_point(aes_string(x = "component_1", y = "component_2", 
                                       color=mapping_str))
    if(color_by == 'pseudotime') {
      plt <- plt + scale_color_continuous(name = color_by) 
    } else {
      plt <- plt + scale_color_discrete(name = color_by)
    }
    
    if("pseudotime" %in% names(M)) {
      ## curve has been fit so plot all
      plt <- plt + geom_path(aes(x = trajectory_1, y = trajectory_2), data=M, color='black',
                             size = 1, alpha = 0.7, linetype=2) 
    }
  } else {
    plt <- plt + geom_point(aes(x = component_1, y = component_2))
  }
  return( plt )
}
