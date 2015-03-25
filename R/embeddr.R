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
  eig <- NULL # eigenvectors/values
  
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
  
  if(sum(eig$values == 0) > 1) warning('More than one non-zero eigenvalue - disjoint clusters')

  colnames(M) <- paste('component_', 1:p, sep='')
  M <- data.frame(M)
  M$cell_id <- colnames(W)
  return( M )
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

#' Cluster the resulting embedding
#' 
#' Cluster the embedded representation using either kmeans or mixture models 
#' from the mclust package
#' 
#' @param M The dataframe containing the embedding (in component_1 and component_2)
#' @param k The number of clusters to find in the data
#' @param method Either 'kmeans' or 'mm' to use \code{mclust}
#' 
#' @return The dataframe M with a new numeric variable `cluster` containing the assigned cluster
cluster_embedding <- function(M, k = 3, method=c('kmeans','mm')) {
  library(dplyr)
  M_xy <- select(M, component_1, component_2)
  method <- match.arg(method)
  if(method == 'kmeans') {
    km <- kmeans(M_xy, k)
    M$cluster <- km$cluster
  } else if(method == 'mm') {
    library(mclust)
    mc <- Mclust(M_xy, G=k)
    M$cluster <- mc$classification
  }
  return( M )
}


fit_pseudotime_thinning <- function(M, clusters=NULL, ...) {
  library(dplyr)
  library(igraph)
  library(Matrix)
  load_all("/net/isi-scratch/kieran/embeddr/curver/")

  Mp <- M
  if(!is.null(clusters)) Mp <- filter(Mp, cluster %in% clusters)
  Y <- curver::reconstruct(select(Mp, component_1, component_2), niter=1, h=0.2)
  #Y <- curver::reconstruct(select(Mp, component_1, component_2), ...)
  
  D <- as.matrix(dist(Y))
  g <- graph.adjacency(D, weighted=TRUE, mode='undirected')
  g_mst <- minimum.spanning.tree(g)
  A <- as.matrix(get.adjacency(g_mst))
  endpoints <- which(rowSums(A) == 1)
  ordering <- get.shortest.paths(g_mst, from=endpoints)
  paths <- ordering$vpath
  best_path <- paths[[ which.max(sapply(paths, length)) ]]
  
  Z <- Y[ordering$vpath[[1]],]
  
  ## now we have the ordering want to work out the arc-length
  n <- dim(Z)[1]
  Z_start <- Z[1:n-1,]
  Z_end <- Z[2:n,]
  Z_diff <- Z_end - Z_start
  pst <- sqrt(rowSums(Z_diff^2))
  pseudotime <- c(0, cumsum(pst))
  Mp$pseudotime <- pseudotime[invPerm(as.integer(ordering$vpath[[1]]))]
  Mp$trajectory_1 <- Y[,1]
  Mp$trajectory_2 <- Y[,2]
  return(Mp)
}

#' Fit the pseudotime curve
#' 
#' Fits the pseudotime curve using principal curves from the princurve library
#' 
#' @param M The dataframe containing the embedding
#' @param clusters The (numeric) clusters to use for the curve fitting. If NULL (default) then
#' all points are used
#' 
#' @return The dataframe \code{M} with three new variables:
#' \describe{
#' \item{pseudotime}{The pseudotime of the cell (arc-length from beginning of curve)}
#' \item{trajectory_1}{The x-coordinate of a given cell's projection onto the curve}
#' \item{trajectory_2}{The y-coordinate of a given cell's projection onto the curve}}
fit_pseudotime <- function(M, clusters = NULL, ...) {
  library(dplyr)
  library(princurve)
  Mp <- M
  if(!is.null(clusters)) Mp <- filter(M, cluster %in% clusters)
  pc <- principal.curve(x = as.matrix(select(Mp, component_1, component_2)), ...)
  Mp$pseudotime <- pc$lambda
  Mp$trajectory_1 <- pc$s[,1]
  Mp$trajectory_2 <- pc$s[,2]
  #Mp <- arrange(Mp, pseudotime)
  return( Mp )
}

#' Plot the cells in the embedding
#' 
#' This function takes a data frame with at least positional (component_0 & component_1) information
#' and plots the resulting embedding. If clusters are assigned it can colour by these, and if a pseudotime
#' trajectory is assigned it will plot this through the embeddding.
#' 
#' @param M The dataframe containing the embedding
#' @param color_by The variable to color the embedding with (defaults to cluster)
#' 
#' @return A \code{ggplot2} plot
plot_embedding <- function(M, color_by = 'cluster') {
  library(ggplot2)
  
  if('pseudotime' %in% names(M)) M <- arrange(M, pseudotime)
  if(!('cluster' %in% names(M))) color_by <- 'pseudotime'
  
  plt <- ggplot(data=M) + theme_bw()
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

#' Plot cells in pseudotime
#' 
#' Plot a set of genes through pseudotime
#' 
#'  @param Mp A dataframe containing the embedding
#'  @param xp The gene-by-cell matrix of log normalised expression counts
#'  @param genes The genes to use for the embedding
#'  @param short_names Short gene names to display; default NULL and \code{genes} is used for gene names
#'  @param nrow Number of rows of plots; passed to \code{facet_wrap}
#'  @param ncol Number of columns of plots; passed to \code{facet_wrap}
#'  
#'  @return A \code{ggplot2} plot
plot_in_pseudotime <- function(Mp, xp, genes, short_names = NULL, nrow = NULL, ncol = NULL) {
  library(reshape2)
  if(ncol(xp) != nrow(Mp)) stop('xp must be gene-by-cell matrix')  

  xp <- data.frame(t(xp))
  xp <- select(xp, one_of(genes))
  if(!is.null(short_names)) names(xp) <- short_names
  xp$pseudotime <- Mp$pseudotime
  df_x <- melt(xp, id.vars='pseudotime', variable.name='gene', value.name='counts')
  ggplot(data=df_x, aes(x=pseudotime, y=counts)) + geom_point(size=1.5, alpha=0.5) +
    theme_bw() + geom_smooth(method='loess', color='firebrick') + facet_wrap(~ gene, nrow = nrow, ncol = ncol) +
    ylab('Normalised log10(FPKM)')
}

#' Reverse pseudotime
#' 
#' Reverse the pseudotimes of cells
#' 
#' @param M A dataframe containing a pseudotime variable to be reversed
#' @return A dataframe with the reversed pseudotime
reverse_pseudotime <- function(M) {
  reverse <- function(x) -x + max(x) + min(x)
  M$pseudotime <- reverse(M$pseudotime)
  return( M )
}

plot_heatmap <- function(M, x, ...) {
  library(gplots)
  xp <- x[,order(M$pseudotime)]
    
  heatmap.2(xp, dendrogram="none", Colv=FALSE,
            col=redblue(256), trace="none", density.info="none", scale="row", ...)
}

plot_graph <- function(M, W) {
  df <- select(M, x = component_1, y = component_2)
  
  diag(W) <- 0
  locs <- which(W > 0, arr.ind = TRUE)
  from_to <- apply(locs, 1, function(xy) {
    xy_from <- df[xy[1],]
    xy_to <- df[xy[2],]
    xx <- c(xy_from, xy_to)
    ##print(class(xx))
    names(xx) <- c('x_from','y_from','x_to','y_to')
    unlist(xx)
  })
  from_to <- data.frame(t(from_to))
  
  plt <- ggplot() + 
    geom_segment(data=from_to, aes(x=x_from, xend=x_to, y=y_from, yend=y_to), 
                 alpha=0.5, color='grey', linetype=2) +
        theme_minimal() 
  if('cluster' %in% names(M)) {
    df$cluster <- M$cluster
    plt <- plt + geom_point(data=df, aes(x=x,y=y,color=as.factor(cluster))) 
  } else {
    plt <- plt + geom_point(data=df, aes(x=x,y=y))
  }
  plt
}






