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
#' @param symmetrize How to make the adjacency matrix symmetric. Note that slightly
#' counterintuitively, node i having node j as a nearest neighbour doesn't guarantee node
#' j has node i. There are several ways to get round this:
#' \itemize{
#' \item \code{mean} If the above case occurs make the link weight 0.5 so the adjacency matrix becomes \eqn{0.5(A + A')}
#' \item \code{ceil} If the above case occurs set the link weight to 1 (ie take the ceiling of the mean case)
#' \item \code{floor} If the above case occurs set the link weight to 0 (ie take the floor of the mean case)
#' }
#'
#' @return An n by n adjacency matrix
weighted_graph <- function(x, kernel=c('nn','dist','heat'), metric=c('euclidean','correlation'),
                           nn = 20, eps = NULL, t = NULL,
                           symmetrize = c('mean','ceil','floor')) {
  ## sanity checking
  kernel <- match.arg(kernel)
  symmetrize <- match.arg(symmetrize)
  metric <- match.arg(metric)

  ## compute distance matrix
  dm <- NULL
  if(metric == 'euclidean') {
    dm <- as.matrix(dist(t(x)))
  } else if(metric == 'correlation') {
    dm <- cor(x)
  }

  W <- NULL
  if(kernel == 'nn' || kernel == 'heat') {
    if(metric == 'euclidean') {
      W <- apply(dm, 2, function(r) {
        ind <- order(r)[1:nn]
        rep_vec <- rep(0, length(r))
        rep_vec[ind] <- 1
        names(rep_vec) <- names(r)
        return( rep_vec )
      })
    } else if(metric == 'correlation') {
      W <- apply(dm, 2, function(r) {
        ind <- order(r, decreasing = TRUE)[1:nn]
        rep_vec <- rep(0, length(r))
        rep_vec[ind] <- 1
        names(rep_vec) <- names(r)
        return( rep_vec )
      })
    }
    ## symmetrise as A -> B => B -> A
    W <- 0.5 * (W + t(W))
    if(symmetrize == 'ceil') {
      W <- ceiling(W)
    } else if(symmetrize == 'floor') {
      W <- floor(W)
    }
  }

  if(kernel == 'heat') {
    if(metric == 'correlation') stop('Heat hernel not supported for correlation metric')
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

  if(sum(eig$values == 0) > 1) warning(paste('More than one non-zero eigenvalue - disjoint clusters. Multiplicity: ', sum(eig$values == 0)))

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
cluster_embedding <- function(M, k = NULL, method=c('kmeans','mm')) {
  library(dplyr)
  M_xy <- dplyr::select(M, component_1, component_2)
  method <- match.arg(method)
  if(method == 'kmeans') {
    km <- kmeans(M_xy, k)
    M$cluster <- km$cluster
  } else if(method == 'mm') {
    library(mclust)
    mc <- Mclust(M_xy, G=k)
    M$cluster <- as.factor(mc$classification)
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
  if(!is.null(clusters)) Mp <- dplyr::filter(M, cluster %in% clusters)
  pc <- principal.curve(x = as.matrix(dplyr::select(Mp, component_1, component_2)), ...)
  pst <- pc$lambda
  pst <- (pst - min(pst)) / (max(pst) - min(pst))
  Mp$pseudotime <- pst
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
  library(dplyr)
  if(ncol(xp) != nrow(Mp)) stop('xp must be gene-by-cell matrix')

  xp <- data.frame(t(xp))
  xp <- dplyr::select(xp, one_of(genes))
  if(!is.null(short_names)) names(xp) <- short_names
  xp$pseudotime <- Mp$pseudotime

  cn <- 'cluster' %in% names(Mp)
  if(cn) xp$cluster <- Mp$cluster
  id_vars <- 'pseudotime'
  if(cn) id_vars <- c(id_vars, 'cluster')

  df_x <- melt(xp, id.vars=id_vars, variable.name='gene', value.name='counts')

  plt <- NULL
  if(cn) {
    plt <- ggplot(data=df_x, aes(x=pseudotime, y=counts, color=cluster))
  } else {
    plt <- ggplot(data=df_x, aes(x=pseudotime, y=counts))
  }

  return( plt + geom_point(size=1.5) +
    theme_bw() + geom_smooth(method='loess', color='firebrick') + facet_wrap(~ gene, nrow = nrow, ncol = ncol) +
    ylab('Normalised log10(FPKM)') )
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
  library(plyr)
  df <- select(M, x = component_1, y = component_2)

  diag(W) <- 0
  locs <- which((1 * lower.tri(W) * W) > 0, arr.ind = TRUE)
  from_to <- apply(locs, 1, function(xy) {
    xy_from <- df[xy[1],]
    xy_to <- df[xy[2],]
    xx <- c(xy_from, xy_to)
    ##print(class(xx))
    names(xx) <- c('x_from','y_from','x_to','y_to')
    unlist(xx)
  })
  colnames(from_to) <- NULL
  from_to <- data.frame(t(from_to))
  from_to$connect <- plyr::mapvalues(W[locs], c(0.5, 1), c('Single','Both'))
  cols <- c('Single' = 'grey', 'Both' = 'red')

  plt <- ggplot() +
    geom_segment(data=from_to, aes(x=x_from, xend=x_to, y=y_from, yend=y_to, color=connect),
                 alpha=0.5, linetype=1) +
        theme_minimal() + scale_color_manual(values = cols) + geom_point(data=df, aes(x=x,y=y))
  plt + xlab('x') + ylab('y')
  plt
}

## expression~VGAM::bs(Pseudotime, df=3)

#' Fit the gene expression profile in pseudotime
#'
#' This function fits the expression profile of y as a function of
#' pseudotime using a natural cubic spline of degree three. A tobit model
#' is used to censor values less than min_expr
#'
#' @return An object of class VGAM
fit_pseudotime_model <- function(y, t, min_expr) {
  b <- bs(t, df=3)
#   fit <- NULL
#   tryCatch({
#     fit <- suppressWarnings(vgam(y ~ b, family = tobit(Lower = min_expr)))
#   }, error = function(e) {
#     fit <- NULL
#   })
#   return( fit )
   lm(y ~ b)
}

#' Fit the null pseudotime model
#'
#' This function fits the null expression profile in y (ie y ~ 1). A tobit model is
#' used to censor values less than min_expr
#'
#'  @return An object of class VGAM
fit_null_model <- function(y, min_expr) {
  # suppressWarnings(vgam(y ~ 1, family = tobit(Lower = min_expr)))
  lm(y ~ 1)
}

#' Plot the fit in pseudotime
#'
#' @param model Model returned by \code{fit_pseudotime_model} or \code{fit_null_model}
#' @param y A vector gene expression of length number of cells
#' @param t The assigned pseudotime
#' @param min_expr The minimum expression detection threshold
#' @clusters Any clustering by cell type
#'
#' @return An plot object from \code{ggplot}
plot_pseudotime_model <- function(model, y, t, min_expr) {
  df <- data.frame(y=y, t=t, p=predict(model), min_expr = min_expr) ## predict(model)[,1]

  plt <- ggplot(df)  + geom_line(aes(x=t,y=p, color='Predicted')) +
    theme_minimal() + geom_line(aes(x=t, y=min_expr, color='Min expr'), linetype=2) +
    scale_color_manual('', values=c('Min expr' = 'grey', 'Predicted'='red'))
  plt <- plt + geom_point(aes(x=t, y=y))
  return( plt )
}

# plot_pseudotime_models <- function(models, x, min_expr) {
#
# }


#' Perform likelihood ratio test
#'
#' @param model The full model y ~ pseudotime
#' @param null_model The null model y ~ 1
#'
#' @return The p-value
compare_models <- function(model, null_model) {
  require(lmtest)
  require(VGAM)
  lrt <- lmtest::lrtest(model, null_model)
  return( lrt$"Pr(>Chisq)"[2] ) # lmtest
  #lrt@Body["Pr(>Chisq)"][2, ]  # VGAM
}

