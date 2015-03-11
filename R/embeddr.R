# 
# R functions for Laplacian Eigenmaps
# See e.g.
#   'A Tutorial on Spectral Clustering', von Luxburg Statistics and Computing, 17 (4), 2007.
# 
# kieran.campbell@sjc.ox.ac.uk
#
#

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


laplacian_eigenmap <- function(W, type=c('norm','unorm'), p=2) {
  type <- match.arg(type)
  if(nrow(W) != ncol(W)) {
    print('Input weight matrix W must be symmetric')
    return( NULL )
  }
  
  L  <- diag(rowSums(W)) - W
  l <- nrow(L)
  
  if(type == 'norm') {
    Ds <- diag(1/sqrt(rowSums(W)))
    L_sym <- Ds %*% L %*% Ds
    eig <- eigen(L_sym, symmetric=T) # values ordered in decreasing order
    return( diag(Ds) * eig$vectors[,(l-1):(l-p)])
  }
  if(type == 'unorm') {
    L  <- diag(rowSums(W)) - W
    eig <- eigen(L, symmetric=T) # values ordered in decreasing order
    l <- nrow(L)
    return( eig$vectors[,(l-1):(l-p)])  
  }
  return( NULL )
}
