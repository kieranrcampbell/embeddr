### embeddr: Laplacian eigenmaps and principal curves for pseudotemporal ordering and clustering of single-cell RNA-seq data

[Preprint on biorxiv](http://biorxiv.org/content/early/2015/09/18/027219)

`Embeddr` creates a reduced dimensional representation of the gene space using a high-variance gene correlation graph and laplacian eigenmaps. It then fits a smooth pseudotime trajectory using principal curves.

### Quick start

#### Install embeddr

First install embeddr and [scater](https://github.com/davismcc/scater):
```r
library(devtools)
install_github('davismcc/scater')
install_github('kieranrcampbell/embeddr')
```

#### Convert data to an SCESet
If you have a gene-by-cell data.frame `X` of single-cell RNA-seq measurements and a cell-by-feature data.frame `PD` of cell descriptors
then you can create an `SCESet` via
```r
pd <- new('AnnotatedDataFrame', data=PD)
sce <- newSCESet(cellData = X, phenoData = pd)
```

#### Use embeddr
```r
library(embeddr)

## Create a cell-cell correlation graph and use it for the reduced embedding:
sce <- embeddr(sce)

## Plot a reduced-dimension embedding
plot_embedding(sce)

## Optionally cluster the embedding. Cluster assignments are stored in pData(sce)$cluster.
## If no number of clusters is designated, the number is chosen using the BIC from package mclust
sce <- cluster_embedding(sce)

## Fit pseudotime using principal curves
sce <- fit_pseudotime(sce)

## Plot genes 1:10 in pseudotime:
plot_in_pseudotime(sce[1:10,])

## Fit differential expression pseudotime model. This will report gene name, p-val and q-val
diff_gene_test <- pseudotime_test(sce)

```

#### Further examples

Fully worked examples using the Monocle dataset and the distal lung epithelium dataset (Quake et al.) can be found in the [vignettes/](https://github.com/kieranrcampbell/embeddr/tree/master/vignettes) directory.
