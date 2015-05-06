
## compare monocle and embeddr's differential expression across pseudotime
set.seed(123)

library(monocle)
library(devtools)
library(gplots)

if(!require(scater)) {
  install_github('davismcc/scater')
  library(scater)
}

if(!require(embeddr)) {
  install_local('/net/isi-scratch/kieran/embeddr/embeddr')
  library(embeddr)
}

data(HSMM)

load('/net/isi-scratch/kieran/embeddr/embeddr/data/sce_23.Rdata')
sce <- sce_23
sce@lowerDetectionLimit <- log10(0.1 + 1)

min_fpkm <- 1

## select genes for cds
cds <- HSMM[,pData(HSMM)$State != 3]
min_cells <- 0.2 * dim(cds)[2]
genes_test_cds <- rowSums(exprs(cds) > min_fpkm) > min_cells

## select genes for sce
min_cells_sce <- 0.2 * dim(sce_23)[2]
genes_test_sce <- rowSums(exprs(sce_23) > log10(min_fpkm + 1)) > min_cells_sce

gene_list <- list(cds = names(which(genes_test_cds)), 
                  sce = names(which(genes_test_sce)))

genes_to_use <- intersect(gene_list[[1]], gene_list[[2]])
## genes_to_use <- sample(genes_to_use, 100)

## downsample cds to same number of cells
cds <- cds[,sample(dim(cds)[2], dim(sce)[2])]

sig_tests <- mclapply(list(cds = cds[genes_to_use,], sce = sce[genes_to_use,]), function(z) {
  if(class(z) == 'CellDataSet') {
    return(monocle::differentialGeneTest(z, fullModelFormulaStr = 'expression~sm.bs(Pseudotime)', cores = 2))
  } else if(class(z) == 'SCESet') {
    return(embeddr::pseudotime_test(z, n_cores = 2))
  }
}, mc.cores = 2)

save(sig_tests, file='/net/isi-scratch/kieran/embeddr/embeddr/data/sig_test_comparison.Rdata')


load('/net/isi-scratch/kieran/embeddr/embeddr/data/sig_test_comparison.Rdata')

cds_sig <- rownames(sig_tests$cds)[sig_tests$cds$qval < 0.05]
sce_sig <- sig_tests$sce$gene[sig_tests$sce$q_val < 0.05]

venn(list(cds=cds_sig, sce=sce_sig))

