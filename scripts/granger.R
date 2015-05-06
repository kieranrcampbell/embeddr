
## Granger causality for pseudotemporal ordering

library(monocle)
library(devtools)
library(biomaRt)

if(!require(scater)) {
  install_github('davismcc/scater')
  library(scater)
}

if(!require(embeddr)) {
  install_local('/net/isi-scratch/kieran/embeddr/embeddr')
  library(embeddr)
}

load('/net/isi-scratch/kieran/embeddr/embeddr/data/sce_23.Rdata')
sce <- sce_23
sce@lowerDetectionLimit <- log10(0.1 + 1)

mrf <- c('MYOD1', 'MYF5', 'MYF6', 'MYOG','ID1','CDK1', 'SPHK1', 'PDGFRA') 
mrf_ind <- sapply(mrf, match, fData(sce_23)$gene_short_name)

plot_in_pseudotime(sce_23[mrf_ind,])

mrf_models <- fit_pseudotime_models(sce_23[mrf_ind,])
plot_pseudotime_model(sce_23[mrf_ind,], models = mrf_models)

prediction_matrix <- predicted_expression(NULL, mrf_models)
pst <- pData(sce_23)$pseudotime

e1 <- as.vector(exprs(sce_23[mrf_ind[1],]))
##df <- data.frame(pst, y=prediction_matrix[,4])

smooth_moving_average <- function(y, pst, w = 0.1, s=0.01, l=1) {
  ## can show n = (l - w) / s + 1
  nt <- (l - w) / s
  window_lower <- (0:nt) * s
  window_upper <- window_lower + w
  window <- data.frame(lower = window_lower, upper = window_upper)
  points_in_interval <- apply(window, 1, function(lu) {
    return(which(lu[1] < pst & pst < lu[2]))
  })
  values <- sapply(points_in_interval, function(ps) {
    mean(y[ps])
  })
  #return(data.frame(t=rowMeans(window), y=values))
  return(values)
}

smas <- data.frame(apply(prediction_matrix, 2, smooth_moving_average, pst))
names(smas) <- mrf

smas_plot <- smas
smas_plot$pseudotime <- 1:nrow(smas)
smas_melted <- melt(smas_plot, variable.name='gene', id.vars='pseudotime', value.name='exprs')
ggplot(smas_melted) + geom_point(aes(x=pseudotime, y=exprs)) +
  facet_wrap(~gene) + theme_bw()
  
library(lmtest)
library(vars)
VARselect(smas$MYOG, lag.max=30)
## grangertest(MYOG ~ MYOD1, data=smas, order=25)


# poisson process ? -------------------------------------------------------
library(MASS)
library(Renext)

t <- pseudotime(sce_23)
ts <- sort(t)
dts <- diff(ts)
dts <- dts[dts < 0.02]

fit <- fitdistr(dts, 'exponential')
gof <- gofExp.test(dts)
print(gof$p.value)

ggplot(data.frame(pseudotime=dts)) + geom_histogram(aes(x=pseudotime)) #+
  stat_function(fun=dexp, arg=list(rate=fit$estimate), color='red')
