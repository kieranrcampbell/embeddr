
library(VGAM)

alpha <- 0
beta <- 1
thresh <- 3

x <- runif(100, 0, 10)
y <- x + rnorm(100, 0, 1)

y[y < thresh] <- thresh

fit <- vgam(y ~ x, family=tobit(Lower=thresh))
fit_null <- vgam(y ~ 1, family=tobit(Lower=thresh))


## now move everything up by 1

vert_dev <- 1
yp <- y + vert_dev 
threshp <- thresh + vert_dev

fitp <- vgam(yp ~ x, family=tobit(Lower=threshp))
fitp_null <- vgam(yp ~ 1, family=tobit(Lower=threshp))

VGAM::lrtest(fit, fit_null)
VGAM::lrtest(fitp, fitp_null) # exactly the same




