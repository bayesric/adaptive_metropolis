# Demo

# Load function
source("https://raw.githubusercontent.com/ericyewang/adaptive_metropolis/master/am_fcn.R")

# Create a multivariate Gaussian target distribution
p = 10 # dimension
Fc = matrix(rnorm(p*p),p,p) # loadings for the covariance matrix
C = Fc%*%t(Fc)+diag(0.1,p) # create a covariance matrix
ldens <- function(x,C){
  return(-t(x)%*%solve(C)%*%x)
}

# Run AM
res = adaptiveMetropolis(p=p,t0=200,t1=500,nm=5000,sd_mh=0.1,ldens=ldens,C=C)

# Plot traceplot
plot(1:5000,res[1,],'l')
