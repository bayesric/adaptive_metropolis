adaptiveMetropolis <- function(p,t0,t1,nm,sd_mh,ldens,...){
  # Adaptive Metropolis algorithm (am; Haario et al. 2001)
  # Args:
  #   p: dimension of the target distribution.
  #   t0: let the sampler run for t0 steps first under naive Metropolis algorithm,
  #       the samples before step t0 will not be used in the adaptation.
  #   t1: time when all samples afterwards will be used in the adaptation.
  #   nm: number of Metropolis samples.
  #   sd_mh: standard deviation for the proposal distribution of the naive Metropolis,
  #           note that the proposal is a homogeneous Gaussian.
  #   ldens: log density function of the target distribution.
  #   ...: optional arguments to ldens.
  
  suppressMessages(require(mvtnorm))
  # preparation
  sd = 2.4^2/p
  eps = 0.1
  sps = NULL
  sp = rep(0,p)
  C0 = diag(sd_mh,p)
  for (iter in 1:nm){
    if (iter%%(nm/10)==0){
      cat("iteration",iter,"\n")
    }
    # propose new sample
    if (iter <= t1){
      nsp = c(rmvnorm(1,sp,C0))
    }else {
      nsp = c(rmvnorm(1,sp,Ct))
    }
    lalpha = ldens(nsp,...)-ldens(sp,...)
    if (log(runif(1))<lalpha){
      sp = nsp
    }
    # store the samples
    sps = cbind(sps,sp)
    # update covariance estimate
    if (iter == t0+1){
      mut = sp
      Ct = sd*eps*diag(1,p)
    }
    if (iter > t0+1){
      citer = iter-t0-1
      tmpmu = (mut*(citer-1)+sp)/citer
      Ct = (citer-1)*Ct/citer+sd/citer*(citer*mut%*%t(mut)-(citer+1)*tmpmu%*%t(tmpmu)+
                                          sp%*%t(sp)+eps*diag(1,p))
      mut = tmpmu
    }
  }
  return(sps)
}
