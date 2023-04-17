detectTimePeriod <- function(d, res, kernel, K_inv) {

  m <- dataMeasures(d)
  z <- 10 #sample freq in per-timestep (Hz for all iap)
  delta <- m$timestep/z

  duration <- (max(d$X) - min(d$X))*10
  sample <- seq(0, duration, by = delta)
  K_z <- covarianceMatrix(res$par, d$X, sample, NULL, kernel)
  mu_z <- t(K_z)%*%K_inv%*%d$Y

  z.spec <- spectrum(mu_z,log="no",span=10,plot=FALSE)
  f <- z.spec$freq/delta
  P <- 2*z.spec$spec

  peak_powrs <- P[localMaxima(P)]
  peak_freqs <- f[localMaxima(P)]

  if (kernel == "stationary") {Q <- stationary.Q(res$par)}
  else if (kernel == "nonstationary") {Q <- nonstationary.Q(res$par)}
  else Q <- 1

  freq_strength <- order(peak_powrs, decreasing = TRUE)
  max_freqs <- peak_freqs[freq_strength]
  cat("Strongest detected time periods are",(1./max_freqs)[1:Q],"\n")
}

calculateStatistic <- function(altMLL, nulMLL) {

  #these are already -log likelihoods
  ratio <- -2*(nulMLL - altMLL)
  statistic <- exp(ratio/2)
  return(statistic)
}

modelSelection <- function(statistic, threshold) {

  cat("Bayes factor is",statistic,"and threshold taken is",threshold,"\n")
  cat("If k trajectories from the same dataset are being analyzed simultaneously, a Bonferroni-type multiple-hypothesis correction can be applied by multiplying each Bayes factor, obtained independently, with (1-0.5^(1/k))/(0.5^(1/k)).\n")
  return(statistic > threshold)
}
