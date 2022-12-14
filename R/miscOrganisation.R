library("MASS")
library("base")
library("mnormt")
library("DEoptim")
library("stats")
library("hilbertSimilarity")
library("ggplot2")
library("mvtnorm")
library("scales")

square <- function(a) {
  return(a*a)
}

calcInverse <- function(mat) {
  tryCatch(
    expr = {
      res <- chol2inv(chol(mat))
      return(res)
    },
    error = function(e){ 
      cat("error of calcInverse is called\n")
      res <- ginv(mat)
      return(res)
    }
  )
}

ystepSize <- function(y) {
  y <- y[order(y)]
  n <- length(y)
  
  max.dist <- 0
  for (i in 1:(n-1)) {
    crt.dist <- abs(y[i] - y[i+1])
    if(crt.dist > max.dist) max.dist <- crt.dist
  }
  
  return(max.dist)
}

dataMeasures <- function(d) {
  
  timestep <- (max(d$X)-min(d$X))/(length(d$X) + 1)
  Y_step <- ystepSize(d$Y)     
  
  X_range <- max(d$X) - min(d$X)
  Y_range <- max(d$Y) - min(d$Y)
  Y_std <- sd(d$Y)
  
  measures <- list("timestep"=timestep, "Y_step"=Y_step,
                   "X_range"=X_range, "Y_range"=Y_range,
                   "Y_std"=Y_std)
  return(measures)
}

stationary.Q <- function(params) {
  return ((length(params)-1)/3)
}

nonstationary.Q <- function(params) {
  return ((length(params)-1)/9)
}
  
nonstat.external.Q <- function(params) {
  return ((length(params)-1)/6)
}