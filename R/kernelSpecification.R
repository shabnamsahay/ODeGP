diagonalFunc <- function(x1, x2, params) {
  
  sigma <- params[1]
  
  if (x1==x2) {
    return(square(sigma))
  } else {
    return(0)
  }
}

diagonalBounds <- function(d, Q=NULL) {
  
  m <- dataMeasures(d)
  lower <- c(0)
  upper <- c(m$Y_range)
  
  bounds <- list("lower"=lower, "upper"=upper, "itrs"=50)
  return(bounds)
}

stationaryFunc <- function(x1, x2, params) {
  
  Q <- stationary.Q(params)
  sum <- 0
  
  for (i in 1:Q) {
    decay <- exp(-2*square(pi)*square(params[(3*i)-1])*square(x1-x2))
    periodic <- cos(2*pi*(x1-x2)/params[(3*i)])
    sum <- sum + (square(params[(3*i)-2])*decay*periodic)
  }
  
  if(x1==x2) sum <- sum + square(params[length(params)])
  return(sum)
}

stationaryBounds <- function(d, Q) {
  
  m <- dataMeasures(d)
  
  single_lower <- c(0, 0, 2*m$timestep)
  single_upper <- c(m$Y_range, m$X_range, m$X_range)
  lower <- c()
  upper <- c()
  
  for(j in 1:Q) {
    lower <- c(lower, single_lower)
    upper <- c(upper, single_upper)
  }
  
  #add lb and ub for global diag noise term
  lower <- c(lower, 0)
  upper <- c(upper, m$Y_std)
  
  bounds <- list("lower"=lower, "upper"=upper, "itrs"=200)
  return(bounds)
}

nonstationaryFunc <- function(x1, x2, params) {
  
  Q <- nonstationary.Q(params)
  sum <- 0
  
  for (i in 1:Q) {
    
    hyp <- nonstationaryHypVals(x1, x2, params[((9*i)-8) : (9*i)])
    
    if(hyp$l1 < 0) cat("l1 is",hyp$l1,"and l2 is",hyp$l2,"\n")
    gibbs1 <- sqrt((2*hyp$l1*hyp$l2)/(square(hyp$l1)+square(hyp$l2)))
    gibbs2 <- exp(-1*square(x1-x2)/(square(hyp$l1)+square(hyp$l2)))
    decay <- gibbs1*gibbs2
    
    periodic <- cos(2*pi*((x1/hyp$p1)-(x2/hyp$p2)))
    sum <- sum + hyp$w1*hyp$w2*decay*periodic
  }

  if(x1==x2) sum <- sum + square(params[length(params)])
  return(sum)
}

nonstationaryHypVals <- function(x1, x2, Qpars) {
  
  #Qpars has exactly 9 elements
  
  w1 <- nonstationaryHypForm(Qpars[1:3], x1)
  w2 <- nonstationaryHypForm(Qpars[1:3], x2)
  
  l1 <- nonstationaryHypForm(Qpars[4:6], x1)
  l2 <- nonstationaryHypForm(Qpars[4:6], x2)
  
  p1 <- nonstationaryHypForm(Qpars[7:9], x1)
  p2 <- nonstationaryHypForm(Qpars[7:9], x2)
  
  hyp <- list("w1"=w1, "w2"=w2, "l1"=l1, "l2"=l2, "p1"=p1, "p2"=p2)
  return(hyp)
}

nonstationaryHypForm <- function(innerparams, t) {
  
  pmax <- innerparams[1]
  pmin <- innerparams[2]
  curvature <- innerparams[3]
  
  term2 <- (pmax-pmin)*exp(-1*curvature*t)
  return(pmax - term2)
}

nonstationaryBounds <- function(d, Q) {
  
  m <- dataMeasures(d)
  
  #func params: maximum val, minimum val, curvature
  single_lower <-  c(m$Y_step, 0, 0,          #w func lower bounds
                     m$timestep, 0, 0,  #l func lower bounds
                     m$X_range/2, 2*m$timestep, 0)   #p func lower bounds
  single_upper <- c(m$Y_range, m$Y_step, m$X_range,   #w func upper bounds
                    m$X_range, m$timestep, m$X_range, #l func upper bounds
                    m$X_range, m$X_range/2, m$X_range) #p func upper bounds
  
  lower <- c()
  upper <- c()
  for(j in 1:Q) {
    lower <- c(lower, single_lower)
    upper <- c(upper, single_upper)
  }
  
  #add lb and ub for diag noise term
  lower <- c(lower, 0)
  upper <- c(upper, m$Y_std)
  
  bounds <- list("lower"=lower, "upper"=upper, "itrs"=200)
  return(bounds)
}

exponentialFunc <- function(x1, x2, params) {
  
  sigma <- params[1]
  l <- params[2]
  
  a <- square(sigma)
  b <- exp(-1*square(x1-x2)/(2*square(l)))
  return(a*b)
}

exponentialBounds <- function(d, Q=NULL) {
  
  m <- dataMeasures(d)
  lower <- c(0, m$timestep) #lbs for sigma and l
  upper <- c(m$Y_range, m$X_range) #ubs for sigma and l
  
  bounds <- list("lower"=lower, "upper"=upper, "itrs"=200)
  return(bounds)
}

periodicFunc <- function(x1, x2, params) {
  
  sigma <- params[1]
  l <- params[2]
  p <- params[3]
  eps <- params[4]
  
  a <- square(sigma)
  b <- -2*square(sin(pi*abs(x1-x2)/p))
  c <- square(l)
  
  term <- a*exp(b/c)
  if (x1==x2) term <- term + square(eps)
  
  return(term)
}

periodicBounds <- function(d, Q=NULL) {
  
  m <- dataMeasures(d)
  lower <- c(0, 0, 2*m$timestep)
  upper <- c(m$Y_range, m$X_range, m$X_range)
  
  #add lb and ub for global diag noise term
  lower <- c(lower, 0)
  upper <- c(upper, m$Y_std)
  
  bounds <- list("lower"=lower, "upper"=upper, "itrs"=200)
  return(bounds)
}

linearFunc <- function(x1, x2, params) {
  
  sigma_b <- params[1]
  sigma <- params[2]
  con <- params[3]
  
  a <- square(sigma_b)
  b <- square(sigma)
  con <- (x1-con)*(x2-con)
  return(a + (b*con))
}

linearBounds <- function(d, Q=NULL) {
  
  m <- dataMeasures(d)
  lower <- c(0, 0, 0)
  upper <- c(m$Y_range, m$Y_range, 2*m$X_range)
  
  bounds <- list("lower"=lower, "upper"=upper, "itrs"=200)
  return(bounds)
}