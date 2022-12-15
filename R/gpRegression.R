covarianceMatrix <- function(params, X1, X2, S, kernel) {

  m <- length(X1)
  n <- length(X2)
  cov_mat <- matrix(0,nrow=m,ncol=n)

  kernelFunc <- paste(kernel,"Func",sep="")

  for (i in 1:m) { for (j in 1:n) {

    #calculate matrix entry using kernel function
    cov_mat[i,j] <- get(kernelFunc)(X1[i],X2[j],params)

    #add error bar terms if S is not NULL
    if(i==j & !is.null(S)) {cov_mat[i,i] <- cov_mat[i,i] + square(S[i])}
  }}

  if(m==n) cov_mat <- (cov_mat + t(cov_mat))/2 #sanity/symmetry check
  return(cov_mat)
}

marginalLikelihood <- function(params, dataset, kernel, print_terms=FALSE) {

  cov_mat <- covarianceMatrix(params, dataset$X, dataset$X, dataset$S, kernel)
  Y <- dataset$Y
  n <- length(Y)

  term1 <- t(Y)%*%calcInverse(cov_mat)%*%Y
  term2 <- log(det(cov_mat))

  if(print_terms) {
    cat(kernel,"datafit term: ", -0.5*term1,"\n")
    cat(kernel,"complexity term: ", -0.5*term2,"\n")
  }

  ll <- -0.5*(term1+term2+(n*log(2*pi)))

  return(ll)
}

negativeMLL <- function(params, dataset, kernel) {
  return(-1*marginalLikelihood(params, dataset, kernel))
}

getHyperparameters <- function(dataset, kernel, Q=NULL) {

  fn.val <- NULL
  res <- NULL
  restarts <- 1

  kernelBounds <- paste(kernel,"Bounds",sep="")
  b <- get(kernelBounds)(dataset, Q)

  verbose <- as.integer(b$itrs/5)
  if(kernel=="diagonal") {verbose <- FALSE}
  set.seed(29) #for reproducibility

  for(i in 1:restarts) {

    outDEoptim <- DEoptim(negativeMLL, b$lower, b$upper,
                          DEoptim.control(NP = 10*length(b$lower),
                          itermax = b$itrs, strategy = 2,
                          F = 1.2, CR = 1, trace = verbose),
                          dataset = dataset, kernel = kernel)

    # parallelType=1,packages=c(),
    # parVar=c('marginalLikelihood','covarianceMatrix',
    #          'diagonalFunc','square','calcInverse',
    #          'stationaryFunc','stationary.Q',
    #          'nonstationaryFunc','nonstationary.Q',
    #          'nonstationaryHypVals','nonstationaryHypForm')

    crt.val <- outDEoptim$optim$bestval
    crt.par <- outDEoptim$optim$bestmem

    if(i==1){
      fn.val <- crt.val
      res <- list("par"=crt.par, "value"=crt.val)
    } else if(crt.val < fn.val) {
      fn.val <- crt.val
      res <- list("par"=crt.par, "value"=crt.val)
    }
  }

  cat(kernel,"hyperparameters: ", res$par,", value: ", res$value,"\n")
  return(res)
}

getPosterior <- function(dataset, res, kernel, plotting=FALSE, Q=NULL) {

  n <- length(dataset$X)
  m <- n*10

  test <- seq(min(0.5,(min(dataset$X)+1)), max(dataset$X)+1, length=m)
  train <- dataset$X

  #plotting samples from the prior distribution :D
  #if(plotting) plotPriorSamples(dataset, kernel, Q)

  K <- covarianceMatrix(res$par, dataset$X, dataset$X, dataset$S, kernel)
  K_inv <- calcInverse(K)
  K_s <- covarianceMatrix(res$par, dataset$X, train, NULL, kernel)
  mu_s <- t(K_s)%*%K_inv%*%dataset$Y

  K_x <- covarianceMatrix(res$par, dataset$X, test, NULL, kernel)
  K_xx <- covarianceMatrix(res$par, test, test, NULL, kernel)

  mu_x <- t(K_x)%*%K_inv%*%dataset$Y
  cov_x = K_xx - (t(K_x)%*%(K_inv)%*%(K_x))

  confidence_bounds <- rep(0,m)
  for (i in 1:m) confidence_bounds[i] <- 2*sqrt(cov_x[i,i])

  # write posterior to a csv file
  df <- data.frame(test, mu_x, mu_x+confidence_bounds, mu_x-confidence_bounds)
  names(df) <- c("Time", "Posterior Mean", "Lower CB", "Upper CB")
  write.csv(x=df, file=paste(kernel,"PosteriorData.csv",sep=""))
  

  if(plotting) {

    #setting up Y co-ordinate limits for plotting
    plot_Ymin <- min(min(mu_x-confidence_bounds), min(dataset$Y))
    plot_Ymax <- max(max(mu_x+confidence_bounds), max(dataset$Y))

    #ggplot plot plot plotting
    ggdf <- data.frame(test, mu_x, mu_x+confidence_bounds, mu_x-confidence_bounds)
    names(ggdf) <- c("Time","Mean","UB","LB")

    origdf <- data.frame(dataset$X, dataset$Y, dataset$S)
    names(origdf) <- c("X", "Y", "S")

    p <- ggplot(ggdf) +
      geom_line(aes(x=Time,y=UB), colour="#619CFF", size=0.81) +
      geom_line(aes(x=Time,y=LB), colour="#619CFF", size=0.81) +
      geom_errorbar(data=origdf, aes(x=X, y=Y, ymin=Y-S, ymax=Y+S),
                    size=0.81, width=0.8, color="#00BA38",
                    position=position_dodge(0.05)) +
      geom_line(aes(x=Time,y=Mean), colour="#F8766D", size=0.81) +
      geom_point(data=origdf, aes(x=X,y=Y), size=2, colour="black") +
      theme_bw() +
      xlab('Time (hrs)')+
      ylab('Arbitrary units') +
      theme(axis.line.x = element_line(color="black", linewidth = 0.5),
            axis.line.y = element_line(color="black", linewidth = 0.5),
            axis.text=element_text(size=20),
            axis.title=element_text(size=30)) +
      # theme(axis.title.x=element_blank(),
      #       axis.text.x=element_blank(),
      #       axis.title.y=element_blank(),
      #       axis.text.y=element_blank())
      theme(legend.position = "none")

    print(p)

  }

  #plotting samples from the posterior distribution C|:D
  #if(plotting) plotPosteriorSamples(dataset, test, mu_x, cov_x, confidence_bounds)

  return(K_inv)
}

plotPriorSamples <- function(dataset, kernel, Q=NULL) {

  num_samples <- 3

  n <- length(dataset$X)
  m <- n*2
  test <- seq(min(1,(min(dataset$X)+1)), max(dataset$X)+1, length=m)

  kernelBounds <- paste(kernel,"Bounds",sep="")
  b <- get(kernelBounds)(dataset, Q)
  init_pars <- 0.5*(b$lower + b$upper)
  #init_pars <- b$upper

  #plotting of prior samples
  #setting up mean, covariance, samples and bound variables

  K_prior <- covarianceMatrix(init_pars, test, test, NULL, kernel)
  mu_prior <- rep(0, m)

  prior_samples <- rmvnorm(n = num_samples+5, mean=mu_prior, sigma=K_prior)
  prior_bounds <- rep(0,m)
  for (i in 1:m) prior_bounds[i] <- 2*sqrt(K_prior[i,i])

  prior <- ggplot()
  prior_colors <- hue_pal()(num_samples+2)

  #plotting the random samples from the prior

  for(i in 1:num_samples) {

    crt_df <- data.frame(test, prior_samples[2*i-1,])
    names(crt_df) <- c("Time", "Sample")
    prior <- prior + geom_line(data=crt_df, aes(x=Time, y=Sample),
                               colour = "#50C878",
                               size=0.81,
                               linetype="dashed")
  }

  #plotting mean, LB and UB of prior samples

  data_df <- data.frame(test, mu_prior, mu_prior-prior_bounds, mu_prior+prior_bounds)
  names(data_df) <- c("Time", "Mean", "LB", "UB")

  prior <- prior + geom_line(data=data_df, aes(x=Time, y=LB), colour = "#619CFF", size=0.81)
  prior <- prior + geom_line(data=data_df, aes(x=Time, y=UB), colour = "#619CFF", size=0.81)
  prior <- prior + geom_line(data=data_df, aes(x=Time, y=Mean), colour = "#F8766D", size=1)

  #adding configs for prior samples plot
  #dY_extra <- 0.5*(max(dataset$Y) - min(dataset$Y))

  prior <- prior + theme_bw() +
    #coord_cartesian(ylim = c(min(dataset$Y)-dY_extra, max(dataset$Y)+dY_extra)) +
    xlab('Time (hrs)')+
    ylab('Arbitrary units') +
    scale_color_manual(name='GP Prior',
                       breaks=c('Mean', 'Samples', 'Confidence Bounds'),
                       values=c('Mean'="#F8766D", 'Samples'="#00BA38", 'Confidence Bounds'="#619CFF")) +
    theme(legend.position = "bottom") +
    theme(axis.line.x = element_line(color="black", linewidth = 0.5),
          axis.line.y = element_line(color="black", linewidth = 0.5),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20)
          # axis.text = element_blank(),
          # axis.title = element_blank(),
          # axis.ticks = element_blank(),
          )

  print(prior)

  #prior samples plotting done
}

plotPosteriorSamples <- function(dataset, test, mu_x, cov_x, bounds_x) {

  num_samples <- 3
  post_samples <- rmvnorm(n = num_samples+5, mean=mu_x, sigma=cov_x)

  post <- ggplot()
  post_colors <- hue_pal()(num_samples+2)

  #plotting random samples from posterior

  for(i in 1:num_samples) {

    crt_df <- data.frame(test, post_samples[2*i-1,])
    names(crt_df) <- c("Time", "Sample")
    post <- post + geom_line(data=crt_df, aes(x=Time, y=Sample),
                             colour = "#50C878",
                             size=0.81,
                             linetype="dashed")
  }

  #plotting posterior mean + bounds

  data_df <- data.frame(test, mu_x, mu_x-bounds_x, mu_x+bounds_x)
  names(data_df) <- c("Time", "Mean", "LB", "UB")

  post <- post + geom_line(data=data_df, aes(x=Time, y=LB), colour = "#619CFF", size=0.81)
  post <- post + geom_line(data=data_df, aes(x=Time, y=UB), colour = "#619CFF", size=0.81)
  post <- post + geom_line(data=data_df, aes(x=Time, y=Mean), colour = "#F8766D", size=1)

  #plotting actual points
  pt_df <- data.frame(dataset$X, dataset$Y, dataset$S)
  names(pt_df) <- c("X", "Y", "S")
  post <- post + geom_point(data=pt_df, aes(x=X, y=Y), size=2, colour="black")# +
          # geom_errorbar(data=pt_df, aes(x=X, y=Y, ymin=Y-S, ymax=Y+S), width=.8, color="black",
          #       position=position_dodge(0.05))

  #adding configs for prior samples plot
  #dY_extra <- 0.5*(max(dataset$Y) - min(dataset$Y))

  post <- post + theme_bw() +

    # coord_cartesian(ylim = c(min(dataset$Y)-dY_extra, max(dataset$Y)+dY_extra)) +
    xlab('Time (hrs)')+
    ylab('Arbitrary units') +

    scale_color_manual(name='GP Prior',
                       breaks=c('Mean', 'Samples', 'Confidence Bounds'),
                       values=c('Mean'="#F8766D", 'Samples'="#00BA38", 'Confidence Bounds'="#619CFF")) +

    theme(legend.position = "bottom") +
    theme(axis.line.x = element_line(color="black", linewidth = 0.5),
          axis.line.y = element_line(color="black", linewidth = 0.5),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20)
          # axis.text=element_blank(),
          # axis.title=element_blank(),
          # axis.ticks=element_blank()
          )
  print(post)
}
