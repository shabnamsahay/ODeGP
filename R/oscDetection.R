oscOrNot <- function(d, threshold=14, detrend=TRUE, nonstat=TRUE, plotting=FALSE, Q=1) {

  #' Oscillatory or not?
  #'
  #' This function is the main interface function that is called to determine whether a dataset is oscillatory or not.
  #' @param d The dataset to be tested, a list with features X (time points), Y (averaged replicates), S (error bars).
  #' @param threshold The cutoff value for the bayes factor for classifying oscillatory data. Defaults to 14.
  #' @param detrend Do you want to detrend the data via linear regression? Defaults to TRUE. If your data is already detrended OR you do not want to detrend the data at all, set this parameter to FALSE.
  #' @param nonstat Do you want to use the nonstationary kernel? Defaults to TRUE.
  #' @param plotting Do you want a plot of the optimised non-stationary kernel posterior on the input data? Defaults to FALSE.
  #' @export
  #' @examples
  #' test_function()

  #handling missing values
  d <- removeMissing(d)

  #plotting raw data
  #if(plotting) plotSimpleData(d)

  #detrending the data
  if(detrend) d <- detrendByRegression(d)

  #centering data so that m(X)=0 (posterior) is applicable
  d$Y <- d$Y - rep(mean(d$Y), length(d$Y))

  #plotting processed data
  if(plotting) plotSimpleData(d)

  #SM kernel or the nonstationary kernel
  if(nonstat) altKern <- "nonstationary"
  else altKern <- "stationary"

  #fitting the models and comparing the lls
  bayesFactor <- kernelComparison(d, altKern, plotting, Q)

  #decision on oscillatory or not
  osc.or.not <- modelSelection(bayesFactor, threshold)
  cat("==================================================\n")
  if(osc.or.not) cat("The data is oscillatory!\n")
  else cat("The data is not oscillatory.\n")
  cat("==================================================\n")

  return(bayesFactor)
}

kernelComparison <- function(d, altKern, plotting=FALSE, Q=1) {

  nulKern <- "diagonal"
  altComponents <- Q #this is Q

  #null hypothesis fitting
  nulResult <- getHyperparameters(d, nulKern)
  nulInv <- getPosterior(d, nulResult, nulKern, FALSE)

  #alternate hypothesis fitting
  altResult <- getHyperparameters(d, altKern, altComponents)
  altInv <- getPosterior(d, altResult, altKern, plotting, altComponents)
  cat("==================================================\n")
  detectTimePeriod(d, altResult, altKern, altInv)

  altMLL <- marginalLikelihood(altResult$par, d, altKern, print_terms=TRUE)
  nulMLL <- marginalLikelihood(nulResult$par, d, nulKern, print_terms=TRUE)

  return(calculateStatistic(altMLL, nulMLL))
}


detrendByRegression <- function(d) {
  relation <- lm(Y~X, data=d)
  tempv <- data.frame(X = d$X)
  result <- predict(relation, newdata=tempv)
  d$Y <- (d$Y - result)
  return(d)
}

removeMissing <- function(d) {

  inds <- which(is.na(d$Y))
  inds <- inds[order(inds)]

  if(length(inds)==0) return(d)

  for (i in 1:length(inds)) {
    d$X <- d$X[-(inds[i]-(i-1))]
    d$Y <- d$Y[-(inds[i]-(i-1))]
    d$S <- d$S[-(inds[i]-(i-1))]
  }

  return(d)
}

plotSimpleData <- function(d) {

  #plotting the raw data
  rawdf <- data.frame(d$X, d$Y,d$S)
  names(rawdf) <- c("X", "Y", "S")

  # Default point plot w/error bars
  p<- ggplot(rawdf, aes(x=X, y=Y)) +
    geom_errorbar(aes(ymin=Y-S, ymax=Y+S), size=0.81, width=0.8,
                  position=position_dodge(0.05), color="#00BA38") +
    # geom_line(colour="#F8766D", size=0.81) +
    geom_point(size=2)

  # miss_df <- data.frame(xmin = d$X,
  #                       xmax = d$X,
  #                       ymin = rep(0, length(d$Y)),
  #                       ymax = d$Y)

  #p <- p + geom_segment(data=miss_df, aes(x = xmin, y = ymin, xend = xmax, yend = ymax))
  #dY_extra <- 0.5*(max(d$Y) - min(d$Y))

  p <- p +  theme_bw() +

    #coord_cartesian(ylim = c(min(d$Y)-dY_extra, max(d$Y)+dY_extra)) +
    xlab('Time (hrs)')+
    ylab('Arbitrary units') +

    theme(#plot.title = element_text(size = 20, hjust = 0.5),
          axis.line.x = element_line(color="black", linewidth = 0.5),
          axis.line.y = element_line(color="black", linewidth = 0.5),
          axis.text=element_text(size=20),
          axis.title=element_text(size=30)
          # axis.text=element_blank(),
          # axis.ticks=element_blank(),
          # axis.title=element_blank()
          ) +
    theme(legend.position = "none")
  print(p)
}
