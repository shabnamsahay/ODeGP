extractData <- function(filename, errorbars=FALSE) {

  #' Get data from a csv file
  #'
  #' This function retrieves and pre-processes the dataset from the given file into a suitable format for ODeGP.
  #' @param filename A csv file containing the time-series data in column-wise format. If n un-merged replicates are present, the column headings should be Time, Column 1, Column 2, ... , Column n (with the ith column containing the data from the ith replicate). If replicates have already been merged to produce mean values and error bars, the column headings should be X, Y, S (where X contains time, Y contains the mean values, and S contains the corresponding error bar values).
  #' @param errorbars TRUE if replicates have been merged in the input data (resulting in errorbars already being present), else FALSE.
  #' @export
  #' @examples
  #' test_function()

  #read the csv into a dataframe
  df <- read.csv(filename)

  #convert nans (non-sampled values) to NAs (missing values)
  df[is.nan(df)] <- NA

  if(errorbars) {

    d <- list("X"=df$X, "Y"=df$Y, "S"=df$S)
    return(d)
  }

  #since there are replicates, plotting them:
  plotReplicates(df)

  #number of replicates
  nreps <- length(df) - 1

  #getting the average and the error bars from replicates
  avg <- apply(df[,-1], MARGIN = 1, FUN = mean, na.rm = TRUE)
  err <- apply(df[,-1], MARGIN = 1, FUN = sd, na.rm = TRUE)

  avg[is.nan(avg)] <- NA #convert nans in avg to NAs
  err <- err*(sqrt(nreps-1)/sqrt(nreps)) #fix sample/population value

  d <- list("X"=df$Time, "Y"=avg, "S"=err)

  save_df <- data.frame(d$X, d$Y, d$S)
  names(save_df) <- c("X","Y","S")
  write.csv(x=save_df, file="RawData.csv")

  return(d)
}

plotReplicates <- function(df) {

  nreps <- length(df) - 1
  interval <- df$Time[2] - df$Time[1]
  missing_times <- c()

  #plotting the replicates

  p <- ggplot()
  col_names <- colnames(df)

  for(i in 2:length(col_names)) {

    nth <- df[, col_names[i]]
    crt_df <- data.frame(df$Time, nth)
    names(crt_df) <- c("Time", "Replicate")

    p <- p + geom_point(data=na.omit(crt_df), aes(x=Time, y=Replicate),
                        size=2, colour="black") +
             geom_line(data=na.omit(crt_df), aes(x=Time, y=Replicate),
                       size = 0.81, colour = hue_pal()(nreps)[i-1])

    missing_times <- c(missing_times, crt_df[is.na(crt_df$Replicate),]$Time)
  }

  #plotting rectanges over the missing time points

  missing_times <- unique(missing_times)
  num_miss <- length(missing_times)


  miss_df <- data.frame(xmin = missing_times - rep(0.5*interval, num_miss),
                        xmax = missing_times + rep(0.5*interval, num_miss),
                        ymin = rep(-Inf, num_miss), ymax = rep(Inf, num_miss))

  p <- p + geom_rect(data=miss_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                     alpha = 2/5, fill="gray")

  p <- p + theme_bw() +
    xlab('Time (hrs)')+
    ylab('Arbitrary units') +
    theme(axis.line.x = element_line(color="black", linewidth = 0.5),
          axis.line.y = element_line(color="black", linewidth = 0.5),
          axis.text=element_text(size=20),
          axis.title=element_text(size=30)) +
    theme(legend.position = "none")
  print(p)

}

is.nan.data.frame <- function(x) {do.call(cbind, lapply(x, is.nan))}
