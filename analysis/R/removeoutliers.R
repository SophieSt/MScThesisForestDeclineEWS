# this function masks out outliers based on quartiles (< Q1; > Q4)
# mildness for outliers can be given in mildness argument
# higher mildness will result in less outliers masked out!
# default mildness is 1.5


removeoutliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}
