sensAnalysisMultPoints <- function(sRandTs, png = FALSE, pngDir = "M:/thesis_MSc/visualizations/sensAnalysis/") {
  sensRandTS <- NULL
  ewsRandTSwins50bw4 <- NULL
  ewsRandTSwins50bw7 <- NULL
  ewsRandTSwins70bw4 <- NULL
  ewsRandTSwins70bw7 <- NULL
  for (randPoint in 1:nrow(sRandTs)) {
    ts <- sRandTs[randPoint,2:ncol(sRandTs)]
    ts <- removeoutliers(ts, mildness = 1.5)
    timeSteps <- (1:length(ts))  # a vector of the time series length to assign as names
    tsNames <- t(rbind(timeSteps, ts))
    colnames(tsNames) <- c('ndvi', 'day16x')
    rownames(tsNames) <- (1:length(ts))
    sensRandTS[[randPoint]] <- sensitivity_ews(tsNames, indicator = 'acf1', interpolate = T,
                                               detrending = 'gaussian',bandwidthrange = c(2, 10), incrbandwidth = 1,
                                               winsizerange = c(10, 75), incrwinsize = 5)
    # print out sensitivity analysis
    if (png == TRUE) {
      if (randPoint <= 15)
        dev.print(png, file = paste0(pngDir, "declineP", randPoint,".png"), width=600, height=600)
      else
        dev.print(png, file = paste0(pngDir, "nDeclineP", randPoint,".png"), width=600, height=600)
      
      # ews to show detreinding effect on TS
      # winsize 50
      ewsRandTSwins50bw4[[randPoint]] <- generic_ews(tsNames, winsize = 50, detrending = 'gaussian', bandwidth = 4, interpolate = T)
      if (randPoint <= 15)
        dev.print(png, file = paste0(pngDir, "declineP", randPoint,"wins50bw4.png"), width=600, height=600)
      else
        dev.print(png, file = paste0(pngDir, "nDeclineP", randPoint,"wins50bw4.png"), width=600, height=600)
      ewsRandTSwins50bw7[[randPoint]] <- generic_ews(tsNames, winsize = 50, detrending = 'gaussian', bandwidth = 7, interpolate = T)
      if (randPoint <= 15)
        dev.print(png, file = paste0(pngDir, "declineP", randPoint,"wins50bw7.png"), width=600, height=600)
      else
        dev.print(png, file = paste0(pngDir, "nDeclineP", randPoint,"wins50bw7.png"), width=600, height=600)
      # winsize 70
      ewsRandTSwins70bw4[[randPoint]] <- generic_ews(tsNames, winsize = 70, detrending = 'gaussian', bandwidth = 4, interpolate = T)
      if (randPoint <= 15)
        dev.print(png, file = paste0(pngDir, "declineP", randPoint,"wins70bw4.png"), width=600, height=600)
      else
        dev.print(png, file = paste0(pngDir, "nDeclineP", randPoint,"wins70bw4.png"), width=600, height=600)
      ewsRandTSwins70bw7[[randPoint]] <- generic_ews(tsNames, winsize = 70, detrending = 'gaussian', bandwidth = 7, interpolate = T)
      if (randPoint <= 15)
        dev.print(png, file = paste0(pngDir, "declineP", randPoint,"wins70bw7.png"), width=600, height=600)
      else
        dev.print(png, file = paste0(pngDir, "nDeclineP", randPoint,"wins70bw7.png"), width=600, height=600)
      graphics.off()
    }
  }
  return(list(sensRandTS, ewsRandTSwins50bw4, ewsRandTSwins50bw7, ewsRandTSwins70bw4, ewsRandTSwins70bw7))
}