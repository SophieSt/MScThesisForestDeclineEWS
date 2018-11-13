rasterEWS <- function(rasterObj, timeStepName, seriesName, doInterpolate = FALSE, hasNames = TRUE,
                      detrending = "no", bandwidth = NULL, winsize = 50, span = NULL, degree = NULL,
                      logtransform = FALSE, AR_n = FALSE, powerspectrum = FALSE) {
  # takes rasterObj with time steps as bands
  tryCatch ({
    nPix <- ncell(rasterObj)  # number of pixels
    timeLen <- nlayers(rasterObj) # length of the time series
    timeInd <- 1 : timeLen  # a vector of the time series length to assign as names
    if (!hasNames) { # but only do that if no names there yet
      names(rasterObj) <- timeInd
    }
  }, error = function(error_condition)  { # in case of no layers (i.e. matrix)
    nPix <- nrow(rasterObj)
    timeLen <- ncol(rasterObj) # length of the time series
    timeInd <- 1 : timeLen  # a vector of the time series length to assign as names
    if (!hasNames) { 
      rownames(rasterObj) <- timeInd
    }
  })
  

  # prepare so that time series: numeric matrix with first column time index, second time series values
  # use of headings
  namesEWS <- c(timeStepName, 
                "ar1", "sd", "sk", "kurt", "cv", "returnrate", "densratio", "acf1") 
  # preallocate matrices to convert to output raster later on
  imageEWS <- matrix(NA, length(namesEWS), nPix) 
  tauEWS <- matrix(NA, length(namesEWS), nPix) 
  rownames(imageEWS) <- namesEWS; rownames(tauEWS) <- namesEWS
  
  for (pix in 1:nPix) {
    y <- t(rasterObj[pix]) # time series per pixel, as vector
    ts <- removeoutliers(y) # quality check
    
    # get time step names
    if (is.raster(rasterObj)){
      nm <- names(rasterObj)
    } else {
      nm <- rownames(rasterObj)
    }
    
    # time series: matrix with time steps & values
    tsMatrix <- data.frame(ts, row.names = nm)
    tsMatrix <- cbind(timeInd, tsMatrix)
    names(tsMatrix) = c(timeStepName, seriesName)
    
    # empty time series: only calculate EWS if more than 2 actual values in time series
    if (sum(!is.na(y)) > 2) {
      pixEWS <- generic_ews(tsMatrix, winsize = winsize, detrending = detrending, bandwidth = bandwidth, interpolate = doInterpolate,
                            span = span, degree = degree, logtransform = logtransform, AR_n = AR_n, powerspectrum = powerspectrum)
      #system("rm Rplots*.pdf")
      ewsMatrix <- as.matrix(pixEWS)
      
      # calculate Kendall's tau per indicator
      timevec <- seq(1, nrow(ewsMatrix))
      for (indic in 1:(length(namesEWS))) {
        kT <- cor.test(timevec, ewsMatrix[,indic], alternative = c("two.sided"), method = c("kendall"), 
                       conf.level = 0.95)
        tau <- kT$estimate
        tauEWS[indic, pix] <- tau # write tau in tau-matrix
      }
      imageEWS[, pix] <- c(t((pixEWS[nrow(pixEWS), ]))) # take last EWS values (last time step) and write in EWS-matrix
      graphics.off()
    }
  }
  
  # spatial EWS
  # rasterinterp0 <- brick(rasterObj)
  # for (pix in 1:nPix) {
  #   y <- c(t(rasterObj[pix])) # time series per pixel, as vector
  #   ts <- removeoutliers(y) # quality check
  #   posNaN <- which(is.na(ts))
  #   if (sum(!is.na(ts)) > 2) {
  #     YY <- approx(timeInd, ts, n = length(ts), method = "linear")
  #     ts[posNaN] <- YY$y[posNaN]
  #   } else {
  #     ts[posNaN] <- rep(0, length(posNaN))
  #   }
  #   rasterinterp0[pix] <- ts
  # }
  # 
  # rastermat <- list()
  # for (name in names(rasterinterp0)) {
  #   rasterband <- rasterinterp0[[name]]
  #   bandmat <- as.matrix(rasterband)
  #   rastermat[[name]] <- bandmat
  # }


  
  
  # rasterize the output EWS matrix
  rasterEWS <- raster(rasterObj) # empty raster, extent, reoslution etc from input raster object
  #rasterTau <- raster(rasterObj)
  for (band in 1:nrow(imageEWS)) { # got through each indicator
    rasterBand <- raster(rasterObj) 
    rasterBand <- setValues(rasterBand, imageEWS[band,])
    rasterEWS <- addLayer(rasterEWS, rasterBand)   # add EWS signal as raster layer
    tauBand <- raster(rasterObj)
    tauBand <- setValues(tauBand, tauEWS[band,])
    rasterEWS <- addLayer(rasterEWS, tauBand)   # add respective trend stat as raster layer
  }
  
  namesEWSTau <- c(timeStepName, timeStepName, 
                   "ar1", "ar1Tau", "sd", "sdTau", "sk", "skTau", "kurt", "kurtTau", 
                   "cv", "cvTau", "returnrate", "returnratetau", "densratio", "densratioTau", "acf1", "acf1tau") 
  
  names(rasterEWS) <- namesEWSTau # watch out, band 1 & 2 are last time step and according tau...
  return(rasterEWS)
}
