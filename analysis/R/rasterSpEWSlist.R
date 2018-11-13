rasterSpEWS <- function(rasterObj, fileout, ncores = 1) {
  # remove outliers from ndvi time series
  rasterOutl <- calc(rasterObj, fun = removeoutliers)
    
  # calculate spatial variance, skewness and autocorrelation per time step per pixel separately, constructing time series of that
  rasterL <- as.list(rasterOutl)
  
  # calculate Spatial Variance & its trend statistic
  spVarTs <- mclapply(rasterL, FUN = spVar, mc.cores = ncores)
  spVarBrick <- stack(spVarTs)
  spVarTau <- calc(spVarBrick, fun = tauL)
  writeRaster(spVarTau, filename = paste0(fileout, 'spVarTau.grd'), datatype = 'INT2S', overwrite = TRUE, maxmemory = 1e+11)
  rm(spVarTs, spVarBrick, spVarTau)
    
  # calculate Spatial Skewness & its trend statistic
  spSkewTs <- mclapply(rasterL, FUN = spSkew, mc.cores = ncores)
  spSkewBrick <- stack(spSkewTs)
  spSkewTau <- calc(spSkewBrick, fun = tauL)
  writeRaster(spSkewTau, filename = paste0(fileout, 'spSkewTau.grd'), datatype = 'INT2S', overwrite = TRUE, maxmemory = 1e+11)
  rm(spSkewTs, spSkewBrick, spSkewTau)
  
  # calculate Spatial Autocorrelation & its trend statistic
  spAutoCorTs <- mclapply(rasterL, FUN = spAutoCor, mc.cores = ncores)
  spAutoCorBrick <- stack(spAutoCorTs)
  spAutoCorTau <- calc(spAutoCorBrick, fun = tauL)
  writeRaster(spAutoCorTau, filename = paste0(fileout, 'spAutoCorTau.grd'), datatype = 'INT2S', overwrite = TRUE, maxmemory = 1e+11)
  rm(spAutoCorTs, spAutoCorBrick, spAutoCorTau)
}



### functions used by the above ###
# calculation of Spatial Variance
spVar <- function(rasterL) {
  spOut <- focal(rasterL, w = matrix(1/9, nc = 3, nr = 3), fun = var, na.rm = TRUE)
  return(spOut)
}
# calculation of Spatial Skewness
spSkew <- function(rasterL) {
  spOut <- focal(rasterL, w = matrix(1/9, nc = 3, nr = 3), fun = skewness, na.rm = TRUE)
  return(spOut)
}
# calculation of Spatial Autocorrelation
spAutoCor <- function(rasterL) {
  spOut <- MoranLocal(rasterL)
  return(spOut)
}

# calculation of trend statistic
tauL <- function(pix) {
  timevec <- (1:length(pix))
  pixTS <- as.vector(pix)
  if ((sum(!is.na(pixTS)) > 100)) { # rewrite into relative length
    y <- na.approx(pixTS, rule = 2)
    kTVar <- cor.test(timevec, y, alternative = c("two.sided"), method = c("kendall"), 
                      conf.level = 0.95)
    tauVar <- as.vector(kTVar$estimate) * 1000
  } else {
    tauVar <- NA
  }
  return(tauVar)
}