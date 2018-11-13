rasterSpEWS <- function(rasterObj) {
  # calculate spatial variance, skewness and autocorrelation per time step per pixel separately, constructing time series of that
  spVarTs <- raster(rasterObj)
  spSkewTs <- raster(rasterObj)
  spAutoCorTs <- raster(rasterObj)
  for (name in names(rasterObj)){
    spVar <- focal(rasterObj[[name]], w = matrix(1/9, nc = 3, nr = 3), fun = var, na.rm = TRUE)
    spSkew <- focal(rasterObj[[name]], w = matrix(1/9, nc = 3, nr = 3), fun = skewness, na.rm = TRUE)
    spAutoCor <- MoranLocal(rasterObj[[name]])
    
    spVarTs <- addLayer(spVarTs, spVar)
    spSkewTs <- addLayer(spSkewTs, spSkew)
    spAutoCorTs <- addLayer(spAutoCorTs, spAutoCor)
  }
  
  tauVar <- raster(rasterObj)
  values(tauVar) <- NA
  tauSkew <- raster(rasterObj) 
  values(tauSkew) <- NA
  tauAutoCor <- raster(rasterObj)
  values(tauAutoCor) <- NA
  timevec <- c(1:nlayers(rasterObj))
  
  nPix <- ncell(rasterObj)
  for (pix in (1:nPix)) {
    pixVar <- c(t(spVarTs[pix]))
    pixSkew <- c(t(spSkewTs[pix]))
    pixAutoCor <- c(t(spAutoCorTs[pix]))
    pixTs <- c(t(rasterObj[pix]))
    
    # calculate trend statistic individually
    if ((sum(!is.na(pixTs)) > 2) && (sum(!is.na(pixVar)) > 2)){
      kTVar <- cor.test(timevec, pixVar, alternative = c("two.sided"), method = c("kendall"), 
                     conf.level = 0.95)
      tauVar[pix] <- kTVar$estimate
    }
    if ((sum(!is.na(pixTs)) > 2) && (sum(!is.na(pixSkew)) > 2)){  
      kTSkew <- cor.test(timevec, pixSkew, alternative = c("two.sided"), method = c("kendall"), 
                        conf.level = 0.95)
      tauSkew[pix] <- kTSkew$estimate
    }
      
    if ((sum(!is.na(pixTs)) > 2) && (sum(!is.na(pixAutoCor)) > 2)){
      kTAutoCor <- cor.test(timevec, pixAutoCor, alternative = c("two.sided"), method = c("kendall"), 
                        conf.level = 0.95)
      tauAutoCor[pix] <- kTAutoCor$estimate
    }
  }
  tauEWS <- brick(tauVar, tauSkew, tauAutoCor)
  names(tauEWS) <- c('Kendall tau Sp Var', 'Kendall tau Sp Skewness', 'Kendall tau Moran I')
  
  return(tauEWS)
}