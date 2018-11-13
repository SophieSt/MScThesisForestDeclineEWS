## script to preprocess reference data
# reference data required as: classified forest pixels, classes: decline/no decline
# get data: https://drive.google.com/open?id=118bBGh7G4W_VEVFQIwCFFo3utuTO-1wr

source('R/sensAnalysisMultPoints.R')

ndviCatTS2012 <- brick("CatNdvi2012Modis.tif")
eviCatTS2012 <- brick("CatEvi2012Modis.tif")
ndmiCatTS2012 <- brick("CatNdmi2012Modis.tif")
for (name in names(ndviCatTS2012)) {
  values(ndviCatTS2012[[name]])[values(ndviCatTS2012[[name]]) == -9999] <- NA
}
for (name in names(eviCatTS2012)) {
  values(eviCatTS2012[[name]])[values(eviCatTS2012[[name]]) == -9999] <- NA
}
for (name in names(ndmiCatTS2012)) {
  values(ndmiCatTS2012[[name]])[values(ndmiCatTS2012[[name]]) == -9999] <- NA
}

# stratified random sampling of pixels from decline/no decline
# what is decline, what is nodecline?
# load vector layer of 2012 decline dataset
# get Catalonia borders
spain <- raster::getData('GADM',country='ESP', level=1)
catalonia <- spain[spain@data$NAME_1 == 'Cataluña',]
# reproject catalonia shape (this is the smaller one, should go faster)
catCRS <- spTransform(catalonia, CRS(proj4string(ndviCatTS2012)))

# write Catalonia into shape file
#writeOGR(obj=catCRS, dsn="M://thesis_MSc/data", layer="catBorders", driver="ESRI Shapefile")

# mask out catalonia
ndvi2012 <- mask(ndviCatTS2012, catCRS)
evi2012 <- mask(eviCatTS2012, catCRS)
ndmi2012 <- mask(ndmiCatTS2012, catCRS)

# save reprojected NA-containing raster
#writeRaster(x=tsCat, filename='CataloniaModis2012reproj.tif')

# get declines
declines2012 <- readOGR('M:/thesis_MSc/data/data/declines2012.kml', layer = 'declines2012')
declines2012CRS <- spTransform(declines2012, CRS(proj4string(ndvi2012)))

# set to 1 where decline, 2 where healthy, NA else
declRaster <- rasterize(declines2012CRS, ndvi2012)
values(declRaster)[is.na(values(ndvi2012[[1]])) & !is.na(values(declRaster))] <- NA # declined but no forest, due to fuzzy polygons
values(declRaster)[!is.na(values(ndvi2012[[1]])) & !is.na(values(declRaster))] <- 1 # declined
values(declRaster)[!is.na(values(ndvi2012[[1]])) & is.na(values(declRaster))] <- 2 # not declined
declFactor <- raster(declRaster)
declFactor$declined <- factor(declRaster, levels = c(1,2), labels = c('decline', 'no decline'))

#writeRaster(x = declRaster, filename = 'CataloniaDeclMask2012.tif', datatype = 'INT2S', overwrite = TRUE)
declRaster <- raster('CataloniaDeclMask2012.tif')

# -> decline / no decline mask
# stratify random sampling
set.seed(1234567)
# take 30 pixels from nPix randomly, stratify
sRandomDecl <- sampleStratified(declRaster, size = 15, na.rm=TRUE, sp=TRUE)
classesSRand <- extract(declRaster, sRandomDecl, df = T)
sRandTs <- extract(ndvi2012, sRandomDecl, method = 'simple', df = TRUE)

# wite points into a shape file
#writeOGR(obj=sRandomDecl, dsn="M://thesis_MSc/data", layer="sRandomPointsSensit", driver="ESRI Shapefile")

# plot overview of declines/no declines with stratified random points
spplot(declFactor$declined, col.regions = c('red4', 'green3'), main = 'Overview of Forest Declines 2012',
       sp.layout=list(list("sp.polygons", catCRS, col="black"),
                      list("sp.points", sRandomDecl, cex = 0.5, pch = 7, col = sRandomDecl$layer, colourkey = T)))

# sensitivity EWS now!
sensAnalysis <- sensAnalysisMultPoints(sRandTs)

# stratify random sampling
set.seed(12345678)
# take 150 pixels from nPix randomly, stratify
sRandomRS <- sampleStratified(declRaster, size = 150, na.rm=TRUE, sp=TRUE)
#classesSRand <- extract(declRaster, sRandomDecl, df = T)
sRandTs <- extract(ndvi2012, sRandomRS, method = 'simple', df = TRUE)

# extract generic ews with above settings, extract their tau as well
sRandEWS <- rasterEWS(t(sRandTs[, 2:274]), timeStepName = 'day16x', seriesName = 'ndvi', doInterpolate = TRUE, hasNames = FALSE,
                        detrending = 'gaussian', bandwidth = 4,
                        winsize = 50)
# extract bds test statistic (eliminate False Positives)
# set up model with random forest to predict decline/no decline based on ews & tau, BDS
# extract training error, optimize?
# variable importance
# predict entire study area

# set up a mask for 2016, extract for 2016 ews & tau & BDS
# perform stratified random sampling? if needed
# assess model performance





# drop years 2012, to avoid biased estimate of EWS
#layersToDrop <- 274:nlayers(lsTimeSeries2012)
#lsTimeSeries2012 <- dropLayer(lsTimeSeries2012, layersToDrop)

rasteredEWS <- rasterEWS(ndvi2012, timeStepName = 'day16x', seriesName = 'ndvi', doInterpolate = TRUE, hasNames = FALSE,
                         detrending = 'gaussian', bandwidth = 7,
                         winsize = 50)
ewsR <- rasteredEWS[[c(1,3,5,7,9,11,1,3,15,17)]]
tauR <- rasteredEWS[[c(2,4,6,8,10,12,14,16,18)]]
#save(rasterEWS, rasterSensit, lsTimeSeries2012, file = "rasterEwsSensFabrega2012Modis.rda")

par(mfrow=c(3,1))
plot(rasterEWS[["ar1"]], main = "EWS = Autoregressive Coefficient [1]
     Fabrega", zlim = c(-1, 1)) # ar1
plot(rasterSensit[["ar1"]], main = "Kendall's tau", zlim = c(-1, 1)) # ar1
plot(ndviCatTS2012[296], main = "NDVI 2012", zlim = c(300, 1000)) # ndvi 2012
dev.copy(pdf, 'Fabrega2012ar1.pdf')
dev.off()

par(mfrow=c(3,1))
plot(rasterEWS[["sd"]], main = "EWS = StDev
     Fabrega", zlim = c(0, 300)) # sd
plot(rasterSensit[["sd"]], main = "Kendall's tau", zlim = c(-1, 1)) # ar1
plot(ndviCatTS2012[["X2012"]], main = "NDVI 2012", zlim = c(300, 1000)) # ndvi 2012
dev.copy(pdf, 'Fabrega2012stdev.pdf')
dev.off()

par(mfrow=c(3,1))
plot(rasterEWS[["sk"]], main = "EWS = Skewness
     Fabrega", zlim = c(0, 3)) # sk
plot(rasterSensit[["sk"]], main = "Kendall's tau", zlim = c(-1, 1)) # ar1
plot(ndviCatTS2012[["X2012"]], main = "NDVI 2012", zlim = c(300, 1000)) # ndvi 2012
dev.copy(pdf, 'Fabrega2012sk.pdf')
dev.off()

par(mfrow=c(3,1))
plot(rasterEWS[["kurt"]], main = "EWS = Kurtosis
     Fabrega", zlim = c(0, 12)) # kurt
plot(rasterSensit[["kurt"]], main = "Kendall's tau", zlim = c(-1, 1)) # ar1
plot(ndviCatTS2012[["X2012"]], main = "NDVI 2012", zlim = c(300, 1000)) # ndvi 2012
dev.copy(pdf, 'Fabrega2012kurtosis.pdf')
dev.off()

par(mfrow=c(3,1))
plot(rasterEWS[["cv"]], main = "EWS = Coefficient of Variation
     Fabrega", zlim = c(0, 0.5)) # cv
plot(rasterSensit[["cv"]], main = "Kendall's tau", zlim = c(-1, 1)) # ar1
plot(ndviCatTS2012[["X2012"]], main = "NDVI 2012", zlim = c(300, 1000)) # ndvi 2012
dev.copy(pdf, 'Fabrega2012cv.pdf')
dev.off()

par(mfrow=c(3,1))
plot(rasterEWS[["returnrate"]], main = "EWS = Return Rate
     Fabrega", zlim = c(0, 1)) # returnrate
plot(rasterSensit[["returnrate"]], main = "Kendall's tau", zlim = c(-1, 1)) # ar1
plot(ndviCatTS2012[["X2012"]], main = "NDVI 2012", zlim = c(300, 1000)) # ndvi 2012
dev.copy(pdf, 'Fabrega2012returnrate.pdf')
dev.off()

par(mfrow=c(3,1))
plot(rasterEWS[["densratio"]], main = "EWS = Density Ratio
    Fabrega", zlim = c(0, 3)) # densratio
plot(rasterSensit[["densratio"]], main = "Kendall's tau", zlim = c(-1, 1)) # ar1
plot(ndviCatTS2012[["X2012"]], main = "NDVI 2012", zlim = c(300, 1000)) # ndvi 2012
dev.copy(pdf, 'Fabrega2012densratio.pdf')
dev.off()

par(mfrow=c(3,1))
plot(rasterEWS[["acf1"]], main = "EWS = Autocorrelation Function [1]
     Fabrega", zlim = c(-1, 1)) # acf1
plot(rasterSensit[["acf1"]], main = "Kendall's tau", zlim = c(-1, 1)) # ar1
plot(ndviCatTS2012[["X2012"]], main = "NDVI 2012", zlim = c(300, 1000)) # ndvi 2012
dev.copy(pdf, 'Fabrega2012acf1.pdf')
dev.off()

# note on the side: ar1 & acf1 look promising
