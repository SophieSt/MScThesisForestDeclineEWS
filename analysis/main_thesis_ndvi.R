## main file for thesis on remotely sensed EWS for forest decline

# if not available install the following package
requiredPackages = c('rgdal', 'raster', 'zoo', 'earlywarnings', 'googledrive', 'parallel')

for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, repos = 'http://cran.uni-muenster.de/')
  library(p,character.only = TRUE)
}


setwd("/home/ubuntu/MSc_Proposal_Report/analysis")

source('R/rasterEWS.R')
source('R/rasterSpEWSlist.R')
source('R/removeoutliers.R')
source('R/generic_ews.R')


#load('.RData')

# download file from drive if not there yet
if (!file.exists('data/CatNdvi2012Modis.tif')){
  drive_auth("utils/ttt.rds")
  dl <- drive_download(
    as_id("12YSy2mYPnv3YExzIRoz-51SxheCUyQb3"), path = 'data/CatNdvi2012Modis.tif', overwrite = TRUE)
}

ndviCatTS2012 <- brick("data/CatNdvi2012Modis.tif")
ndviCatTS2012 <- calc(ndviCatTS2012, fun = function(x) {x[x <= 0] <- NA; return(x)})
ndviNA <- calc(ndviCatTS2012, fun = removeoutliers)


# get Catalonia borders
spain <- raster::getData('GADM', path = 'data', country='ESP', level=1)

catalonia <- spain[6,] # referring to province of Catalonia. indexing by name not possible because unix unable to understand spanish enye letter
# reproject catalonia shape (this is the smaller one, should go faster)
catCRS <- spTransform(catalonia, CRS(proj4string(ndviCatTS2012)))

# write Catalonia into shape file
#writeOGR(obj=catCRS, dsn="M://thesis_MSc/data", layer="catBorders", driver="ESRI Shapefile")

# mask out catalonia
catExt <- crop(ndviCatTS2012, catCRS) # crop to extent first for speed reasons
ndvi2012 <- mask(catExt, catCRS)

# rasterize EWS outcome
# rasteredEWSndvi <- rasterEWS(ndvi2012, timeStepName = 'day16x', seriesName = 'ndvi', doInterpolate = TRUE, hasNames = FALSE,
#                              detrending = 'gaussian', bandwidth = 4,
#                              winsize = 50)
# #ewsR <- rasteredEWSndvi[[c(1,3,5,7,9,11,1,3,15,17)]]
# #tauR <- rasteredEWSndvi[[c(2,4,6,8,10,12,14,16,18)]]
# writeRaster(x = rasteredEWSndvi, filename = 'output/rasteredEWSndvi2012.tif', datatype = 'FLT4S', overwrite = TRUE)

# rasterize Spatial EWS
rasteredSpEWSndvi <- rasterSpEWS(ndvi2012, fileout = 'output/ndvi2012', ncores = 5)
#writeRaster(x = rasteredSpEWSndvi, filename = 'output/rasteredSpEWSndvi2012.tif', datatype = 'FLT4S', overwrite = TRUE)


# compare EWS for different species
# ndviEWS <- brick('rasteredEWSndvi2012.tif')
# declines2012 <- readOGR('M:/thesis_MSc/data/data/DBC12_especie.shp')
# declines2012CRS <- spTransform(declines2012, CRS(proj4string(ndviEWS)))
# perpoly <- extract(ndviEWS, declines2012CRS, fun = mean, na.rm = T, df = T, layer = 3, nl = 16)
# perpolyspec <- cbind(perpoly, declines2012CRS$DBC12_espe)
# perpolyspec <- perpolyspec[complete.cases(perpolyspec[,2:16]),2:18]
# perspecies <- aggregate(perpolyspec, by=list(perpolyspec$`declines2012CRS$DBC12_espe`), FUN=mean)[1:17]
# write.csv(perspecies, file = "ndvi2012ewsPerSpecies.csv")
