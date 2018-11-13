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

# download file
if (!file.exists('data/CatEvi2012Modis.tif')){
  drive_auth("utils/ttt.rds")
  dl <- drive_download(
    as_id("1N4S_-onpG37u37IEhOlsXWQFbNcXGtjG"), path = 'data/CatEvi2012Modis.tif', overwrite = TRUE)
}

eviCatTS2012 <- brick("data/CatEvi2012Modis.tif")
eviCatTS2012 <- calc(eviCatTS2012, fun = function(x) {x[x <= 0] <- NA; return(x)})

# get Catalonia borders
spain <- raster::getData('GADM', path = 'data', country='ESP', level=1)

catalonia <- spain[6,]
# reproject catalonia shape (this is the smaller one, should go faster)
catCRS <- spTransform(catalonia, CRS(proj4string(eviCatTS2012)))

# write Catalonia into shape file
#writeOGR(obj=catCRS, dsn="M://thesis_MSc/data", layer="catBorders", driver="ESRI Shapefile")

# mask out catalonia
catExt <- crop(eviCatTS2012, catCRS)
evi2012 <- mask(catExt, catCRS)

# rasterize EWS outcome
rasteredEWSevi <- rasterEWS(evi2012, timeStepName = 'day16x', seriesName = 'evi', doInterpolate = TRUE, hasNames = FALSE,
                             detrending = 'gaussian', bandwidth = 4,
                             winsize = 50)
writeRaster(x = rasteredEWSevi, filename = 'output/rasteredEWSevi2012.tif', datatype = 'FLT4S', overwrite = TRUE)

# rasterize Spatial EWS
rasteredSpEWSevi <- rasterSpEWS(evi2012, fileout = 'output/evi2012', ncores = 5)

