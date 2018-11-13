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
if (!file.exists('data/CatNdmi2012Modis.tif')){
  drive_auth("utils/ttt.rds")
  dl <- drive_download(
    as_id("1Mky5WdGFDsWj8srsbFWAdu8rpRg6DbQY"), path = 'data/CatNdmi2012Modis.tif', overwrite = TRUE)
}

ndmiCatTS2012 <- brick("data/CatNdmi2012Modis.tif")
ndmiCatTS2012 <- calc(ndmiCatTS2012, fun = function(x) {x[x <= 0] <- NA; return(x)})

# get Catalonia borders
spain <- raster::getData('GADM', path = 'data', country='ESP', level=1)

catalonia <- spain[6,]
# reproject catalonia shape (this is the smaller one, should go faster)
catCRS <- spTransform(catalonia, CRS(proj4string(ndmiCatTS2012)))

# write Catalonia into shape file
#writeOGR(obj=catCRS, dsn="M://thesis_MSc/data", layer="catBorders", driver="ESRI Shapefile")

# mask out catalonia
catExt <- crop(ndmiCatTS2012, catCRS)
ndmi2012 <- mask(catExt, catCRS)

# rasterize EWS outcome
# rasteredEWSndmi <- rasterEWS(ndmi2012, timeStepName = 'day16x', seriesName = 'ndmi', doInterpolate = TRUE, hasNames = FALSE,
#                              detrending = 'gaussian', bandwidth = 4,
#                              winsize = 50)
# #ewsR <- rasteredEWSndmi[[c(1,3,5,7,9,11,1,3,15,17)]]
# #tauR <- rasteredEWSndmi[[c(2,4,6,8,10,12,14,16,18)]]
# writeRaster(x = rasteredEWSndmi, filename = 'output/rasteredEWSndmi2012.tif', datatype = 'FLT4S', overwrite = TRUE)

# rasterize Spatial EWS
rasteredSpEWSndmi <- rasterSpEWS(ndmi2012, fileout = 'output/ndmi2012', ncores = 5)
#writeRaster(x = rasteredSpEWSndmi, filename = 'output/rasteredSpEWSndmi2012.tif', datatype = 'FLT4S', overwrite = TRUE)



# ndmiEWS <- brick('rasteredEWSndmi2012.tif')
# declines2012 <- readOGR('M:/thesis_MSc/data/data/DBC12_especie.shp')
# declines2012CRS <- spTransform(declines2012, CRS(proj4string(ndmiEWS)))
# perpoly <- extract(ndmiEWS, declines2012CRS, fun = mean, na.rm = T, df = T, layer = 3, nl = 16)
# perpolyspec <- cbind(perpoly, declines2012CRS$DBC12_espe)
# perpolyspec <- perpolyspec[complete.cases(perpolyspec[,2:16]),2:18]
# perspecies <- aggregate(perpolyspec, by=list(perpolyspec$`declines2012CRS$DBC12_espe`), FUN=mean)[1:17]
# write.csv(perspecies, file = "ndmi2012ewsPerSpecies.csv")
