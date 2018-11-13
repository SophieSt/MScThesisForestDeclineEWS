# load ndvi EWS
requiredPackages = c('rgdal', 'raster', 'zoo', 'earlywarnings', 'googledrive', 'parallel')

for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, repos = 'http://cran.uni-muenster.de/')
  library(p,character.only = TRUE)
}


source('R/removeoutliers.R')

ndvi <- brick('data/CatNdvi2012Modis.tif')
ndvi@file@nodatavalue <- -32768

spain <- raster::getData('GADM', path = 'data', country='ESP', level=1)
catalonia <- spain[6,] # referring to province of Catalonia. indexing by name not possible because unix unable to understand spanish enye letter
catCRS <- spTransform(catalonia, CRS(proj4string(ndvi)))
treecover <- raster('data/treeCover2011Modis.tif')
treecover@file@nodatavalue <- 0
forestcover <- crop(treecover, extent(catCRS)) %>%
  mask(catCRS)

ndviNA <- ndvi %>%
  calc(fun = removeoutliers) %>%
  calc(fun = is.na) 
save(ndviNA, file = 'ndviNA.RData')

timestep_na <- cellStats(ndviNA, stat = 'sum', na.rm = TRUE)
timestep_na_perc <- (timestep_na / ncell(ndvi)) * 100
barplot(timestep_na_perc, xaxt = 'n', xlab = '', ylab = '', cex.axis = 1.05) 
axis(side = 1, at = c(0, 23, 46, 69, 92, 115, 138, 161, 184, 207, 230, 253, 272) * 1.2,
     labels = c(2000:2012), cex.axis = 1.05, line = 0.2)
title(main = '', xlab = 'year', ylab = "No of NA's per time step [%]", line = 2.3, font.lab = 2)

save(timestep_na, timestep_na_perc, file = 'ndvi_na_timestep.RData')




