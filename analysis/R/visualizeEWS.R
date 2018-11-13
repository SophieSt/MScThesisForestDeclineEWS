## nice plots
requiredPackages = c('prettymapr', 'latticeExtra', 'colorspace')

for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, repos = 'http://cran.uni-muenster.de/')
  library(p,character.only = TRUE)
}

vegind = 'NDVI'
ewslist <- c('acf1', 'ar1', 'densrat', 'sd', 'sk', 'kurt', 'SpSkew', 'SpVar', 'SpAutoCor')
spewslist <- c('SpSkew', 'SpVar', 'SpAutoCor')

spain <- raster::getData('GADM', path = 'data', country='ESP', level=1)

catalonia <- spain[6,]
catCRS <- spTransform(catalonia, CRS(proj4string(ndviEWS)))

declines2012 <- readOGR('M:/thesis_MSc/data/data/declines2012.kml', layer = 'declines2012')
declines <- spTransform(declines2012, CRS(proj4string(ndviEWS)))


ewsplot <- function(ewsraster, ews, vegind, catCRS) {
  ewstau <- paste0(ews, 'Tau')
  
  png(paste0('M://thesis_MSc/repo_thesis/figures/', vegind, ewstau, '.png'), width = 850, height = 800)
  plot(ewsraster, ewstau, col = diverge_hcl(16), main = '', axes = FALSE, maxpixels = ncell(ndvi1km))
  plot(catCRS, add = TRUE, border = 'grey')
  #plot(declines, add = TRUE, border = 'grey27', lwd = 0.05)
  prettymapr::addnortharrow(scale = 0.6, text.col = 'black', cols = c('black', 'black'))
  title(paste0(vegind, '-based trend in ', ews))
  dev.off()
}

for (ews in ewslist){
  ewsplot(ndvi1km, ews, 'NDVI', catCRS, declines)
}
for (ews in spewslist){
  ewsplot(ndviSpEWS1km, ews, 'NDVI', catCRS, declines)
}
for (ews in ewslist){
  ewsplot(ndmi1km, ews, 'NDMI', catCRS, declines)
}
for (ews in spewslist){
  ewsplot(ndmiSpEWS1km, ews, 'NDMI', catCRS, declines)
}
for (ews in ewslist){
  ewsplot(evi1km, ews, 'EVI', catCRS, declines)
}
for (ews in spewslist){
  ewsplot(eviSpEWS1km, ews, 'EVI', catCRS, declines)
}

ndvitau <- subset(ndvi_preds, subset = c(2,4,6,8,10,12,13,14,15))
scale <- list("SpatialPolygonsRescale", layout.scale.bar(), 
             offset = c(400000, 4500000), scale = 110000, fill=c("transparent","black"), which = 9)
text1 = list("sp.text", c(400000, 4518000), "0", which = 9, cex = 0.6)
text2 = list("sp.text", c(490000, 4518000), "100 km", which = 9, cex = 0.6)
arrow <- list("SpatialPolygonsRescale", layout.north.arrow(), 
              offset = c(485000, 4530000), scale = 35000, which = 9)
cat_bound <- list("sp.polygons", catCRS)
spplot(ndvitau, col.regions = rev(diverge_hcl(n = 100)), zlim = c(-0.55, 0.75), 
       sp.layout = list(scale, text1, text2, arrow, cat_bound), 
       key.space = list(x=0.1,y=.95,corner=c(0,1)))

