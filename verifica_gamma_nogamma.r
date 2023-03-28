library(raster)
install.packages("rasterVis")
library(rasterVis)
try<-raster("D:/Desktop/prova gamma.tif")
try1<-raster("D:/Desktop/no gamma.tif")

image(try)
levelplot(try, 
          col.regions = rev(terrain.colors(255)),
          at = seq(0, 255, length = 256),
          colorkey = list(space = "bottom", 
                          width = 1, height = 0.5, 
                          labels=list(at=seq(0, 255, by=50))),
          margin = FALSE)
par(mfrow=c(1,2))
