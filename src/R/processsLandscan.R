# setwd("~/Dropbox/git/id_spatial_sim/src/R")
require("raster")
require("rgdal")
landscanfile="/Users/sriley/Dropbox/dataplain/gis/landscan/africa03/africa03/w001001.adf" 
r <- raster(landscanfile)
studyArea <- extent(-17.021052,-7.100397,4.519773,12.906347) 
pop <- crop(r,studyArea)
sum(as.matrix(pop),na.rm=TRUE)
sum()
writeRaster(pop,filename="../../data/westAfricaAscii.txt",format="ascii",overwrite=TRUE)

