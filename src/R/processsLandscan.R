# setwd("~/Dropbox/git/id_spatial_sim/src/R")
require("raster")
require("rgdal")
landscanfile="/Users/sriley/Dropbox/dataplain/gis/landscan/africa03/africa03/w001001.adf" 
r <- raster(landscanfile)

extWide <- extent(-17.021052,-7.100397,4.519773,12.906347)
extConackry <- extent(-14.113223,-13.113467,9.362858,9.947714)
extMonrovia <- extent(-10.948165,-10.551283,6.183227,6.508067)

sum(as.matrix(crop(r,extWide)),na.rm=TRUE)

writeRaster(crop(r,extWide),filename="../../data/westAfricaAscii.txt",format="ascii",overwrite=TRUE)
writeRaster(crop(r,extConackry),filename="../../data/conakryAscii.txt",format="ascii",overwrite=TRUE)
writeRaster(crop(r,extMonrovia),filename="../../data/monroviaAscii.txt",format="ascii",overwrite=TRUE)
