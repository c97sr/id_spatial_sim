# setwd("~/Dropbox/02_git/id_ss/src/R")
require("raster")
require("rgdal")
landscanfile="E:/Dropbox/data/gis/landscan/africa03/africa03/w001001.adf" 
r <- raster(landscanfile)

# as (x,y) in degrees
extWide <- extent(-17.021052,-7.100397,4.519773,12.906347)
extConackry <- extent(-14.113223,-13.113467,9.362858,9.947714)
extMonrovia <- extent(-10.948165,-10.551283,6.183227,6.508067)
extBikoroSmall <- extent(17.3015,19.1747,-1.60866,0.863016)

sum(as.matrix(crop(r,extBikoro)),na.rm=TRUE)

writeRaster(crop(r,extWide),filename="../../data/westAfricaAscii.txt",format="ascii",overwrite=TRUE)
writeRaster(crop(r,extConackry),filename="../../data/conakryAscii.txt",format="ascii",overwrite=TRUE)
writeRaster(crop(r,extMonrovia),filename="../../data/monroviaAscii.txt",format="ascii",overwrite=TRUE)
writeRaster(crop(r,extBikoroSmall),filename="../../data/bikoroSmallAscii.txt",format="ascii",overwrite=TRUE)
