# setwd("~/Dropbox/git/id_spatial_sim/src/R")
require("raster")
require("rgdal")
landscanfile="/Users/sriley/Dropbox/dataplain/gis/landscan/africa03/africa03/w001001.adf" 
r <- raster(landscanfile)
studyArea1 <- extent(-17.021052,-7.100397,4.519773,12.906347)
studyArea2 <- extent(-14.113223,-13.113467,9.362858,9.947714)
pop1 <- crop(r,studyArea1)
pop2 <- crop(r,studyArea2)
writeRaster(pop1,filename="../../data/westAfricaAscii.txt",format="ascii",overwrite=TRUE)
writeRaster(pop2,filename="../../data/conakryAscii.txt",format="ascii",overwrite=TRUE)


 
