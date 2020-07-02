library(rgdal)
library(raster)
library(ggplot2)
library(sp)
library(rgeos)
library(maptools)

setwd("/home/james/hadean/id_spatial_sim/scenarios/covid_uk/")


e <- matrix(nrow = 2, ncol = 2)
e[1,1] = -10.8
e[1,2] = 2
e[2,1] = 49.1
e[2,2] = 60.9

c2 <- read.asciigrid("Densities/gpw_v4_population_density_rev11_2020_30_sec_2.asc")
summary(c2)
c2 <- as(crop(raster(c2), extent(e)), 'SpatialGridDataFrame')
summary(c2)
plot(c2)

c3 <- read.asciigrid("Densities/gpw_v4_population_density_rev11_2020_30_sec_3.asc")
summary(c3)
c3 <- as(crop(raster(c3), extent(e)), 'SpatialGridDataFrame')
summary(c3)
plot(c3)

c <- as(merge(raster(c2), raster(c3)), 'SpatialGridDataFrame')
plot(c)

m <- matrix(nrow = 6, ncol = 2)
m[1,1] = e[1,1]
m[1,2] = e[2,2]

m[2,1] = e[1,1]
m[2,2] = e[2,1]

m[3,1] = -4
m[3,2] = 50

m[4,1] = 1
m[4,2] = 50.6


m[5,1] = e[1,2]
m[5,2] = 51.1


m[6,1] = e[1,2]
m[6,2] = e[2,2]

p <- Polygon(m)
ps <- Polygons(c(p), "j")
sp <- SpatialPolygons(c(ps))


grid <- c[!is.na(over(c, sp)),]

plot(grid)
summary(grid)

writeRaster(raster(grid), "data/UK_Ireland_2020_30_sec.asc", overwrite = TRUE)


