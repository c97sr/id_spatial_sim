rm(list=ls(all=TRUE))
# setwd("~/Dropbox/git/id_spatial_sim/testsuite")
require("raster")
require("scales")
require("sp")
require("rgdal")
source("~/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")

dat0 <- read.table(file="./output/ft_sp_pset_0_Events.out",header=TRUE)
dimDat0 <- dim(dat0)
noevents <- dimDat0[1]
nocols <- dimDat0[2]
names(dat0)
tabInfs0 <- dat0[dat0$Event==0,]
table(tabInfs0$Run,tabInfs0$Generation)
table(tabInfs0$Generation)[2]/table(tabInfs0$Generation)[1]

dat1 <- read.table(file="./output/ft_sp_pset_1_Events.out",header=TRUE)
dim(dat1)
tabInfs1 <- dat1[dat1$Event==0,]
estR0 <- table(tabInfs1$Generation)[2] / table(tabInfs1$Generation)[1]
estR0
tabRunGen1 <- table(tabInfs1$Run,tabInfs1$Generation)
hist( tabRunGen1[,2]/tabRunGen1[,1],breaks=seq(0,7,0.25),
    xlab="Ratio first to second gens",main="")
abline(v=estR0,lwd=2,col="red")	
vecTimesG1 <- tabInfs1[tabInfs1$Generation==2,"Day"] 
mean(vecTimesG1)
hist(vecTimesG1,breaks=seq(0,50,2),main="",
    xlab="Time from exposure to infection (2 day bins)")


popgrid <- read.asciigrid("../data/westAfricaAscii_agg100.asc",as.image=TRUE)
sum(popgrid$z,na.rm=TRUE)
popgrid$x <- popgrid$x * pi / 180
popgrid$y <- popgrid$y * pi / 180
epiImage <- eventImage(dat1,popgrid,0,0,1000,0,1000)
image(popgrid$x,popgrid$y,log(popgrid$z+0.5),col=rev(grey_pal()(100)),xlab="",ylab="")
image(epiImage$x,epiImage$y,epiImage$z,add=TRUE)

shapeDir <- "/Users/sriley/srileytmp/sfs"
dists <- readOGR(dsn=shapeDir,layer="ThreeCountries")
summary(dists)
plot(dists)

incidence <- overlay(SpatialPoints(cbind(lon=tabInfs1$X*180/pi,lat=tabInfs1$Y*180/pi)),dists)

