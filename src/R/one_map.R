rm(list=ls(all=TRUE))

# Set you present working directory to find the output
# setwd("~/Dropbox/git/id_spatial_sim/src/R")

# Load required packages
require("raster")
require("sp")
require("scales")
require("rgdal")
source("~/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")
source("../../../sr_source/ebola/ebolaFuncs.R")

weeksUsed <- 0:30

# Load up the real data
data.to.use <- paste("~/Dropbox/shares/neil_Ebola/",readLines("~/Dropbox/shares/neil_Ebola/Data_to_use.txt"),".RData",sep="")
load(data.to.use, verbose=TRUE)
dat_full <- dat
dat <- dat_full[,c("EpiCaseDef","DateOnsetInferred","district")]
dat[, "DateOnsetInferred"] <- as.Date(dat[, "DateOnsetInferred"])
dat$district <- as.character(dat$district)
dat <- spatial.prune.v2(dat,0,c(1,2,3),useEpiWeek=FALSE)

# Trim the data down to the first 1000
dat <- dat[order(dat$EpiDay)[1:1000],]

# Load the shape files
shapeDir <- "/Users/sriley/srileytmp/sfs"
dists <- readOGR(dsn=shapeDir,layer="ThreeCountries")

y <- (
  make.incidence.from.linelist(
      as.character(dists$ADM2_NAME),dat$district,dat$EpiWeek,DTs=weeksUsed
  )
)$inctab

inc.heat.chart.pdf.v2(y,vecCountries=as.character(dists$ADM0_NAME),outstem=paste("~/srileytmp/data_heat_chart",sep=""))

fnPopdata = "../../data/westAfricaAscii.asc"
fnOutStem = "~/Dropbox/projects/ebola/"

# Make an aggregated version of the rater file
popgrid <- read.asciigrid(fnPopdata,as.image=TRUE)
sum(popgrid$z,na.rm=TRUE)
# popgrid_tmp <- raster(popgrid)
# popgrid_coarse_tmp <- aggregate(popgrid_tmp,factor=1000,fun=sum)
# writeRaster(popgrid_coarse_tmp,file="../data/westAfricaAscii_agg1000.asc",format="ascii",
#    overwrite=TRUE)

# popgrid <- read.asciigrid(fnPopdata,as.image=TRUE)
sum(popgrid$z,na.rm=TRUE)

# Preconditions for the batch runs, remember weeksUsed
allBatches <- 1:9
allReals <- 0:99
maxoff <- 5
noff <- 2*maxoff+1
arrSumSq <- array(dim=c(noff,length(allReals),length(allBatches)))
arrAllInc <- array(
    dim=c(dim(y)[1],dim(y)[2],noff,length(allReals),length(allBatches)),
    dimnames=list(rownames(y),colnames(y),1:noff,allReals,allBatches)
)

for (j in allBatches) {
  fnEpidata = paste("~/srileytmp/event_files/20141007/batch_",j,"_pset_0_Events.out",sep="")
  datSim <- read.table(file=fnEpidata,header=TRUE)
  datSim <- datSim[datSim$Event==0,]
  datSim$X <- datSim$X * 180 / pi
  datSim$Y <- datSim$Y * 180 / pi
  datSim <- cbind(datSim,dist_code=overlay(SpatialPoints(cbind(lon=datSim$X,lat=datSim$Y)),dists))
  datSim <- cbind(datSim,district=dists$ADM2_NAME[datSim$dist_code])
  datSim$EpiWeek <- ceiling(datSim$Day/7.0)
  
  # Up to here. Need to make this bit work
  noReals <- max(datSim$Run)
  for (i in 0:noReals) {
    dat_one_r <- datSim[datSim$Run==i,]
    tmp <- make.incidence.from.linelist(
        as.character(dists$ADM2_NAME),
        dat_one_r$district,
        dat_one_r$EpiWeek,
        DTs=((min(weeksUsed)-maxoff)):(max(weeksUsed)+maxoff))
    for (k in 1:(noff)) {
      # XXXX up to here
      x <- tmp$inctab[k:(k+length(weeksUsed)-1),]
      arrAllInc[,,k,i+1,j] <- x
      
      # Need to run a sliding window for 
      arrSumSq[k,i+1,j] <- sum((x-y)^2)    
      
    }
    
  }
  cat("batch",j,"complete\n")
}

log10(min(arrSumSq))
ibest <- which(arrSumSq==min(arrSumSq),arr.ind=TRUE) 

x <- y
x[] <- arrAllInc[,,ibest[1],ibest[2],ibest[3]]
inc.heat.chart.pdf.v2(
    x,
    vecCountries=as.character(dists$ADM0_NAME),
    outstem=paste("~/srileytmp/best",sep="")
)
sum(x)

for (b in allBatches) {
  x[] <- arrAllInc[,,index_mins[b],b] 
  cat(sum(x),"\n")
  inc.heat.chart.pdf.v2(
      x,
      vecCountries=as.character(dists$ADM0_NAME),
      outstem=paste("~/srileytmp/best_batch_",b,sep="")
  )  
}

# epiImage <- eventImage(datSim,popgrid,0,0,1000,0,0)
# sum(epiImage$z,na.rm=TRUE)

# png(file=paste(fnOutStem,"figs/map.png",sep=""),width=12,height=12,units="cm",res=300)
#  image(popgrid$x,popgrid$y,log(popgrid$z+0.5),col=rev(grey_pal()(100)))
#  image(epiImage$x,epiImage$y,epiImage$z,add=TRUE)
# dev.off()
