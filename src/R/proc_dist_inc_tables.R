# Start from a clean environment
rm(list=ls(all=TRUE))

# Set you present working directory to find the output
# setwd("~/Dropbox/git/id_spatial_sim/src/R")
# setwd("/home/sriley/git/id_spatial_sim/src/R")
# setwd("~/Dropbox/git/id_spatial_sim/src/R")

# Load required packages
require("raster")
require("sp")
require("scales")
require("rgdal")
require("data.table")
require("adegenet")
require("NMOF")

source("idSimFuncs.R")

# Main parameters
# Change values here mainly to get the script to run
weeksUsed <- 0:40
tabledir <- "~/srileytmp/event_files/20141027/"
paramfname <- "paramscan27Oct2014.txt"
outputfile <- paste("~/srileytmp/postProc_",gsub(" ","_",gsub(":","_",date())),".Rdata",sep="")
cat("Output file will be",outputfile,"\n")

# Load up the real data
data.to.use <- paste("~/Dropbox/shares/neil_Ebola/",readLines("~/Dropbox/shares/neil_Ebola/Data_to_use.txt"),".RData",sep="")
load(data.to.use, verbose=TRUE)
dat_full <- dat
dat <- dat_full[,c("EpiCaseDef","DateOnsetInferred","district")]
dat[, "DateOnsetInferred"] <- as.Date(dat[, "DateOnsetInferred"])
dat$district <- as.character(dat$district)
dat <- spatial.prune.v2(dat,0,c(1,2,3),useEpiWeek=FALSE)

# Load the shape files and make a utility df
shapeDir <- "/Users/sriley/srileytmp/sfs"
dists <- readOGR(dsn=shapeDir,layer="ThreeCountries")
countryLatOrder <- c("LIBERIA","SIERRA LEONE","GUINEA")
distsDF <- data.frame(
    dist=as.character(dists$ADM2_NAME),
    country=as.character(dists$ADM0_NAME),
    country_order=match(dists$ADM0_NAME,countryLatOrder),
    lats=dists$CENTER_LAT,
    stringsAsFactors=FALSE)
distsDF <- distsDF[order(distsDF$country_order,distsDF$lats),]

setConakry <- c("KALOUM","DIXINN","MATAM","MATOTO","RATOMA")
setFreetown <- c("WESTERN RURAL","WESTERN")

# Temporary fix for the shape files
distsDF$dist[distsDF$dist=="RATOMA"] <- "CONAKRY"
distsDF$dist[distsDF$dist=="WESTERN"] <- "FREETOWN"
distsDF <- distsDF[-match(c("KALOUM","DIXINN","MATAM","MATOTO","WESTERN RURAL"),distsDF$dist),]
dat$district <- ifelse(dat$district %in% setConakry,"CONAKRY",dat$district)
dat$district <- ifelse(dat$district %in% setFreetown,"FREETOWN",dat$district)

# Need to CHECK districts agree here
tableDistricts <- read.csv(paste(tabledir,"../DistrictLookup.csv",sep=""),stringsAsFactors=FALSE)

# Reorder the districts by country, then latitude
distsLatOrder <- distsDF$dist
distsCountryLO <- distsDF$country

y <- (
      make.incidence.from.linelist(
          distsLatOrder,dat$district,dat$EpiWeek,DTs=weeksUsed
      )
      )$inctab

# Load the parameter table
paramTab <- read.table(file=paste(tabledir,paramfname,sep=""),header=TRUE,sep=" ")
names(paramTab) <- c("R0_Spatial","Offset_Transmit_Spatial","Decay_Transmit_Spatial")

paramsToLoad <- 1:400
realsToLoad <- 0:399

arrAllInc <- array(
    dim=c(dim(y)[1],dim(y)[2],length(realsToLoad),length(paramsToLoad)),
    dimnames=list(1:(dim(y)[1]),colnames(y),realsToLoad,paramsToLoad)
)

# Read in the files using known numbers of params and realizations
system.time(
for (p in 1:length(paramsToLoad)) {
  for (r in 1:length(realsToLoad)) {
    pnum <- paramsToLoad[p]
    rnum <- realsToLoad[r]
    sTT <- fread(paste(tabledir,
        "WestAfrica_paramset",pnum,".",rnum,".adunit.csv",sep=""),
        select=64:125,header=FALSE)
    setnames(sTT,names(sTT),tableDistricts[,"WHO_districts_v2"])
    tmp <- convert.spatialsim.WHO.V1(sTT,colnames(y),tableDistricts,length(weeksUsed))
    arrAllInc[,,r,p] <- tmp
  }
  cat("Completed",p,"of",length(paramsToLoad),"\n")
}
)

save(arrAllInc,y,paramsToLoad,paramTab,distsCountryLO,file=outputfile)
