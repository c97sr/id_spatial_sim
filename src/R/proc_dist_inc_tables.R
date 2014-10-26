rm(list=ls(all=TRUE))

# Set you present working directory to find the output
# setwd("~/Dropbox/git/id_spatial_sim/src/R")
# setwd("/home/sriley/git/id_spatial_sim/src/R")

# Load required packages
require("raster")
require("sp")
require("scales")
require("rgdal")
# source("~/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")
# source("../../../sr_source/ebola/ebolaFuncs.R")
source("idSimFuncs.R")

weeksUsed <- 0:40

# Load up the real data
data.to.use <- paste("~/Dropbox/shares/neil_Ebola/",readLines("~/Dropbox/shares/neil_Ebola/Data_to_use.txt"),".RData",sep="")
load(data.to.use, verbose=TRUE)
dat_full <- dat
dat <- dat_full[,c("EpiCaseDef","DateOnsetInferred","district")]
dat[, "DateOnsetInferred"] <- as.Date(dat[, "DateOnsetInferred"])
dat$district <- as.character(dat$district)
dat <- spatial.prune.v2(dat,0,c(1,2,3),useEpiWeek=FALSE)

# Trim the data down to the first 1000
# dat <- dat[order(dat$EpiDay)[1:1000],]

# Load the shape files
shapeDir <- "/home/sriley/srileytmp/sfs/"
dists <- readOGR(dsn=shapeDir,layer="ThreeCountries")

# Reorder the districts by country, then latitude
distsLatOrder <- as.character(dists$ADM2_NAME[rev(order(dists$CENTER_LAT))])
distsCountryLO <- as.character(dists$ADM0_NAME[rev(order(dists$CENTER_LAT))])

y <- (
      make.incidence.from.linelist(
          distsLatOrder,dat$district,dat$EpiWeek,DTs=weeksUsed
      )
      )$inctab

# Start to look at the table files
tabledir <- "/home/sriley/srileytmp/events/gemma_20141021/"
paramfname <- "paramscan.txt"

# Load the parameter table
paramTab <- read.table(file=paste(tabledir,paramfname,sep=""),header=TRUE,sep=" ")
names(paramTab) <- c("R0_Spatial","Decay_Transmit_Spatial","Offset_Transmit_Spatial")

# Load the place lookups
tableDistricts <- read.table("../../data/WestAfricaDistricts.txt")
names(tableDistricts) <- c("Code","Country","Districts")

# Read in a single table
dataDir <- "/home/sriley/srileytmp/events/gemma_20141021"
testTab <- read.csv("/home/sriley/srileytmp/events/gemma_20141021/WestAfrica_paramset111.121.adunit.csv")

save(arrAllInc,chosen_params,params,file="~/srileytmp/events/gemma_20141017/postProc.Rdata")
