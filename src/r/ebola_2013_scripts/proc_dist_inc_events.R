rm(list=ls(all=TRUE))

# Set you present working directory to find the output
# setwd("~/Dropbox/git/id_spatial_sim/src/R")
# setwd("~/git/id_spatial_sim/src/R")

# Load required packages
require("raster")
require("sp")
require("scales")
require("rgdal")
# source("~/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")
source("../../../sr_source/ebola/ebolaFuncs.R")
source("idSimFuncs.R")

weeksUsed <- 0:40
maxNoEvents <- 60000

# Load up the real data
data.to.use <- paste("~/Dropbox/shares/neil_Ebola/",readLines("~/Dropbox/shares/neil_Ebola/Data_to_use.txt"),".RData",sep="")
load(data.to.use, verbose=TRUE)
dat_full <- dat
dat <- dat_full[,c("EpiCaseDef","DateOnsetInferred","district")]
dat[, "DateOnsetInferred"] <- as.Date(dat[, "DateOnsetInferred"])
dat$district <- as.character(dat$district)
dat <- spatial.prune.v2(dat,0,c(1,2,3),useEpiWeek=FALSE)
dat <- dat[!is.na(dat$district),]

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

# Preconditions for the batch runs, remember weeksUsed
R0s <- c("1.80","1.60","1.40")
# R0s <- c("1.80")
chosen_params <- 1:100
chosenReals <- sample(0:399,50)
# chosenReals <- 0:399
# chosen_params <- sample(1:100,1)
totalBatches <- length(R0s)*length(chosen_params)
arrAllInc <- array(
    dim=c(dim(y)[1],dim(y)[2],length(chosenReals),totalBatches),
    dimnames=list(rownames(y),colnames(y),1:length(chosenReals),1:totalBatches)
)

for (i in 0:(length(R0s)-1)) {
  for (j in 1:length(chosen_params)) {
    
    tmp <- make.incidence.from.batch(
        y,
        paste("~/srileytmp/events/gemma_20141017/WestAfrica_R",R0s[i+1],"_paramset",chosen_params[j],".infevents.csv",sep=""),
        dists,
        fileformat="EbolaSim",
        weeksUsed,
        distsLatOrder,
        chosenReals,
        maxNoEvents
    )
    
    arrAllInc[,,,i*length(chosen_params)+j] <- tmp$inc
    
    cat("batch",i*length(chosen_params)+j," of ",totalBatches,"complete\n")
    
    gc()
    
  }
}

# Read in the parameter file that goes with the runs
params <- read.table("~/srileytmp/events/gemma_20141017/paramscan.txt",header=TRUE)

# Save the results as a binary
save(arrAllInc,chosen_params,params,file="~/Dropbox/tmp/postProc20oct_opt.Rdata")
