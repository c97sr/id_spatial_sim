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
source("../../../sr_source/ebola/ebolaFuncs.R")
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
shapeDir <- "/Users/sriley/srileytmp/sfs/"
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
R0s <- c("1.40","1.60","1.8")
chosenReals <- sample(0:399,100)

load("~/Dropbox/tmp/postproc.Rdata",verbose=TRUE)
# save(arrAllInc,chosen_params,params,file="~/srileytmp/events/gemma_20141017/postProc.Rdata")

diff_stats <- c("sum_sq")

sumArr <- array(dim=c(length(diff_stats),dim(arrAllInc)[3],dim(arrAllInc[4])))

totalBatches <- length(R0s)*length(chosen_params)

# Up to here. Need to trawl through for best fit

dim(arrAllInc)

log10(min(arrSumSq))
plot(apply(arrSumSq,c(3),min))
ibest <- which(arrSumSq==min(arrSumSq),arr.ind=TRUE) 

x <- y
x[] <- arrAllInc[,,ibest[1],ibest[2],ibest[3]]
rownames(x) <- rownames(y)
inc.heat.chart.pdf.v2(
    x,
    vecCountries=distsCountryLO,
    outstem=paste("~/srileytmp/best",sep=""),
    xlabs=seq(0,40,5)
)

inc.heat.chart.pdf.v2(
    y,
    vecCountries=distsCountryLO,
    outstem=paste("~/srileytmp/data_heat_chart",sep=""),
    xlabs=seq(0,40,5)
)

