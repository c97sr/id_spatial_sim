# Start from a clean environment
rm(list=ls(all=TRUE))

# Set you present working directory to find the output
# Just to check the Dropbox solution works
# And it seems to work on the nexus as well
# setwd("~/Dropbox/git/id_spatial_sim/src/R")
# setwd("/home/sriley/git/id_spatial_sim/src/R")

# Load required packages
# require("raster")
# require("sp")
# require("scales")
# require("rgdal")
# require("data.table")
require("adegenet")
# qrequire("NMOF")

source("idSimFuncs.R")

reportdir <- "~/Dropbox/projects/ebola/sim_reports/report_0.91/"

load("~/srileytmp/postProc_Mon_Oct_27_14_04_57_2014.RData",verbose=TRUE)

# Below here to become a separate file pretty soon
noweeks <- dim(arrAllInc)[1]
nodists <- dim(arrAllInc)[2]
noreals <- dim(arrAllInc)[3]
noparams <- dim(arrAllInc)[4]
testedStats <- c("sumSq","sumSqReport","reportLevel","largeSumSq")
compStats <- array(
    dim=c(length(testedStats),dim(arrAllInc)[3],dim(arrAllInc)[4]),
    dimnames=list(testedStats,1:dim(arrAllInc)[3],paramsToLoad)
)

for (r in 1:noreals) {
  for (p in 1:noparams) {
    
    # Simple sun of squares
    x <- arrAllInc[,,r,p]
    compStats["sumSq",r,p] <- sum((x-y)^2)
    reportingCorr <- max(1,sum(x)/sum(y))
    compStats["sumSqReport",r,p] <- sum((x-y*reportingCorr)^2)
    compStats["reportLevel",r,p] <- reportingCorr
    if (reportingCorr > 1) {
      compStats["largeSumSq",r,p] <- compStats["sumSqReport",r,p]
    } else {
      compStats["largeSumSq",r,p] <- 1e100
    }
    
  }
}

# Playing here with particle mcmc in discrete parameter space
modLike <- function(y,x,pen=0) {
  nrow <- dim(y)[1]
  ncol <- dim(y)[2]
  r <- 1
  c <- 1
  rtn <- 0
  while (r < nrow) {
    while (c < ncol) {
      if (y[r,c] > 0 && x[r,c] < 1) {
        rtn <- rtn + pen 
      } else {
        rtn <- rtn + dpois(y[r,c],x[r,c],log=TRUE)
      }
      c <- c + 1
    }
    r <- r + 1
  }
  rtn
}

modLike(y,arrAllInc[,,100,10])
