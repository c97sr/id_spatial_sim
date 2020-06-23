#' # Setup possible networks for covid transmission studies
#'
#' This script is designed to run pretty quickly for a small population.
#' 
#' First clear the environment of variables
rm(list=ls(all=TRUE))

#' Then load up stable packages we might need
require("raster")
require("scales")
require("sp")
require("rgdal")
require("devtools")

#' And SR's sandbox package if needed
## Install with install_github("c97sr/idd") if needed
require("idd")

#' Set top-level dir and file stems for this script. Assuming it
#' will be run from the top dir of the scenario.
dirTop <- "../../"
fnStemNetworks <- "ncov_small_v2"
fnStemRuns <- paste(fnStemNetworks,"_run",sep="")

#' Source the R functions used here. Should be in a package, that can be
#' installed using github_install with subdir option.
source(paste(dirTop,"src/rcode/idSimFuncs.R",sep=""))

#' ## Build synthetic population
#' 
#' This call to the executable should make a single network of 142000 nodes
#' each assigned to a wrokplace of 300 people at which each node makes 30
#' connections. The commuting behaviour is assumed to be similar to people
#' going to work and the propoensity for long journeys is far too high
#' for schools, especially primary
system(
    paste(
        paste(dirTop,"/src/cpp/ebola_build.exe",sep=""),
        "./params/b_params.in",
        paste("./output/",fnStemNetworks,sep="")
    )
)

#' This chart is the initial random distribution of distances to the workplace
#' (school), the distribution after 67 million MCMC updates and the
#' distribution after 100 million updates. The latter 2 are so similar you
#' cannot see the 67 million update line.
seqn <- c(0,6,9)
x1 <- read.csv(paste(
    "./output/",fnStemNetworks,"_commute_dist_",seqn[1],".csv",sep=""))
x2 <- read.csv(paste(
    "./output/",fnStemNetworks,"_commute_dist_",seqn[2],".csv",sep=""))
x3 <- read.csv(paste(
    "./output/",fnStemNetworks,"_commute_dist_",seqn[3],".csv",sep=""))
plot(x1,xlim=c(0,40),ylim=c(0,15000),type="l",lwd=2,col="red")
points(x2,type="l",lwd=2,col="green")
points(x3,type="l",lwd=2,col="blue")

#' Run an epidemic on this population. This is not used for anthing other
#' than to make sure the network build code revisions have not broken the
#' simulation code in an obvious way.
system(
    paste(
        paste(dirTop,"/src/cpp/ebola_run.exe",sep=""),
        "./params/r_params.in",
        paste("./output/",fnStemNetworks,sep=""),
        paste("./output/",fnStemRuns,sep="")
    )
)
