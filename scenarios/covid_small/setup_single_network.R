#' # Setup possible networks for covid transmission studies
#'
#' This script is designed to run pretty quickly for a small population.
#' 
#' First clear the environment of variables
rm(list=ls(all=TRUE))

#' Then load up packages maintained by others
require("raster")
require("scales")
require("sp")
require("rgdal")
require("devtools")

#' And packages maintained by the id_spatial_sim team (or close to)
## Install with install_github("c97sr/idd") if needed
require("idd")

#' Set top-level dir and load currently useful R functions
source("../../src/rcode/idSimFuncs.R")

#' ## Build synthetic population
#' 
#' The first line of the batch file builds a synthetic population with density
#' proportional the ebola affected region in west Africa, but much smaller. 
#' With a total population of only 100,000. Each person has, one average, 10 network 
#' links but these links are distributed entirely randomly in space. This takes a 
#' while because the average population is very low and there is high variability. 
#' Hence the accept-reject method for assinging nodes has many rejection steps. We assume
#' that only one individual lives in a household for this population. 
system(paste(	"../../build/ebola_build.exe",
				"./params/b_params.in",
				"./output/ncov_small"))

#' ## Run an epidemic on the population, to make sure that the thing will run
system(paste(	"../../build/ebola_run.exe",
				"./params/r_params.in",
				"./output/ncov_small",
				"./output/ncov_small_runs"))

#' We first load the linelist of events from all the realizations. And check the
#' dimensions of the output. The output was designed before csvs became so
#' dominant!
dat0 <- read.table(file="./output/pop1_test_pset_0_Events.out",header=TRUE)
dimDat0 <- dim(dat0)
noevents <- dimDat0[1]
nocols <- dimDat0[2]

#' The column headings describe the information captured in the event file
names(dat0)

#' We can subset these 'data' to look at only infection. Then examine the number 
#' of infections by generation for each realization.
tabInfs0 <- dat0[dat0$Event==0,]
table(tabInfs0$Run,tabInfs0$Generation)
