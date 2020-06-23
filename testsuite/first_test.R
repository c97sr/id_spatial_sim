#' # Fast test for the infectious disease spatial simulation library
#' 
#' This needs fleshing out, but it should be independent of actual scenarios
#' as much as possible. As this is currently designed, it will only run on 
#' linux-like machines. More specifically, it won't run on windows.
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

#' And finally load up functions held elsewhere in this repo
source("../src/rcode/idSimFuncs.R")

#' ## Build synthetic population
#' 
#' The first line of the batch file builds a synthetic population with density
#' proportional the ebola affected region in west Africa, but much smaller. 
#' With a total population of only 100,000. Each person has, one average, 10 network 
#' links but these links are distributed entirely randomly in space. This takes a 
#' while because the average population is very low and there is high variability. 
#' Hence the accept-reject method for assinging nodes has many rejection steps. We assume
#' that only one individual lives in a household for this population. 
system(paste(	"../src/cpp/ebola_build.exe",
				"./params/fast_test_build_params.in",
				"./output/pop1"))

#' ## Run an epidemic on the population
#' 
#' The second line of the batch file runs an outbreak of only two generations 20
#' times. The outbreak is seeded in the same area as the reported patient zero 
#' for the 2014 Ebola outbreak. There are 4 iniitally infectious individuals 
#' at time $t=0$. Transmission is only via the spatial kernel and thus allows 
#' us to test that the basic reproductive number is parameterized correctly. 
#' We can also report the serial interval. 
system(paste(	"../src/cpp/ebola_run.exe",
				"./params/fast_test_run_params.in",
				"./output/pop1",
				"./output/pop1_test"))

#' We first load the linelist of events from all the realizations. And check the
#' dimensions of the output. The output was designed before csvs became so
#' dominant!
dat0 <- read.table(file="./output/pop1_test_pset_2_Events.out",header=TRUE)
dimDat0 <- dim(dat0)
noevents <- dimDat0[1]
nocols <- dimDat0[2]

#' The column headings describe the information captured in the event file
names(dat0)

#' We can subset these 'data' to look at only infection. Then examine the number 
#' of infections by generation for each realization.
tabInfs0 <- dat0[dat0$Event==0,]
table(tabInfs0$Run,tabInfs0$Generation)

