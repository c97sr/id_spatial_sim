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

# Read in the parameter file that goes with the runs
params <- read.table("~/srileytmp/events/gemma_20141017/paramscan.txt",header=TRUE)

# Preconditions for the batch runs, remember weeksUsed
R0s <- c("1.40","1.60","1.8")
chosen_params <- sample(1:100,50)
chosenReals <- sample(0:399,100)
totalBatches <- length(R0s)*length(chosen_params)
arrAllInc <- array(
    dim=c(dim(y)[1],dim(y)[2],length(chosenReals),totalBatches),
    dimnames=list(rownames(y),colnames(y),1:length(chosenReals),1:totalBatches)
)

for (i in 0:(length(R0s)-1)) {
  for (j in 1:length(chosen_params)) {
    
    tmp <- make.incidence.from.batch(
        y,
#   paste("~/srileytmp/event_files/20141007/batch_",j,"_pset_0_Events.out",sep=""),
#   paste("~/srileytmp/event_files/gemma_20141017/WestAfrica_R1.40_paramset",j,".infevents.csv",sep="")      
        paste("~/srileytmp/events/gemma_20141017/WestAfrica_R",R0s[i+1],"_paramset",chosen_params[j],".infevents.csv",sep=""),
        dists,
        fileformat="EbolaSim",
        weeksUsed,
        distsLatOrder,
        chosenReals
    )
    
    gc()
    
    arrAllInc[,,,i*length(chosen_params)+j] <- tmp$inc
    
    cat("batch",i*length(chosen_params)+j," of ",totalBatches,"complete\n")
    
  }
}

save(arrAllInc,chosen_params,params,file="~/srileytmp/events/gemma_20141017/postProc.Rdata")

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

# epiImage <- eventImage(datSim,popgrid,0,0,1000,0,0)
# sum(epiImage$z,na.rm=TRUE)

# png(file=paste(fnOutStem,"figs/map.png",sep=""),width=12,height=12,units="cm",res=300)
#  image(popgrid$x,popgrid$y,log(popgrid$z+0.5),col=rev(grey_pal()(100)))
#  image(epiImage$x,epiImage$y,epiImage$z,add=TRUE)
# dev.off()

# Make a latin hyper cube for the analysis
# nosamples <- 100
# param_scan <- data.frame(
#    R0_Spatial = srg.hyper.vector(nosamples,1.4,2.0,FALSE),
#    Decay_Transmit_Spatial = srg.hyper.vector(nosamples,1,100,TRUE),
#    Offset_Transmit_Spatial = srg.hyper.vector(nosamples,1.5,5.5,FALSE)
#)
#write.table(param_scan,file="~/srileytmp/paramscan.txt",row.names=FALSE,sep=" ")
