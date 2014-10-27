rm(list=ls(all=TRUE))
# Set you present working directory to find the output
# setwd("~/Dropbox/git/id_spatial_sim/src/R")
# setwd("/home/sriley/git/id_spatial_sim/src/R")

# Load required packages
require("raster")
require("sp")
require("scales")
require("rgdal")
require("data.table")
require("adegenet")
require("NMOF")
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

# Merge away the Conakry issues temporarily by making all conakry be ratoma
setNonRatomaConakry <- c("KALOUM","DIXINN","MATAM","MATOTO")
dat$district <- ifelse(dat$district %in% setNonRatomaConakry,"RATOMA",dat$district)

# Trim the data down to the first 1000
# dat <- dat[order(dat$EpiDay)[1:1000],]

# Load the shape files
shapeDir <- "/Users/sriley/srileytmp/sfs"
dists <- readOGR(dsn=shapeDir,layer="ThreeCountries")

# Reorder the districts by country, then latitude
distsLatOrder <- as.character(dists$ADM2_NAME[rev(order(dists$CENTER_LAT))])
distsCountryLO <- as.character(dists$ADM0_NAME[rev(order(dists$CENTER_LAT))])
write.csv(data.frame(WHO_countries=distsCountryLO,WHO_districts=distsLatOrder),file="~/Dropbox/tmp/WHO_Dists.csv",row.names=FALSE)

y <- (
      make.incidence.from.linelist(
          distsLatOrder,dat$district,dat$EpiWeek,DTs=weeksUsed
      )
      )$inctab

# Temp fix to remove non ratoma conakry
# checked with sum that doesn't remove cases
y <- y[,-match(setNonRatomaConakry,colnames(y))]
distsCountryLO <- distsCountryLO[-match(setNonRatomaConakry,distsLatOrder)] 

# Start to look at the table files
tabledir <- "~/srileytmp/event_files/20141021/"
paramfname <- "paramscan.txt"

# Load the parameter table
paramTab <- read.table(file=paste(tabledir,paramfname,sep=""),header=TRUE,sep=" ")
names(paramTab) <- c("R0_Spatial","Offset_Transmit_Spatial","Decay_Transmit_Spatial")

# Load the place lookups
tableDistricts <- read.csv("../../data/DistrictLookup.csv",stringsAsFactors=FALSE)

paramsToLoad <- 1:200
realsToLoad <- 0:199

arrAllInc <- array(
    dim=c(dim(y)[1],dim(y)[2],length(realsToLoad),length(paramsToLoad)),
    dimnames=list(1:(dim(y)[1]),colnames(y),realsToLoad,paramsToLoad))

# Read in the files using known numbers of params and realizations
system.time(
for (p in 1:length(paramsToLoad)) {
  for (r in 1:length(realsToLoad)) {
    pnum <- paramsToLoad[p]
    rnum <- realsToLoad[r]
    sTT <- fread(paste(tabledir,
        "WestAfrica_paramset",pnum,".",rnum,".adunit.csv",sep=""),
        select=65:127,header=FALSE)
    setnames(sTT,names(singleTestTab),tableDistricts[,"WHO_districts"])
    tmp <- convert.spatialsim.WHO.V1(sTT,colnames(y),tableDistricts,length(weeksUsed))
    arrAllInc[,,r,p] <- tmp
  }
  cat("Completed",p,"of",length(paramsToLoad),"\n")
}
)

save(arrAllInc,y,paramsToLoad,paramTab,file=paste("~/srileytmp/postProc_",gsub(" ","_",gsub(":","_",date())),".Rdata",sep=""))

reportdir <- "~/Dropbox/projects/ebola/sim_reports/"

load("~/srileytmp/postProc_Sun_Oct_26_14_14_34_2014.RData",verbose=TRUE)
# Quick debug
singleTestTab <- fread(paste(tabledir,
        "WestAfrica_paramset",92,".",22,".adunit.csv",sep=""),header=FALSE)


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

bestSumSq <- min(compStats["sumSqReport",,])
ibest <- which(compStats[2,,]==bestSumSq,arr.ind=TRUE)
multiBest <- compStats["reportLevel",ibest[1],ibest[2]]
paramTab[22,]
sum(apply(compStats["sumSqReport",,],c(2),min) < 5e5)

pdf(file=paste(reportdir,"figure1_R0.pdf",sep=""))
plot(paramTab$R0_Spatial,apply(compStats["sumSqReport",,],c(2),min),ylim=c(3e5,5e5))
dev.off()
plot(paramTab$R0_Spatial,apply(compStats["sumSqReport",,],c(2),min))
pdf(file=paste(reportdir,"figure2_power.pdf",sep=""))
plot(paramTab$Decay_Transmit_Spatial,apply(compStats["sumSqReport",,],c(2),min),log="y",ylim=c(3e5,5e5))
dev.off()
plot(paramTab$Decay_Transmit_Spatial,apply(compStats["sumSqReport",,],c(2),min),log="y")
pdf(file=paste(reportdir,"figure3_offset.pdf",sep=""))
plot(paramTab$Offset_Transmit_Spatial,apply(compStats["sumSqReport",,],c(2),min),log="y",ylim=c(3e5,5e5))
dev.off()
plot(paramTab$Offset_Transmit_Spatial,apply(compStats["sumSqReport",,],c(2),min),log="y")
hist(log10(apply(arrAllInc,c(3,4),sum)))

x <- arrAllInc[,,ibest[1],ibest[2]]/multiBest

sum(x)
sum(y)

# Use the updated version of the heat chart function to make the log rainbow (lego) version
logbreaks <- c(
    #log10(0+0.5)-1e-10,log10(0.5+0.5)-1e-10,
    log10(0+0.5),
    seq(log10(1+0.5)-1e-10,2.8,length.out=100))
legendscale <- c(0,1,2,3,4,5,10,20,50,100,200,400)
logcols_seasun <- seasun(length(logbreaks)-2)

inc.heat.chart.pdf.v3(
    log10(y+0.5),    
    cols = c(colors()[1],logcols_seasun),
    breaks = logbreaks,
    legendlabs = legendscale,
    legendats = log10(legendscale+0.5),
    vecCountries=distsCountryLO,
    outstem=paste(reportdir,"figure4_data_",sep=""),
    xlabs=seq(0,40,5)
)

inc.heat.chart.pdf.v3(
    log10(x+0.5),    
    cols = c(colors()[1],logcols_seasun),
    breaks = logbreaks,
    legendlabs = legendscale,
    legendats = log10(legendscale+0.5),
    vecCountries=distsCountryLO,
    outstem=paste(reportdir,"figure5_best_ASOS_",sep=""),
    xlabs=seq(0,40,5)
)

bestLargeSumSq <- min(compStats["largeSumSq",,])
ibest2 <- which(compStats["largeSumSq",,]==bestLargeSumSq,arr.ind=TRUE)

paramTab[147,]
sum(apply(compStats["sumSqReport",,],c(2),min) < 5e5)
urFactor <- compStats["reportLevel",ibest2[1],ibest2[2]]
x2 <- arrAllInc[,,ibest2[1],ibest2[2]]/multiBest

inc.heat.chart.pdf.v3(
    log10(x2+0.5),    
    cols = c(colors()[1],logcols_seasun),
    breaks = logbreaks,
    legendlabs = legendscale,
    legendats = log10(legendscale+0.5),
    vecCountries=distsCountryLO,
    outstem=paste(reportdir,"figure6_best_ASOS_no_under_report_",sep=""),
    xlabs=seq(0,40,5)
)

# optSS <- function(vecTheta) {
#   
#   offset <- vecTheta[1]
#   under.report <- vecTheta[2]
#   
# }

# xtmp <- x
# for (offset in 1:maxoffset) {
#   xtmp[(1+offset):noweeks,] <- x[1:(noweeks-offset),]
#   xtmp[1:offset,] <- 0
#   offmeasures <- 
# }
