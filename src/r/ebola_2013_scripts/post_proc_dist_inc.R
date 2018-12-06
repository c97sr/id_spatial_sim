rm(list=ls(all=TRUE))

# Set the present working directory to find the output
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

load("~/Dropbox/tmp/postProc19oct.Rdata",verbose=TRUE)
# save(arrAllInc,chosen_params,params,file="~/srileytmp/events/gemma_20141017/postProc.Rdata")

dim(arrAllInc)
npsets <- dim(arrAllInc)[4]
nreals <- dim(arrAllInc)[3]
diff_stats <- c("sum_sq","thesh_sum_sq")
nstats <- length(diff_stats)

sumArr <- array(dim=c(nstats,nreals,npsets))

for (p in 1:npsets) {
  for (r in 1:nreals) {
    sumArr[1,r,p] <- sum((arrAllInc[,,r,p]-y)^2)
    sumArr[2,r,p] <- sum(ifelse(arrAllInc[,,r,p]>5,(arrAllInc[,,r,p]-y)^2,y^2))
  }
}

# Up to here. Need to trawl through for best fit
log10Sum <- log10(sumArr)

plot(params$power,apply(sumArr[2,,101:200],c(2),min),ylim=c(3e5,6e5))
minval <- sort(sumArr[1,,101:200])[2]
ibest <- which(sumArr[1,,101:200]==minval,arr.ind=TRUE) 
apply(sumArr[1,,101:200],c(2),min)
params[21,]

x <- y
x[] <- arrAllInc[,,ibest[1],ibest[2]]
sum(x)
sum(y)

rownames(x) <- rownames(y)
inc.heat.chart.pdf.v2(
    x,
    vecCountries=distsCountryLO,
    outstem=paste("~/srileytmp/best_sim",sep=""),
    xlabs=seq(0,40,5)
)

inc.heat.chart.pdf.v2(
    y,
    vecCountries=distsCountryLO,
    outstem=paste("~/srileytmp/data_",sep=""),
    xlabs=seq(0,40,5)
)

#logbreaks <- seq(-0.30102-1e-10,2.7,0.176+0.30103)
logbreaks <- c(
    #log10(0+0.5)-1e-10,log10(0.5+0.5)-1e-10,
    log10(0+0.5),
    seq(log10(1+0.5)-1e-10,2.8,length.out=100))
legendscale <- c(0,1,2,3,4,5,10,20,50,100,200,400)
logcols <- rainbow(length(logbreaks)-2,start=0,end=0.9) 
inc.heat.chart.pdf.v2(
    log10(x+0.5),    
    cols = c(colors()[1],logcols),
    breaks = logbreaks,
    legendlabs = legendscale,
    legendats = log10(legendscale+0.5),
    vecCountries=distsCountryLO,
    outstem=paste("~/srileytmp/best_sim_log",sep=""),
    xlabs=seq(0,40,5)
)
inc.heat.chart.pdf.v2(
    log10(y+0.5),    
    cols = c(colors()[1],logcols),
    breaks = logbreaks,
    legendlabs = legendscale,
    legendats = log10(legendscale+0.5),
    vecCountries=distsCountryLO,
    outstem=paste("~/srileytmp/data_log",sep=""),
    xlabs=seq(0,40,5)
)

write.table(t(x),file="~/Dropbox/tmp/sim_data_for_wes.txt",
    sep="\t",row.names=TRUE,col.names=FALSE)

