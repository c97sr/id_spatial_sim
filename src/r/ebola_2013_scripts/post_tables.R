# Start from a clean environment if not using immediately after run generation
# File for processing output from gemma's code
rm(list=ls(all=TRUE))
load("~/srileytmp/postProc_Wed_Oct_29_11_54_01_2014.RData",verbose=TRUE)
remakeSummaries <- FALSE
remakeAnneal <- FALSE

# Set you present working directory to find the output
# setwd("~/Dropbox/git/id_spatial_sim/src/R")
# setwd("/home/sriley/git/id_spatial_sim/src/R")

# Load required packages
# require("raster")
# require("sp")
# require("scales")
# require("rgdal")
# require("data.table")
require("adegenet")
require("GenSA")
# require("NMOF")

source("idSimFuncs.R")

reportdir <- "~/Dropbox/projects/ebola/sim_reports/report_0.91/"

# Below here to become a separate file pretty soon
noweeks <- dim(arrAllInc)[1]
nodists <- dim(arrAllInc)[2]
noreals <- dim(arrAllInc)[3]
noparams <- dim(arrAllInc)[4]
testedStats <- c("sumSq","sumSqReport","reportLevel","largeSumSq","thresh5SumSq","thresh5SumSqReport","time20")
compStats <- array(
    dim=c(length(testedStats),dim(arrAllInc)[3],dim(arrAllInc)[4]),
    dimnames=list(testedStats,1:dim(arrAllInc)[3],paramsToLoad)
)

cumweek <- function(mat,window=3) {
  x_tmp <- apply(mat,2,cumsum)
  x_tmp2 <- x_tmp
  x_tmp2[(1+window):(dim(x_tmp)[1]),] <- x_tmp2[(1+window):(dim(x_tmp)[1]),] - x_tmp2[1:(dim(x_tmp)[1]-window),]
  x_tmp2
}

timethresh <- function(vec,thresh=20){
  sum(vec<thresh)
}

# Next try week of x cases
if (remakeSummaries) {
  for (r in 1:noreals) {
    for (p in 1:noparams) {

      # Simple sun of squares
      x <- arrAllInc[,,r,p]
      x_tmp <- apply(x,2,cumsum)
      y_tmp <- apply(y,2,cumsum)
      x_tmp2 <- x_tmp
      x_tmp2[4:(dim(x_tmp)[1]),] <- x_tmp2[4:(dim(x_tmp)[1]),] - x_tmp2[1:(dim(x_tmp)[1]-3),]
      compStats["sumSq",r,p] <- sum((x-y)^2)
      reportingCorr <- max(1,sum(x)/sum(y))
      compStats["sumSqReport",r,p] <- sum((x-y*reportingCorr)^2)
      compStats["reportLevel",r,p] <- reportingCorr
      if (reportingCorr > 1) {
        compStats["largeSumSq",r,p] <- compStats["sumSqReport",r,p]
      } else {
        compStats["largeSumSq",r,p] <- 1e100
      }
      compStats["thresh5SumSq",r,p] <- sum(ifelse(x_tmp2>5,(x-y)^2,y^2))
      compStats["thresh5SumSqReport",r,p] <- sum(ifelse(x_tmp2/reportingCorr > 5,(x-y*reportingCorr)^2,(y*reportingCorr)^2))
      compStats["time20",r,p] <- sqrt(sum((apply(x_tmp,2,timethresh)-apply(y_tmp,2,timethresh))^2))

    }

    cat(paste("completed",r,"of",noreals,"\n"))

  }
  save(compStats,file="./compStats.Rdat")
} else {
  load("./compStats.Rdat",verbose=TRUE)
}

hist((compStats["time20",,]))
plot(paramTab$R0_Spatial,apply(compStats["time20",,],c(2),min))

bestSumSq <- min(compStats["sumSqReport",,])
ibest <- which(compStats[2,,]==bestSumSq,arr.ind=TRUE)
multiBest <- compStats["reportLevel",ibest[1],ibest[2]]
sum(apply(compStats["sumSqReport",,],c(2),min) < 5e5)

pdf(file=paste(reportdir,"fig_R0_sumSqReport.pdf",sep=""))
plot(paramTab$R0_Spatial*1.1,apply(compStats["sumSqReport",,],c(2),min),ylim=c(3e5,5e5))
dev.off()

pdf(file=paste(reportdir,"fig_R0_thresh5SumSqReport.pdf",sep=""))
plot(paramTab$R0_Spatial*1.1,apply(compStats["thresh5SumSqReport",,],c(2),min),ylim=c(0,10e5))
dev.off()

pdf(file=paste(reportdir,"fig_power_sumSqReport.pdf",sep=""))
plot(paramTab$Decay_Transmit_Spatial,apply(compStats["sumSqReport",,],c(2),min),log="y",ylim=c(3e5,5e5))
dev.off()

pdf(file=paste(reportdir,"fig_offset_sumSqReport.pdf",sep=""))
plot(paramTab$Offset_Transmit_Spatial,apply(compStats["sumSqReport",,],c(2),min),log="y",ylim=c(3e5,5e5))
dev.off()

# Use the updated version of the heat chart function to make the log rainbow (lego) version
logbreaks <- c(
    #log10(0+0.5)-1e-10,log10(0.5+0.5)-1e-10,
    log10(0+0.5),
    seq(log10(1+0.5)-1e-10,2.8,length.out=100))
legendscale <- c(0,1,2,3,4,5,10,20,50,100,200,400)
logcols_seasun <- seasun(length(logbreaks)-2)

x <- arrAllInc[,,ibest[1],ibest[2]]/multiBest
sum(x)
sum(y)

inc.heat.chart.pdf.v3(
    log10(y+0.5),
    cols = c(colors()[1],logcols_seasun),
    breaks = logbreaks,
    legendlabs = legendscale,
    legendats = log10(legendscale+0.5),
    vecCountries=distsCountryLO,
    outstem=paste(reportdir,"fig_data_",sep=""),
    xlabs=seq(0,40,5)
)

inc.heat.chart.pdf.v3(
    log10(x+0.5),
    cols = c(colors()[1],logcols_seasun),
    breaks = logbreaks,
    legendlabs = legendscale,
    legendats = log10(legendscale+0.5),
    vecCountries=distsCountryLO,
    outstem=paste(reportdir,"fig_sumSqReport_",sep=""),
    xlabs=seq(0,40,5)
)

bestLargeSumSq <- min(compStats["largeSumSq",,])
ibest2 <- which(compStats["largeSumSq",,]==bestLargeSumSq,arr.ind=TRUE)
x2_all <- arrAllInc[,,ibest2[1],]

paramTab[147,]
sum(apply(compStats["sumSqReport",,],c(2),min) < 5e5)
urFactor <- compStats["reportLevel",ibest2[1],ibest2[2]]
x2 <- arrAllInc[,,ibest2[1],ibest2[2]]/multiBest
x2_cum <- cumweek(x2)
inc.heat.chart.pdf.v3(
    log10(ifelse(x2_cum>5,x2,0)+0.5),
    cols = c(colors()[1],logcols_seasun),
    breaks = logbreaks,
    legendlabs = legendscale,
    legendats = log10(legendscale+0.5),
    vecCountries=distsCountryLO,
    outstem=paste(reportdir,"fig_largeSumSq_",sep=""),
    xlabs=seq(0,40,5)
)

tmp <- min(compStats["time20",,])
ibest3 <- which(compStats["time20",,]==tmp,arr.ind=TRUE)
paramTab[ibest3[1],]
x3 <- arrAllInc[,,ibest3[1],ibest3[2]]
x3_cum <- cumweek(x3)
sum(x3)
inc.heat.chart.pdf.v3(
  log10(ifelse(x3_cum > 5,x3,0) + 0.5),
  cols = c(colors()[1],logcols_seasun),
  breaks = logbreaks,
  legendlabs = legendscale,
  legendats = log10(legendscale+0.5),
  vecCountries=distsCountryLO,
  outstem=paste(reportdir,"fig_time20_",sep=""),
  xlabs=seq(0,40,5)
)

# Plot some goodness of fits charts
fit_x3 <- ifelse(x3>0 & y>1,x3/y,NA)
fit_x2 <- ifelse(x2>0 & y>1,x2/y,NA)

pdf(file=paste(reportdir,"fig_gof_ss_time20_control.pdf",sep=""))
  hist(ifelse(x2 > 0 & y > 0,x2-y,NA),breaks=seq(-300,300,5)-2.5,xlim=c(-100,100),
      main="Sum of square fit",xlab="Sim - data (where either non zero)")
  hist(ifelse(x3 > 0 & y > 0,x3-y,NA),breaks=seq(-300,300,5)-2.5,xlim=c(-100,100),
      main="Time to 20 cases fit",xlab="Sim - data (where either non zero)")
dev.off()

hist(ifelse(x2_all > 0 & replicate(400,x2) > 0,replicate(400,x2)-x2_all,NA),breaks=seq(-10000,10000,5)-2.5,xlim=c(-100,100),
    main="Best time fit - all sims (same parameters)",
    xlab="Sim random - sim data (where either non zero)")


image(log10(fit_x2),,col=seasun(100))
image(fit_x3,breaks=0:100,col=seasun(100))

topLogSumSq <- apply(log10(compStats["sumSqReport",,]),c(2),function(vec){mean(sort(vec)[1:5])})
topThreshLogSumSq <- apply(log10(compStats["thresh5SumSq",,]),c(2),function(vec){mean(sort(vec)[1:5])})
pdf(file=paste(reportdir,"fig_R0_thresh5SumSq.pdf",sep=""))
  plot(paramTab$R0_Spatial,topLogSumSq,ylim=c(5.5,6))
dev.off()

pdf(file=paste(reportdir,"fig_R0_topThreshLogSumSq.pdf",sep=""))
  plot(paramTab$R0_Spatial,topThreshLogSumSq,ylim=c(5.5,6))
dev.off()

maskLow <- topLogSumSq < 5.7
sum(topLogSumSq < 5.7)
pdf(file=paste(reportdir,"fig_correlation_good_power_offsets.pdf",sep=""))
plot(paramTab$Decay_Transmit_Spatial[maskLow],paramTab$Offset_Transmit_Spatial[maskLow])
dev.off()
plot(paramTab$Offset_Transmit_Spatial,topLogSumSq,ylim=c(5.5,6.0))
pdf(file=paste(reportdir,"fig_ratio_offset_over_power.pdf",sep=""))
plot(paramTab$Offset_Transmit_Spatial/paramTab$Decay_Transmit_Spatial,topLogSumSq,ylim=c(5.5,6.0))
dev.off()

# Define a nearest neighbor matrix for simmulated annealing
nnearNeigh <- 5
names(paramTab)
pMins <- c(1.2,1,1.5)
pMaxes <- c(1.6,100,5.5)
pLogs <- c(FALSE,TRUE,FALSE)
unitParams <- paramTab
unitParams[] <- NA
unitParams$R0_Spatial <- (paramTab$R0_Spatial-1.2)/(1.6-1.2)
unitParams$Decay_Transmit_Spatial <- (paramTab$Decay_Transmit_Spatial-1.5)/(5.5-1.5)
unitParams$Offset_Transmit_Spatial <- log10(paramTab$Offset_Transmit_Spatial)
unitParams$Offset_Transmit_Spatial <- unitParams$Offset_Transmit_Spatial/2
distFunc <- function(X,Y,unitVals){
  rtn <- 	(unitVals[X,1]-unitVals[Y,1])^2 +
      (unitVals[X,2]-unitVals[Y,2])^2 +
      (unitVals[X,3]-unitVals[Y,3])^2
  rtn <- sqrt(rtn)
  rtn
}
distLookup <- outer(1:400,1:400,distFunc,unitVals=unitParams)
indexNeigh <- t(apply(distLookup,c(1),function(x){order(x)[2:(2+nnearNeigh)]}))

# Simulate simulated annealing!
# Needs rerunning for about the number of runs we could do overnight
# Needs a temperature cooling mechanism to be tested and also to be tested from
# different starting points
simsimAnneal <- function(n) {
  nsamples 	<- n
  cur_index <- sample.int(400,1)
  cur_stat 	<- topLogSumSq[cur_index]
  rtn 			<- matrix(nrow=nsamples,ncol=1+length(names(paramTab)))
  csamp 		<- 1
  finalTmp	<- log10(0.001)
  initTmp 	<- log10(0.1)
  while (csamp <= nsamples) {
    logTemp 	<- 10^(initTmp + (csamp/nsamples)*(finalTmp-initTmp))
    rtn[csamp,] <- t(as.numeric(c(paramTab[cur_index,],cur_stat)))
    whichNeigh <- sample.int(5,1)
    next_index <- indexNeigh[cur_index,whichNeigh]
    next_stat <- topLogSumSq[next_index]
    accept <- FALSE
    if (next_stat < cur_stat) {
      accept <- TRUE
    } else {
      acceptProb <- exp(-(next_stat-cur_stat)/logTemp)
      if (runif(1) < acceptProb) {
        accept <- TRUE
      }
    }
    if (accept) {
      cur_index <- next_index
      cur_stat <- next_stat
    }
    csamp <- csamp + 1
  }
  rtn[nsamples,]
}
simsimAnneal(10)

repSimSimAnneal <- function(n,m) {
  rtn <- matrix(nrow=n,ncol=4)
  for (i in 1:n) {
    rtn[i,] <- simsimAnneal(m)
  }
  rtn
}

if (remakeAnneal) {
  system.time(annealRes <- repSimSimAnneal(100,2000))
  save(annealRes,file="annealRes.Rdat")
} else {
  load("annealRes.Rdat",verbose=TRUE)
}

annealRes <- annealRes[order(annealRes[,4]),]

pdf(file=paste(reportdir,"fig_simanneal_R0.pdf",sep=""))
plot(
    jitter(annealRes[,1],amount=0.005),
    jitter(annealRes[,4],amount=0.03),
    pch=19,
    cex=0.5,
    xlim=c(1.2,1.6),
    ylim=c(5.5,6.2))
dev.off()

pdf(file=paste(reportdir,"figure_simanneal_power_offset_ratio.pdf",sep=""))
plot(
    jitter(annealRes[,2]/annealRes[,3],amount=0.1),
    jitter(annealRes[,4],amount=0.03),
    pch=19,
    cex=0.5,
    xlim=c(1,30),
    ylim=c(5.5,6.2))
dev.off()

