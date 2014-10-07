rm(list=ls(all=TRUE))

# Set you present working directory to find the output
# setwd("~/Dropbox/projects/ebola/sims/runs_20140710")
require("raster")

# options(error=recover)
# options(error=NULL)
# require("fields")
# require("lattice")
require("sp")
require("scales")
source("~/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")

fnPopdata = "../popdata/westAfricaAscii.asc"
fnEpidata = "./output/batch_1__pset_0_Events.out"
fnOutStem = "~/Dropbox/tmp/debg"

popgrid <- read.asciigrid(fnPopdata,as.image=TRUE)
coarse_popgrid <- aggregate(popgrid,factor=2)
sum(popgrid$z,na.rm=TRUE)
popgrid$x <- popgrid$x * pi / 180
popgrid$y <- popgrid$y * pi / 180
dat <- read.table(file=fnEpidata,header=TRUE)
epiImage <- eventImage(dat,popgrid,0,0,1000,0,0)

runs <- as.numeric(names(table(dat$Run)))

for (g in c(1)) {
  png(file=paste(fnOutStem,"run_",g,".png",sep=""),width=7,height=7,units="in",res=300)
  image(popgrid$x,popgrid$y,log(popgrid$z+0.5),col=rev(grey_pal()(100)))
  image(epiImage$x,epiImage$y,epiImage$z,add=TRUE)
  dev.off()
}

infs <- dat[dat$Event == 0,]
mean(table(infs$Run,infs$Generation)[,2]/table(infs$Run,infs$Generation)[,1])







