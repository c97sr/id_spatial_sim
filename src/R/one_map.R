rm(list=ls(all=TRUE))
# setwd("~/Dropbox/git/id_spatial_sim/src/R")
require("raster")

# options(error=recover)
# options(error=NULL)
# require("fields")
# require("lattice")
require("sp")
require("scales")
source("~/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")

fnPopdata = "../../data/monroviaAscii.asc"
fnEpidata = "../../scenarios/ebola/output/ebola_pset_0_Events.out"
fnOutStem = "~/Dropbox/tmp/debg"

popgrid <- read.asciigrid(fnPopdata,as.image=TRUE)
popgrid$x <- popgrid$x * pi / 180
popgrid$y <- popgrid$y * pi / 180
dat <- read.table(file=fnEpidata,header=TRUE)
epiImage <- eventImage(dat,popgrid,0,0,1000,0,0)

pdf(file=paste(fnOutStem,"pop.pdf",sep=""))
	image(popgrid$x,popgrid$y,log(popgrid$z+0.5),col=rev(grey_pal()(100)))
	image(epiImage$x,epiImage$y,epiImage$z,add=TRUE)
dev.off()

tmp <- dat[dat$Event == 0,]
names(tmp)
x <- table(tmp$Run,tmp$Event)
x






