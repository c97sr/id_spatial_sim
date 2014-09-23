rm(list=ls(all=TRUE))
# setwd("~/Dropbox/svneclipse/spatialsim/spatialsim/src/rcode")

# options(error=recover)
# options(error=NULL)
# require("fields")
# require("lattice")
require("sp")
require("scales")
source("../../../../idsource/R/stevensRfunctions.R")

fnPopdata = "../../data/prd_sim_ascii_n_44384655.txt"
fnEpidata = "../../runtemplates/ebola/output/flusars_2_pset_0_Events.out"
fnOutStem = "~/Dropbox/tmp/debg"

popgrid <- read.asciigrid(fnPopdata,as.image=TRUE)
popgrid$x <- popgrid$x * pi / 180
popgrid$y <- popgrid$y * pi / 180
data <- read.table(file=fnEpidata,header=TRUE)
epiImage <- eventImage(data,popgrid,0,0,1000,0,0)

pdf(file=paste(fnOutStem,"pop.pdf",sep=""))
	image(popgrid$x,popgrid$y,log(popgrid$z+0.5),col=grey_pal()(100))
	image(epiImage$x,epiImage$y,epiImage$z,add=TRUE)
dev.off()
