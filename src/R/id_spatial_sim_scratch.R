# setwd("~/Dropbox/02_git/id_ss/src/R")
# setwd("/home/sriley/git/id_spatial_sim/src/R")
# source("~/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")

fnPopdata = "../../data/westAfricaAscii.asc"
fnOutStem = "~/Dropbox/projects/ebola/"

# Make an aggregated version of the rater file
popgrid <- read.asciigrid(fnPopdata,as.image=TRUE)
sum(popgrid$z,na.rm=TRUE)
# popgrid_tmp <- raster(popgrid)
# popgrid_coarse_tmp <- aggregate(popgrid_tmp,factor=1000,fun=sum)
# writeRaster(popgrid_coarse_tmp,file="../data/westAfricaAscii_agg1000.asc",format="ascii",
#    overwrite=TRUE)

# epiImage <- eventImage(datSim,popgrid,0,0,1000,0,0)
# sum(epiImage$z,na.rm=TRUE)

# png(file=paste(fnOutStem,"figs/map.png",sep=""),width=12,height=12,units="cm",res=300)
#  image(popgrid$x,popgrid$y,log(popgrid$z+0.5),col=rev(grey_pal()(100)))
#  image(epiImage$x,epiImage$y,epiImage$z,add=TRUE)
# dev.off()

# Make a latin hyper cube for the analysis
nosamples <- 1000
param_scan <- data.frame(
   R0_Spatial = srg.hyper.vector(nosamples,1.3,1.9,FALSE),
   Offset_Transmit_Spatial = srg.hyper.vector(nosamples,1,100,TRUE),
   Decay_Transmit_Spatial = srg.hyper.vector(nosamples,1.5,5.5,FALSE),
   Contact_Trace_Capacity = srg.hyper.vector(nosamples,10,100,TRUE)
)
write.table(param_scan,file="~/srileytmp/paramscan14Nov2014.txt",row.names=FALSE,sep=" ")
