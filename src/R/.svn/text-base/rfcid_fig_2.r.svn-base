# This fiel now superceeded by sps.plot.snaposhots in stevensRfunctions

# Assumes current working directory to be spatialsim or in a folder in the same directory
# as spatial sim. If required change working directory below.
# setwd(XXX)

rm(list=ls(all=TRUE))
options(error=NULL)
require("sp")
require("fields")
# source("/Users/sriley/Dropbox/svneclipse/spatialsim/spatialsim/src/rcode/utilities.r")
source("/Users/sriley/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")

fnPopdata = "~/Dropbox/svneclipse/spatialsim/spatialsim/data/prd_sim_ascii_n_44384655.txt"
fnEpidata = "~/Dropbox/svneclipse/spatialsim/spatialsim/runtemplates/rfcid_rep/output/flusars_2_pset_1_Events.out"
fnOutfile = "~/Dropbox/tmp/sev_flu_cont"

popgrid <- read.asciigrid(fnPopdata,as.image=TRUE)
data <- read.table(file=fnEpidata,header=TRUE)
data$X <- data$X * 180 / pi
data$Y <- data$Y * 180 / pi

winwidth <- 20/cm(1)
winheight <- 10/cm(1)
resolution=300
popaxis <- grey(0.9-0.6*0:100/100)
infaxis <- tim.colors(256)
ximax=4
yimax=2
chtlet <- c("a)","b)","c)","d)","e)","f)","g)","h)","i)","j)","k)","l)","m)","n)","o)","p","q")
chttimes <- c(2,4,6,8,10,12,14,16,18,20,22)
# chttimes <- c(6,12,18,24,30,36,42,48,54,60,55)

# epiImage <- eventImage(data,popgrid,0,0,chttimes[1],0,0)

# X11(width=winwidth,height=winheight)
# jpeg(file=paste(fnOutfile,".jpg",sep=""),width=1024,height=768)
pdf(file=paste(fnOutfile,".pdf",sep=""),width=winwidth,height=winheight)

par(cex=0.8)
par(mai=c(0/cm(1),0/cm(1),0/cm(1),0/cm(1)))

for (yi in yimax:1) {
	for (xi in 1:ximax) {	
		if (xi==1 && yi==yimax) blnew=FALSE 
		else blnew=TRUE
		cur_pos <- srg.chart.pos(xi,yi,ximax,yimax,xg=0.05/winwidth,yg=0.05/winheight) 
		par(fig=cur_pos,new=blnew)
		cur_index <- xi+(yimax-yi)*ximax
		epiImage <- sps.event.image(data,popgrid,0,0,chttimes[cur_index],0,0)
		image(popgrid,col=popaxis,axes=FALSE)
		image(epiImage,col=infaxis,add=TRUE)
		mtext(paste(chtlet[cur_index],chttimes[cur_index],"days"),side=3,line=-1)
	}	
}

dev.off()