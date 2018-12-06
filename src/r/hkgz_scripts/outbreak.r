rm(list=ls(all=TRUE))
# setwd("/Users/sriley/Dropbox/git/id_spatial_sim/src/R")

# options(error=recover)
# options(error=NULL)
require("fields")
require("lattice")
require("sp")
require("scales")
source("../../../../svneclipse/idsource/R/stevensRfunctions.R")

fnPopdata = "../../data/prd_sim_ascii_n_44384655.txt"
fnEpidata = "../../runtemplates/ebola/output/flusars_2_pset_0_Events.out"
fnOutfile = "~/Dropbox/tmp/ebola_dbg"
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

pdf(file=paste(fnOutStem,"test.pdf",sep=""))
image(epiImage$x,epiImage$y,epiImage$z)
dev.off()
# image.plot(epiImage)

map_xmin = min(popgrid$x)
map_xmax = max(popgrid$x)
map_ymin = min(popgrid$y)
map_ymax = max(popgrid$y)

x_limits = c(map_xmin,map_xmax)
y_limits = c(map_ymin,map_ymax)

w_width = 8
w_height = 768/1024*w_width
r_margin = 0.3*w_width
m_width = w_width-r_margin

display_x_mid = m_width/2/w_width
display_y_mid = 0.5

display_aspect = m_width/w_height
map_aspect = (map_xmax-map_xmin)/(map_ymax-map_ymin)

if (map_aspect > display_aspect) {
	display_height = 1-0.03
	display_width = display_height * map_aspect
} else {
	display_width = m_width/w_width-0.03
	display_height = display_width / map_aspect
}

map_location=c(display_x_mid-display_width/2,display_x_mid+display_width/2,display_y_mid-display_height/2,display_y_mid+display_height/2)
epiImage <- eventImage(data,popgrid,0,0,1000,0,0)

windows(width=w_width,height=w_height)
par(fig=c(0,1,0,1),mai=(c(0,0,0,r_margin)))
legend1_location=c(0.7,0.75,0.1,0.9)
legend2_location=c(0.81,0.86,0.1,0.9)
image.plot(epiImage,legend.only=TRUE,col=tim.colors(),smallplot=legend1_location)
image.plot(popgrid,col=grey(0.9-0.6*0:100/100),add=TRUE,legend.only=TRUE,smallplot=legend2_location)

# par(fig=map_location,mai=c(0,0,0,0),new=TRUE)
par(fig=c(0,0.5,0,0.5),mai=c(0,0,0,0),new=TRUE)

plot(1:2,xlim=x_limits,ylim=y_limits,type="n",axes=FALSE,new=FALSE)
image(popgrid,col=grey(0.9-0.6*0:100/100),add=TRUE)
image(epiImage,col=tim.colors(256),add=TRUE)

savePlot(filename=fnOutfile,type="bmp")