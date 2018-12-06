rm(list=ls(all=TRUE))
options(error=NULL)
#source("SR_misc.r")

fnEpidata = "~/Dropbox/svneclipse/spatialsim/spatialsim/runtemplates/rfcid_rep/output/flusars_2_pset_0_Events.out"
fnOutfile = "~/Dropbox/tmp/fig3.pdf"

# This needs to be tuned to be a function that gives incidence based on the event files from large simulations

data <- read.table(file=fnEpidata,header=TRUE)

decLongLatDist <- function(x1,y1,x2,y2,translate=FALSE) {
	if (translate) {
		x1 <- x1*pi/180
		x2 <- x2*pi/180
		y1 <- y1*pi/180
		y2 <- y2*pi/180
	}
	rtn <- 6378.7*acos(sin(y1)*sin(y2)+cos(y1)*cos(y2)*cos(x2-x1))
	rtn
}

discinc <- function(data,event=0,dist=10,lat=0,long=0,sttime=0,dt=1,notime=10,ordreals=0:0) {
	rtn <- array(0,c(notime,length(ordreals)))
	rc 	<- 1
	norows <- dim(data)[1]
	maxtime <- sttime + dt * notime
	target_real_ind <- 1
	maxrealind <- length(ordreals)
	while (rc <= norows) {
		if (data$Run[rc] > ordreals[target_real_ind] && target_real_ind <= maxrealind) {
			target_real_ind <- target_real_ind + 1
		}
		if (data$Event[rc]==event &&
			data$Run[rc] == ordreals[target_real_ind] &&
			data$Day[rc] >= sttime && data$Day[rc] < maxtime &&
			decLongLatDist(data$X[rc],data$Y[rc],long,lat) < dist ) {
			time_index <- floor(data$Day[rc] - sttime) + 1
			rtn[time_index,ordreals[target_real_ind]+1] <- rtn[time_index,ordreals[target_real_ind]+1]+ 1
		}
		rc <- rc + 1
	}
	rtn
}

gzlat <- 23.0797309875 
gzlong <- 113.203125
hklat <- 22.3964
hklong <- 114.1095

xvals <- discinc(data,long=gzlong*pi/180,lat=gzlat*pi/180,dist=20,notime=400,dt=0.25)
yvals <- discinc(data,long=hklong*pi/180,lat=hklat*pi/180,dist=20,notime=400,dt=0.25)
for (i in 2:length(xvals)) xvals[i] <- xvals[i-1] + xvals[i]
for (i in 2:length(yvals)) yvals[i] <- yvals[i-1] + yvals[i] 
plot(xvals,yvals,type="l",xlim=c(0,max(v=xvals,yvals)),ylim=c(0,max(xvals,yvals)))

# Example for debugging
# discinc(data,long=114.278*pi/180,lat=22.413*pi/180,dist=10000)

winwidth <- 10/cm(1)
winheight <- 20/cm(1)
par(cex=0.8)
par(mai=c(0/cm(1),0/cm(1),0/cm(1),0/cm(1)))

# windows(width=winwidth,height=winheight)
# jpeg(file=paste(fnOutfile,".jpg",sep=""),width=1024,height=768)
pdf(file=paste(fnOutfile,".pdf",sep=""),width=winwidth,height=winheight)

par(fig= sr_chart_pos(1,1,1,2))
plot(1:2)
par(fig= sr_chart_pos(1,2,1,2),new=TRUE)
plot(1:2)

dev.off()