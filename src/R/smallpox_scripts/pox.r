workdir <- "D:\\share\\files\\projects\\pox2002\\results\\060425_trick"

dataA <- read.table(file=paste(workdir,"\\wps_24.00_wpp_0.5000_pow_2.29_erm_2.90786_era_16.2468_r0s_0.132726_ar_0.573833_hrf_8.88887_hhn_7.74194_pset_0_Events.out",sep=""),header=TRUE,sep="\t")
dataB <- read.table(file=paste(workdir,"\\wps_24.00_wpp_0.5000_pow_2.29_erm_3.34983_era_5.07597_r0s_0.702689_ar_0.647774_hrf_6.38485_hhn_5.09581_pset_0_Events.out",sep=""),header=TRUE,sep="\t")
dataC <- read.table(file=paste(workdir,"\\wps_48.00_wpp_0.2500_pow_2.29_erm_3.66869_era_16.7357_r0s_0.568637_ar_0.851063_hrf_2.4653_hhn_1.93623_pset_0_Events.out",sep=""),header=TRUE,sep="\t")

# test line here

test <- function(data) {
	currentrow=1
	currentinf=currentrow
	currenttest=currentinf
	indexinf=0
	indexvacc=1
	maxrow=dim(data)[1]
	maxsteps=7
	currentreal=0
	## longlist = array(0,c(0))
	returnval = array(0,c(maxsteps+2))
	while (currentrow != maxrow) {
		if (data$Event[currentrow]==indexinf && data$Day[currentrow]!= 0) {
			inftime=data$Day[currentrow]
			currentnode=data$Index[currentrow]
			currentrun=data$Run[currentrow]
			currenttime=data$Time[currentrow]
			endrow=currentrow
			while (endrow != maxrow && data$Run[endrow] == currentrun) endrow=endrow+1
			val = -1
			for (i in currentrow:endrow) {
				if (	data$Index[i] == currentnode
					&& data$Event[i] == indexvacc) val = data$Day[i]-inftime
			}
			## longlist = c(longlist,val)
			if (val < 0) returnval[1] = returnval[1]+1
			else if (val <= maxsteps) returnval[val+1] = returnval[val+1]+1
			else returnval[maxsteps+2] = returnval[maxsteps+2]+1
		}
		currentrow=currentrow+1
	}
	returnval
}

a = test(dataA)
b = test(dataB)
c = test(dataC)

write.table(a,file=paste(workdir,"\\ScenA_tab.out",sep=""),sep="\t")
write.table(b,file=paste(workdir,"\\ScenB_tab.out",sep=""),sep="\t")
write.table(c,file=paste(workdir,"\\ScenC_tab.out",sep=""),sep="\t")

table(dataA$Event)
table(dataB$Event)
table(dataC$Event)
