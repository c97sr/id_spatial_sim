
rm(list=ls())
ct_data1 <- read.csv(file="C:/Users/kkokwok/Dropbox/id_spatial_sim/scenarios/hk_contacts/output/ebola_1_pset_0_Events.out", header=TRUE, sep="")

a<- table(ct_data1$Day,ct_data1$Event)/10

b <- sum(a[,1])
plot(a[,1],main="Simulated Epidemic Curve of SARS (contact tracing at day 50)",sub=bquote("Network Ro=3; n=7.5; million total number of infection"==.(b))
,xlab="Day", ylab="number of infection")

ct_data2 <- read.csv(file="C:/Users/kkokwok/Dropbox/id_spatial_sim/scenarios/hk_contacts/output/ebola_2_pset_0_Events.out", header=TRUE, sep="")


c<- table(ct_data2$Day,ct_data2$Event)/10

d <- sum(c[,1])
plot(c[,1],main="Simulated Epidemic Curve of SARS (contact tracing at day 50)",sub=bquote("Network Ro=3; n=7.5; million total number of infection"==.(b))
,xlab="Day", ylab="number of infection")


ct_data3 <- read.csv(file="C:/Users/kkokwok/Dropbox/id_spatial_sim/scenarios/hk_contacts/output/ebola_3_pset_0_Events.out", header=TRUE, sep="")
ct_data4 <- read.csv(file="C:/Users/kkokwok/Dropbox/id_spatial_sim/scenarios/hk_contacts/output/ebola_4_pset_0_Events.out", header=TRUE, sep="")

e<- table(ct_data3$Day,ct_data3$Event)/10

f <- sum(e[,3])


g<- table(ct_data4$Day,ct_data4$Event)/10

h <- sum(g[,3])

ct_data5 <- read.csv(file="C:/Users/kkokwok/Dropbox/id_spatial_sim/scenarios/hk_contacts/output/ebola_5_pset_0_Events.out", header=TRUE, sep="")
ct_data6 <- read.csv(file="C:/Users/kkokwok/Dropbox/id_spatial_sim/scenarios/hk_contacts/output/ebola_6_pset_0_Events.out", header=TRUE, sep="")

i<- table(ct_data5$Day,ct_data5$Event)/10

j <- sum(i[,3])


k<- table(ct_data6$Day,ct_data6$Event)/10

l <- sum(k[,3])


plot(c(as.numeric(rownames(a))),a[,1],xlim=c(0,100),ylim=c(0,1000),xlab="Day",ylab="Number of infection",main="Black: Vaccine Probability susceptible=1;Red: Vaccine Probability susceptible=0")
plot(c(as.numeric(rownames(a))),a[,1],xlim=c(0,100),xlab="Day",ylab="Number of infection",main="Black: Vaccine Probability susceptible=1;Red: Vaccine Probability susceptible=0")

points(c(as.numeric(rownames(c))), c[,1], col="red", pch="*")
points(c(as.numeric(rownames(e))), e[,1], col="blue", pch="#")
points(c(as.numeric(rownames(g))), g[,1], col="orange", pch="?")

#points(c(as.numeric(rownames(i))), i[,1], col="yellow", pch="?")
#points(c(as.numeric(rownames(k))), k[,1], col="green", pch="?")
