rm(list=ls())        
ct_data1a <- read.csv(file="C:/Users/kkokwok/Dropbox/id_spatial_sim/scenarios/hk_contacts/output/ncov_1a_pset_0_Events.out", header=TRUE, sep="")
n<-200
a<-0
b<-0
temp<-rep(0,times=n)
count<-rep(0,times=n)
l <-n-1
for(k in 0:l) {
test_data <- ct_data1a[((ct_data1a$Event==0)&(ct_data1a$Run==k)),]


for(i in 1:dim(test_data)[1]) {
    for(j in 1:dim(test_data)[1]) {
     
if (test_data$Index[i]==test_data$infector[j]) {

        temp[k+1]= temp[k+1] + test_data$Day[j]-test_data$Day[i]
        count[k+1]= count[k+1] + 1

      }

}
print(k+1)
print(temp[k+1]/count[k+1])
}
}
c<-temp/count
mean(c)	