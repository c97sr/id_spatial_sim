library("matrixStats")
library(readr)

#setwd(".../id_spatial_sim/scenarios/ebola/output")
nruns <-20

# creates the list of all the events files in the directory
list = dir(pattern = "*_pset_0_Events.out") ##########Change letter here - +1 below

#Ensure files are in numerical order
list

generation_total= as.numeric()
generation_average=as.numeric()
generation_average_final=as.numeric()
daily_incidence_total=as.numeric()
output_matrix=matrix(nrow=length(list)*nruns, ncol=10)
output_mean_variance_matrix=matrix(nrow=length(list), ncol=13)

for (network in 1:length(list)){
  
  events = read_delim(list[network], 
                      "\t", escape_double = FALSE, col_names = TRUE, 
                      trim_ws = TRUE)
  
  str(events)
  
  ##For calcuating generation reproduction number
  events$Event[events$Event != 0] = NA
  events = na.omit(events)
  
  for(run in 0:(nruns-1)){
    
    subset=subset(events,Run == run)    
    
    ##For calculation of the generational ratio
    generation_per_run_table = cbind(rep(network,length(table(subset$Generation))),
                                     rep(run,length(table(subset$Generation))),
                                     c(1:length(table(subset$Generation))),
                                     table(subset$Generation))
    
    
    ratio_per_run=as.numeric()
    
    
    for(row in 1:nrow(generation_per_run_table)){
      if(row==1){
        ratio_per_run[row]=NA
      }else{
        ratio_per_run[row]=generation_per_run_table[row,4]/generation_per_run_table[row-1,4]
      }
    }  
    
    #Calculating R0 from the 4th to the 7th generation pair in each runs
    avg_R0_per_run  = sum(ratio_per_run[5:8])/4
    output_matrix[((network-1)*nruns+(run+1)),1] = network
    output_matrix[((network-1)*nruns+(run+1)),2] = run
    output_matrix[((network-1)*nruns+(run+1)),3] = avg_R0_per_run
    
    generation_per_run_table=cbind(generation_per_run_table,ratio_per_run) 
    
    #Combines generation data from all runs together
    generation_total=rbind(generation_total,generation_per_run_table)
    
    ##For obtaining the daily incidence in each run
    daily_incidence_matrix=matrix(nrow=max(subset$Day)+1,ncol=4)
    
    for(day in 0:max(subset$Day)){
      daily_incidence_matrix[(day+1),1] = network
      daily_incidence_matrix[(day+1),2] = run
      daily_incidence_matrix[(day+1),3] = day
      daily_incidence_matrix[(day+1),4] = sum(subset$Day == day)
      
    }
    
    daily_incidence_total=rbind(daily_incidence_total,daily_incidence_matrix)
    
    
    ##For obtaining epidemic properties
    ##Peak size
    output_matrix[((network-1)*nruns+(run+1)),4] = max(daily_incidence_matrix[,4])
    
    ##Time of peak
    output_matrix[((network-1)*nruns+(run+1)),5] = which.max(daily_incidence_matrix[,4])-1
    
    ##Time to extinction
    output_matrix[((network-1)*nruns+(run+1)),6] = max(subset$Day)
    
    ##Final size of epidemic = nrow(subset)
    output_matrix[((network-1)*nruns+(run+1)),7] = nrow(subset)
    
    ##To obtain the time between half of epidemic peak size
    half_max=max(daily_incidence_matrix[,4])/2
    lower_half=which.max(daily_incidence_matrix[,4])-1
    upper_half=which.max(daily_incidence_matrix[,4])-1
    
    if(lower_half == 0){
      output_matrix[((network-1)*nruns+(run+1)),8]=lower_half
      lower_half=0
    }else{
      while(lower_half > 0){
        
        lower_diff=daily_incidence_matrix[lower_half+1,4]-half_max
        if(lower_diff>0) {lower_half=lower_half-1
        } else{
          output_matrix[((network-1)*nruns+(run+1)),8]=lower_half
          lower_half=0
        }
      }  
    }
    
    while(upper_half < max(subset$Day)){
      upper_diff=daily_incidence_matrix[upper_half+1,4]-half_max
      if(upper_diff>0) {upper_half=upper_half+1
      } else{
        output_matrix[((network-1)*nruns+(run+1)),9]=upper_half
        upper_half=max(subset$Day)+1
      }
    }
    
    ##Time interval between peaks
    output_matrix[((network-1)*nruns+(run+1)),10] = output_matrix[((network-1)*nruns+(run+1)),9]-output_matrix[((network-1)*nruns+(run+1)),8]
    
    ##closeloop for run 0:9
  }
  
}

write.table(output_matrix, file ="outSAwp.csv", row.names=FALSE, ##########Change letter here - +1 above
            col.names = c("Network","Run","Avg_R0_per_run", "Peak_size",
                          "Peak_time","Time_to_extinction","Final_size",
                          "Lower_half","Upper_half","Time_interval_half_peak"), 
            sep=",")
