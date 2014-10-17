make.incidence.from.batch <- function(
    y,
    filen,
    shapes,
    fileformat,
    wksUsed,
    noff,
    distsInOrder,
    realsUsed
) {
  
  # Load different file formats
  if (fileformat=="id_spatial_sim") {
    datSim <- read.table(file=filen,header=TRUE)
    datSim <- datSim[datSim$Event==0,]
    datSim$X <- datSim$X * 180 / pi
    datSim$Y <- datSim$Y * 180 / pi
  } else if (fileformat=="EbolaSim") {
    datSim <- read.csv(
        file=filen,header=FALSE,
        col.names=c("Run","Day","Index","X","Y","t_infector","infector"))    
  } else stop("fileformat not recognised")
  
  datSim <- datSim[datSim$Run %in% realsUsed,]
  
  datSim <- cbind(datSim,dist_code=overlay(SpatialPoints(cbind(lon=datSim$X,lat=datSim$Y)),shapes))
  datSim <- cbind(datSim,district=shapes$ADM2_NAME[datSim$dist_code])
  datSim$EpiWeek <- ceiling(datSim$Day/7.0)
  
  noReals <- length(realsUsed)
  
  # Define the returnvalues
  rtnSumSq <- array(dim=c(noff,noReals))
  rtnAllInc <- array(
      dim=c(dim(y)[1],dim(y)[2],noff,noReals),
      dimnames=list(rownames(y),colnames(y),1:noff,1:noReals)
  )
  
  # Search through realizations and offsets
  for (i in 1:length(realsUsed)) {
    dat_one_r <- datSim[datSim$Run==realsUsed[i],]
    tmp <- make.incidence.from.linelist(
        distsLatOrder,
        dat_one_r$district,
        dat_one_r$EpiWeek,
        DTs=((min(wksUsed)-noff)):(max(wksUsed)+noff))
    for (k in 1:(noff)) {
      x <- tmp$inctab[k:(k+length(wksUsed)-1),]
      rtnAllInc[,,k,i] <- x 
      rtnSumSq[k,i] <- sum((x/sum(x)*sum(y)-y)^2)    
    }
  }
  
  list(inc=rtnAllInc,stats=rtnSumSq)
  
}

