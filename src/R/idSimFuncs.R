make.incidence.from.batch <- function(
    y,
    filen,
    shapes,
    fileformat,
    wksUsed,
    distsInOrder,
    realsUsed,
    maxNoEvents
) {
  
  require("data.table")
  
  # Load different file formats
  if (fileformat=="id_spatial_sim") {
    datSim <- fread(input=filen,header=TRUE)
    datSim <- datSim[datSim$Event==0,]
    datSim$X <- datSim$X * 180 / pi
    datSim$Y <- datSim$Y * 180 / pi
  } else if (fileformat=="EbolaSim") {
    datSim <- fread(input=filen,header=FALSE,showProgress=FALSE)
    setnames(datSim,names(datSim),c("Run","Day","Index","X","Y","t_infector","infector"))
  } else stop("fileformat not recognised")
    
  datSim <- datSim[datSim$Run %in% realsUsed,]
  
  datSim <- cbind(datSim,dist_code=overlay(SpatialPoints(cbind(lon=datSim$X,lat=datSim$Y)),shapes))
  datSim <- cbind(datSim,district=shapes$ADM2_NAME[datSim$dist_code])
  datSim <- datSim[!is.na(datSim$district)]
  datSim$EpiWeek <- ceiling(datSim$Day/7.0)
  
  noReals <- length(realsUsed)
  
  # Define the returnvalues
  rtnAllInc <- array(
      dim=c(dim(y)[1],dim(y)[2],noReals),
      dimnames=list(rownames(y),colnames(y),1:noReals)
  )
  
  # Search through realizations and offsets
  for (i in 1:noReals) {
    dat_one_r <- datSim[datSim$Run==realsUsed[i],]
    if (dim(dat_one_r)[1] > maxNoEvents) {
      dat_one_r <- dat_one_r[1:maxNoEvents,]
    }
    
    tmp <- make.incidence.from.linelist(
        distsLatOrder,
        dat_one_r$district,
        dat_one_r$EpiWeek,
        DTs=min(wksUsed):max(wksUsed))
    rtnAllInc[,,i] <- tmp$inctab
  }
  
  # Return the incidence tables
  list(inc=rtnAllInc)
  
}

