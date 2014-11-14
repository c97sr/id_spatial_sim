pdf.figure.proof <-
    function(findex=1,file=paste("./pdf_sub_",findex,".pdf",sep=""),pw=7,ph=7,textscale=0.6,xpd=NA) {
  plotwidth <- pw/cm(1)
  plotheight <- ph/cm(1)
  margin <- 5/cm(1)
  pdfwidth <- plotwidth+2*margin
  pdfheight <- plotheight+2*margin
  posplot = c(
      margin/pdfwidth,
      (margin+plotwidth)/pdfwidth,
      margin/pdfheight,
      (margin+plotheight)/pdfheight
  )
  pdf(file=file,height=pdfheight,width=pdfwidth,useDingbats=FALSE)
  par(
      cex = textscale,
      fig = posplot,
      mai = c(0,0,0,0),
      xpd = xpd)
  posplot
}

ebola.detrend <- function(incTab,df=0) {
  
  # Make some detrended data
  rtn_fit <- as.matrix(incTab)
  rtn_marg <- as.matrix(incTab)
  dists <- colnames(incTab)
  tsnames <- rownames(incTab)
  nots <- dim(incTab)[1]
  rtn_fit[] <- 0
  rtn_marg[] <- 0
  colnames(rtn_fit) <- dists
  colnames(rtn_marg) <- dists
  rownames(rtn_fit) <- tsnames
  rownames(rtn_marg) <- tsnames
  for (d in dists) {
    xvals <- 1:nots
    yvals <- incTab[,d]
    if (df==0) {
      trend <- rep(min(yvals),nots)
    } else if (df==1) {
      trend <- rep(mean(yvals),nots)
    } else {
      trend_tmp <- (smooth.spline(xvals,yvals,df=df))
      trend <- as.numeric(trend_tmp$y)
    }
    rtn_fit[,d] <- trend
    rtn_marg[,d] <- as.numeric(yvals) - trend
  }
  
  list(mod=rtn_fit,marg=rtn_marg)
  
}

ebola.pdf.cor.fig <- function(tab,file=NULL,pw=7,ph=7,textscale=10/12) {
  if (is.null(file)) error("File name must be specified")
  pdf.figure.proof(file=file,pw=pw,ph=ph,textscale=textscale)
  corlistlog <- rcorr(tab,type="pearson")
  image(corlistlog$r,axes=FALSE,col=cols,breaks=breaks)
  ebola.dist.axis(1)
  ebola.dist.axis(2)
  for (i in 1:nodists) {
    for (j in 1:nodists) {
      if (i != j && corlistlog$P[i,j] < 0.05/ ((nodists^2-nodists)/2)) {
        points((i-1)/(nodists-1),(j-1)/(nodists-1),pch=21)
      }
    }
  }
  dev.off()  
}

# Make an incidence table
make.incidence.from.linelist <- function(
    vecAllDists,
    vecIncDists,
    vecTimePeriod,
    DTs = min(vecTimePeriod):max(vecTimePeriod)) {
  
  early <- 0
  late <- 0
  outside <- 0
  firstDT <- min(DTs)
  lastDT <- max(DTs)
  noDT <- length(DTs)
  rtn <- matrix(data=0,nrow=noDT,ncol=length(vecAllDists),dimnames=list(firstDT:(lastDT),vecAllDists))
  rowcounter <- 1
  maxrow <- length(vecIncDists)  
  
  # This is the section that needs optimizing
  # Work on the server, sort twice, first by day and then by district and do a single sweep through
  # Should be able to reduce the processing time vastly
  while (rowcounter < maxrow) {
    currentdist <- vecIncDists[rowcounter]
    if (currentdist %in% vecAllDists) {
      currentweek <- as.numeric(vecTimePeriod[rowcounter])
      if (currentweek < firstDT) {
        early <- early + 1
      } else if (currentweek > lastDT) {
        late <- late + 1
      } else {
        rtn[currentweek-firstDT+1,as.character(currentdist)] <- 
            rtn[currentweek-firstDT+1,as.character(currentdist)] + 1
      }
    } else {
      outside <- outside + 1
    }
    rowcounter <- rowcounter+1
  }
  
  # Return the final values
  list(inctab=rtn,early=early,late=late,outside=outside)
}

# Function to output a heat chart of districts with legend
# To be run in two different ways
inc.heat.chart.pdf.v2 <- function(
    z,
    vecCountries,
    cols = c(colors()[1],rev(heat.colors(200)[1:180])),
    breaks=c(0,0.5,1.5,1.5+(1:(length(cols)-2))/(length(cols)-2)*zmax),
    zmax=ceiling(max(z)/100)*100,
    filen="district_inc",
    xlabs=c(1,seq(20,40,5)),
    chttype="",
    pdfWidth=10,
    leftMarginFraction=0.25,
    zmaxabsolute=200,
    outstem="./figs/",
    legendlabs=NULL,
    legendats=NULL) {
  
  # Define local functions
  ebola.dist.axis.v2 <- function(side,vecDists,vecCountries,vecColCountries) {
    noDists <- length(vecDists)
    labellocs <- (0:(noDists-1))/(noDists-1)
    reps <- names(table(vecCountries))
    axis(side,at=c(0,1),labels=c("",""))
    for (ctry in reps) {
      axis(	side,at=labellocs[vecCountries==ctry],
          labels=vecDists[vecCountries==ctry],
          col.axis=vecColCountries[ctry],
          las=2)
    }
  }
  
  # Test for log values
  col.countries <- c("GUINEA"="#4f81bd","SIERRA LEONE"="#c0504d", "LIBERIA"="#9bbb59","NIGERIA"="#a6a6a6")
  
  # Define local aux variables
  minDT <- min(as.numeric(rownames(z)))
  maxDT <- max(as.numeric(rownames(z)))
  vecDists <- colnames(z)
  
  noDists <- length(vecDists)
  useCountries <- TRUE
  if(length(vecCountries)==0) {
    useCountries <- FALSE
  } else if (length(vecCountries)!=noDists) {
    stop("Vector of country names not the same as vector of district names")
  }
  if (useCountries) {
    tabCountries <- table(vecCountries)
  }
  
  # Open the device and set the usual parameters
  # pdf(file=paste(outstem,filen,chttype,".pdf",sep=""),height=15,width=pdfWidth,useDingbats=FALSE)
  png(file=paste(outstem,filen,chttype,".png",sep=""),height=15,width=pdfWidth,units="in",res=300)
  par(	cex=12/12,
      mai=c(0,0,0,0),
      xpd=FALSE
  )
  
  # Check that the full scale is being drawn (otherwise might get white instead of high val)
  if (max(breaks) < max(z)) stop("Full range of z not plotted")
  
  # Plot the main chart
  par(fig=c(leftMarginFraction,0.8,0.075,0.975),new=FALSE)
  image(z,breaks=breaks,axes=FALSE,col=cols,xlab="",ylab="")
  xaxilabels <- xlabs 
  axis(1,at=(xaxilabels-minDT)/(maxDT-minDT),labels=xaxilabels)
  mtext("Epidemiological week",1,line=3)
  if (useCountries) {
    mtext("Liberia",side=2,line=10.5,at=tabCountries["LIBERIA"]/2/noDists)
    mtext("Sierra Leone",side=2,line=10.5,at=(tabCountries["LIBERIA"]+tabCountries["SIERRA LEONE"]/2)/noDists)
    mtext("Guinea",side=2,line=10.5,at=(tabCountries["LIBERIA"]+tabCountries["SIERRA LEONE"]+tabCountries["GUINEA"]/2)/noDists)
  } 
  ebola.dist.axis.v2(2,vecDists,vecCountries,col.countries)	
  
  # Plot the legend
  par(fig=c(0.9,0.95,0.5,0.9),new=TRUE)
  
  plot(1:2,type="n",
      axes=FALSE,col=cols,
      ylab="Incidence",xlab="",
      ylim=c(min(breaks),max(breaks)),
      xlim=c(0,1)
  )
  for (i in 1:length(cols)) {
    polygon(c(0,1,1,0),c(breaks[i],
            breaks[i],breaks[i+1],breaks[i+1]),border=NA,col=cols[i])
  }
  
#  image(1,breaks,
#      t(as.matrix(breaks)),
#      axes=FALSE,col=cols,
#      ylab="Incidence",xlab="",
#      ylim=c(min(breaks),max(breaks)),
#      oldstyle=FALSE)
  if (is.null(legendlabs)) {
    axis(2,las=1)
  } else {
    axis(2,las=1,labels=legendlabs,at=legendats)
  }
  mtext("incidence",side=3,line=0.5)
  if (chttype=="_absolute") {
    mtext("Absolute",side=3,line=2.5)
    mtext("weekly",side=3,line=1.5)
  } else if (chttype=="_percap") {
    mtext("Per capita",side=3,line=2.5)
    mtext("weekly",side=3,line=1.5)
  } else {
    mtext("Weekly",side=3,line=1.5)
  }
  
  # Close the pdf device
  dev.off()
}

spatial.prune.v2 <- function(matDat,firstDT,casesinc, useEpiWeek=TRUE) {
  
  # Prune data for:Nigerian districts, EpiCaseDef and NAs in date
  nigerian_districts <- c("ALIMOSHO","ETI OSA","IBEJU LEKKI","KOSOFE","OSHODI/ISOLO","PORT-HARCOURT","SURULERE")
  rtn <- matDat[!(matDat$district %in% nigerian_districts),]
  rtn <- rtn[rtn$EpiCaseDef %in% casesinc,]
  rtn <- rtn[!is.na(rtn$DateOnsetInferred),]
  
  # Make week variable	
  if (useEpiWeek) {
    swone <- as.Date("30/12/13", "%d/%m/%y")
  } else {
    swone <- min(rtn$DateOnsetInferred)
  }
  rtn$EpiWeek <- as.numeric(ceiling((rtn$DateOnsetInferred - swone)/7))
  rtn$EpiDay <- as.numeric((rtn$DateOnsetInferred - swone))
  
  # Trim weeks used as required and make utility variables for bounds
  rtn <- rtn[rtn$EpiWeek >= firstDT,]
  rtn
  
}


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

testInlineCpp <- function() {
  
  library(inline)  # use version 0.3.10 for rcpp() wrapper
  addMat <- rcpp(signature(ms="numeric"), body='
          +    Rcpp::NumericMatrix M(ms);
          +    Rcpp::NumericMatrix N = Rcpp::clone(M);
          +    for (int i=0; i<M.nrow(); i++)
          +       for (int j=0; j<M.ncol(); j++) 
          +          N(i,j) = M(i,j) + 1;
          +    return(N);
          + ')
  
}

    # R> addMat(matrix(1:9,3,3))
# [,1] [,2] [,3]
# [1,]    2    5    8
# [2,]    3    6    9
# [3,]    4    7   10

# Function to take output of the spatial simulation and make it comparable with
# WHO data.
# tabss - the raw csv read in from spatial sim output
convert.spatialsim.WHO.V1 <- function(tabss,distNames,tabLookup,max_week) {
  
  # Extract the correct names
  #  ssNames <- colnames(tabss)
  #  ssNames <- sub("C","",ssNames)
  #  ssNames <- ifelse(is.na(match(ssNames,tabLookup$Code)),ssNames,tabLookup$WHO_districts[match(ssNames,tabLookup$Code)])
  #  setnames(tabss,names(tabss),ssNames)
  
  # Re-prders the tabss names to be in the distNames order
  rtn_day <- as.matrix(subset(tabss,select=distNames))
  
  # Setup the return type
  rtn_week <- matrix(
      nrow=max_week,ncol=dim(rtn_day)[2],
      dimnames=list(1:max_week,colnames(rtn_day)))
  
  rtn_week[] <- 0
  nodays <- dim(rtn_day)[1]
  
  current_day <- 1
  current_week <- 1
  day_next_week <- 8
  while (current_week <= max_week && current_day <= nodays) {    
    rtn_week[current_week,] <- rtn_week[current_week,] + rtn_day[current_day,]
    current_day <- current_day + 1
    if (current_day == day_next_week) {
      current_week <- current_week + 1
      day_next_week <- day_next_week + 7
    }
  } 
  
  rtn_week
  
}

# Function to output a heat chart of districts with legend
# To be run in two different ways
inc.heat.chart.pdf.v3 <- function(
    z,
    vecCountries,
    cols = c(colors()[1],rev(heat.colors(200)[1:180])),
    breaks=c(0,0.5,1.5,1.5+(1:(length(cols)-2))/(length(cols)-2)*zmax),
    zmax=ceiling(max(z)/100)*100,
    filen="district_inc",
    xlabs=c(1,seq(20,40,5)),
    chttype="",
    pdfWidth=10,
    leftMarginFraction=0.25,
    zmaxabsolute=200,
    outstem="./figs/",
    legendlabs=NULL,
    legendats=NULL,
    addCountryBreak=FALSE) {
  
  # Add a break in for the first two changes in country
  if (addCountryBreak) {
  }
  
  # Define local functions
  ebola.dist.axis.v2 <- function(side,vecDists,vecCountries,vecColCountries) {
    noDists <- length(vecDists)
    labellocs <- (0:(noDists-1))/(noDists-1)
    reps <- names(table(vecCountries))
    axis(side,at=c(0,1),labels=c("",""))
    for (ctry in reps) {
      axis(	side,at=labellocs[vecCountries==ctry],
          labels=vecDists[vecCountries==ctry],
          col.axis=vecColCountries[ctry],
          las=2)
    }
  }
  
  # Test for log values
  col.countries <- c("GUINEA"="#4f81bd","SIERRA LEONE"="#c0504d", "LIBERIA"="#9bbb59","NIGERIA"="#a6a6a6")
  
  # Define local aux variables
  minDT <- min(as.numeric(rownames(z)))
  maxDT <- max(as.numeric(rownames(z)))
  vecDists <- colnames(z)
  
  noDists <- length(vecDists)
  useCountries <- TRUE
  if(length(vecCountries)==0) {
    useCountries <- FALSE
  } else if (length(vecCountries)!=noDists) {
    stop("Vector of country names not the same as vector of district names")
  }
  if (useCountries) {
    tabCountries <- table(vecCountries)
  }
  
  # Open the device and set the usual parameters
  # pdf(file=paste(outstem,filen,chttype,".pdf",sep=""),height=15,width=pdfWidth,useDingbats=FALSE)
  png(file=paste(outstem,filen,chttype,".png",sep=""),height=15,width=pdfWidth,units="in",res=300)
  par(	cex=12/12,
      mai=c(0,0,0,0),
      xpd=FALSE
  )
  
  # Check that the full scale is being drawn (otherwise might get white instead of high val)
  if (max(breaks) < max(z)) stop("Full range of z not plotted")
  
  # Plot the main chart
  par(fig=c(leftMarginFraction,0.8,0.075,0.975),new=FALSE)
  image(z,breaks=breaks,axes=FALSE,col=cols,xlab="",ylab="")
  xaxilabels <- xlabs 
  axis(1,at=(xaxilabels-minDT)/(maxDT-minDT),labels=xaxilabels)
  mtext("Epidemiological week",1,line=3)
  if (useCountries) {
    mtext("Liberia",side=2,line=10.5,at=tabCountries["LIBERIA"]/2/noDists)
    mtext("Sierra Leone",side=2,line=10.5,at=(tabCountries["LIBERIA"]+tabCountries["SIERRA LEONE"]/2)/noDists)
    mtext("Guinea",side=2,line=10.5,at=(tabCountries["LIBERIA"]+tabCountries["SIERRA LEONE"]+tabCountries["GUINEA"]/2)/noDists)
  } 
  ebola.dist.axis.v2(2,vecDists,vecCountries,col.countries)	
  
  # Plot the legend
  par(fig=c(0.9,0.95,0.5,0.9),new=TRUE)
  
  plot(1:2,type="n",
      axes=FALSE,col=cols,
      ylab="Incidence",xlab="",
      ylim=c(min(breaks),max(breaks)),
      xlim=c(0,1)
  )
  for (i in 1:length(cols)) {
    polygon(c(0,1,1,0),c(breaks[i],
            breaks[i],breaks[i+1],breaks[i+1]),border=NA,col=cols[i])
  }
  
#  image(1,breaks,
#      t(as.matrix(breaks)),
#      axes=FALSE,col=cols,
#      ylab="Incidence",xlab="",
#      ylim=c(min(breaks),max(breaks)),
#      oldstyle=FALSE)
  if (is.null(legendlabs)) {
    axis(2,las=1)
  } else {
    axis(2,las=1,labels=legendlabs,at=legendats)
  }
  mtext("incidence",side=3,line=0.5)
  if (chttype=="_absolute") {
    mtext("Absolute",side=3,line=2.5)
    mtext("weekly",side=3,line=1.5)
  } else if (chttype=="_percap") {
    mtext("Per capita",side=3,line=2.5)
    mtext("weekly",side=3,line=1.5)
  } else {
    mtext("Weekly",side=3,line=1.5)
  }
  
  # Close the pdf device
  dev.off()
}


