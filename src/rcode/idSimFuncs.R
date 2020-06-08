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


both.points.in.extent <- function(x1,y1,x2,y2,ext) {
  rtn <- ifelse (	
      x1 >= ext@xmin &
          x1 <  ext@xmax &
          y1 >= ext@ymin &
          y1 <  ext@ymax &
          x2 >= ext@xmin &
          x2 <  ext@xmax &
          y2 >= ext@ymin &
          y2 <  ext@ymax, TRUE, FALSE)
  rtn
}

util_chart_pos <- function(	xindex,yindex,xn,yn,xlm=0,xrm=0,
		xg=0.01,ybm=0,ytm=0,yg=0.01) {
	
	x_left_mar = xlm
	x_right_mar = xrm
	x_gap = xg
	x_n_charts = xn
	x_width = (1-x_left_mar-x_right_mar-x_gap*(x_n_charts-1))/x_n_charts
	x_charts_l=rep(0,x_n_charts)
	for (i in 1:x_n_charts) x_charts_l[i]=x_left_mar+(i-1)*x_width+(i-1)*x_gap
	x_charts_r=rep(0,x_n_charts)
	for (i in 1:x_n_charts) x_charts_r[i]=x_charts_l[i]+x_width
	
	y_top_mar = ytm
	y_bottom_mar = ybm
	y_gap = yg
	y_n_charts = yn
	y_height = (1-y_top_mar-y_bottom_mar-y_gap*(y_n_charts-1))/y_n_charts
	y_charts_b = rep(0,y_n_charts)
	for (i in 1:y_n_charts) y_charts_b[i] = y_bottom_mar+(i-1)*y_height+(i-1)*y_gap
	y_charts_t = rep(0,y_n_charts)
	for (i in 1:y_n_charts) y_charts_t[i] = y_charts_b[i]+y_height
	
	rtn <- c(x_charts_l[xindex],x_charts_r[xindex],y_charts_b[yindex],y_charts_t[yindex])
	
	rtn
	
}

util_eventImage <- function(data,popgrid,ev,st,et,sr,er) {
	noEvents = dim(data)[1]
	rtnval = popgrid
	rtnval$z[]=NA
	xmin = popgrid$x[1]
	xgap = popgrid$x[2]-xmin
	nx = length(popgrid$x)-1
	xmax = xmin+nx*xgap
	ymin = popgrid$y[1]
	ygap = popgrid$y[2]-ymin
	ny = length(popgrid$y)-1
	ymax = ymin + ny*ygap
	for (i in 1:noEvents) {
		if (data$Event[i]==ev && data$Day[i] >= st && data$Day[i] <= et && data$Run[i] >= sr && data$Run[i] <= er) {
			xReal = data$X[i]
			yReal = data$Y[i]
			if (xReal >= xmin && xReal < xmax && yReal >= ymin && yReal < ymax) {
				xcoord = (xReal-xmin)/xgap+1
				# ycoord = ny - (yReal-ymin)/ygap+1
				ycoord = (yReal-ymin)/ygap+1
				if (is.na(rtnval$z[xcoord,ycoord])) rtnval$z[xcoord,ycoord]=1
				else rtnval$z[xcoord,ycoord]=rtnval$z[xcoord,ycoord]+1
			}
		}
	}
	noRuns = er-sr+1
	rtnval$z=rtnval$z/noRuns
	rtnval	
}

load.plot.commute <- function("fnCommutes",plot=TRUE) {
  
}

