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

