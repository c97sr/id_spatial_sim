# Try fitting a simple skyspline model to the baseline sim in epiHistory2tree0
source('epiHistory2tree0.R')

require(skyspline) 

#~ ebovtab <- read.table( '../../scenarios/ebola/output/ebola_pset_0_Events.out', header = T )
ebovtab <- read.table( '../../epiHistory2tree/ebola_pset_0_Events.out', header = T )

## Call to observation function
obsebovtab <- obsprocess(ebovtab,obsevents=c(0),rateSymp=0.5,weekCap=20)
obsebovtab2 <- obsebovtab[ obsebovtab$actObs>0 , ]
# plot cumulative infections with or without obs process
X11(); with( obsebovtab, plot( sort(Day), 1:length(Day), log = 'y' ) )
X11(); with( obsebovtab2, plot( sort(Day), 1:length(Day), log = 'y', main='obs process' ) )

ebovtab$Index <- as.character(ebovtab$Index)
ebovtab$infector <- as.character( ebovtab$infector )
inf_ebovtab <- ebovtab[ ebovtab$Event==0, ]
ao_ebovtab <- ebovtab[ !ebovtab$Event==0, ]

donor <- inf_ebovtab$infector
pids = recip <- inf_ebovtab$Index
donor[ !(donor %in% pids )] <- 'src'
t <- inf_ebovtab$Day
tremoval <- sapply( pids, function(pid) {
	x <- ao_ebovtab$Day[ ao_ebovtab$Index==pid]
	if (any(!is.na(x))) return( max(x , na.rm=T) ) #TODO some vals too small
	NA
})
tremoval[is.na(tremoval)] <- max( tremoval, na.rm=T)

tr_tab <- data.frame( donor = donor, recip = recip, t = t)
rem_tab <- data.frame( pids = pids, tremoval = tremoval)

tre <- epi2tree( tr_tab, rem_tab )

## downsample to observed cases
rem_ebovtab <- ebovtab[ ebovtab$Event %in% c(2,6 ), ]
pidstokeep <- intersect( obsebovtab2$Index, rem_ebovtab$Index ) #only keep pids that have a removal 
obs_tips <- paste0( pidstokeep, '_0')
tre2 <- drop.tip( tre, setdiff( tre$tip.label,obs_tips ) )

## plot
plot( tre2, show.tip=F, no.mar=T)


## skyspline

n <- length(tre2$tip.label)
sts <- setNames( node.depth.edgelength( tre2 )[1:n], tre2$tip.label)


# com12 with lognormal death rate prior
if(T)
{
	bdt <- DatedTree( tre2, sts , minEdgeLength = .5, tol = 2)
	fit5 <- fit.skyspline.ml( bdt , death_rate_guess = 1/16, R0guess = 1.5, maxHeight=85, y0_guess=50, t0 = max(sts) - bdt$maxH 
	 , priors = list(death_rate = function(x)  dlnorm( x, log(1/16), 1, log=T)) )
	print(fit5$par )
	print( exp(fit5$par ) )
	print(head(fit5$demo.history))
	print(tail(fit5$demo.history))
	pbfit5 <- parboot.skyspline.mle(fit5, nreps = 100, maxHeight=89,t0 = max(sts) - bdt$maxH )
	#pbfit5 <- parboot.skyspline.mle(fit5, nreps = 10, maxHeight=89,t0 = max(sts) - bdt$maxH )
	
	plot.pbfit( pbfit5 , type = 'cumulative') -> cpl 
	#cpl <- cpl + geom_hline( aes(yintercept = nrow( obsebovtab ))
	cpl.df <- with(obsebovtab, data.frame(x = sort(Day), y = 1:length(Day)))
	cpl <- cpl + geom_point( aes(x=x,y=y), data= cpl.df )
	ggsave('testEpiHistory2tree0.0.cumulative.png', cpl)
	
	ypl <- plot.pbfit( pbfit5, type = 'size')
	
	ggsave('testEpiHistory2tree0.0.size.png', ypl)
}


if (F)
{ # generation time from sims
pids <- unique( ebovtab$Index )
xtab <- ebovtab[ ebovtab$Event %in% c(0,2,6 ),]
xtab$Event[ xtab$Event %in% c(2,6) ] <- -2
#~ x <- sapply( pids, function(pid) diff( range(xtab$Day[ebovtab$Index==pid ]) ))
x <- sapply( pids, function(pid) {
	xtab2 <- xtab[ xtab$Index==pid , ]
	if ( !( -2 %in% xtab2$Event ) ) return(NA)
	xtab2$Day[ xtab2$Event == -2 ] - xtab2$Day[ xtab2$Event == 0 ]
})
x <- x[ !is.na(x) ]
#~ > summary( x ) 
#~    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#~    2.00   11.00   16.00   17.85   23.00   70.00 
#~ >
}
