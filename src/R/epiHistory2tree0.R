require(ape) 


epi2tree <- function( transmission_table, removal_table )
{
#  vectors donor, recip, t same length
# vectors pids, tremoval same length; all pids unique, matches donor and recip
# should make this two data frames
# all recip has exactly one donor
# earliest donor should be src
#~ 	donor, recip, t
	attach( transmission_table)
	donor[is.na(donor)] <- 'src'
#~ 	pids, tremoval
	attach( removal_table )
	
	#TODO should validate input
	
	N <- length(donor)
	Nnode <- N - 1 + 1 #include source +1
	edge <- matrix( NA, nrow = N + Nnode, ncol = 2)
	edge.length <- rep(NA, N + Nnode)
	recip2donor <- setNames( recip, donor )
	
	i <- order(t) 
	donor <- donor[i] 
	recip <- recip[i]
	t <- t[i] 
	
	recip2donorNode <- setNames( rep(NA, length(pids)), pids)
	recip2donorNode[ recip2donor[pids]=='src' ] <- 'src'
	
	names(tremoval) <- pids 
	
	# add transm edges
	donor2trCounter <- setNames(rep(0, length(pids)+1), c('src', pids) )
	donor2ntr <- sapply( pids, function(pid) sum( donor== pid ))
	node2time <- list()
	node2time[['src']] <- t[1]-1
	k <- 1
	for ( i in 1:length(t)){
		u <- donor[i]
		v <- recip[i]
		donor2trCounter[u] <- donor2trCounter[u] + 1
		trcount <- donor2trCounter[u] 
		if (u == 'src'){
			donornode <- u
		} else{
			donornode <- paste(u, trcount,  sep='_')
		}
		recip2donorNode[v] <- donornode
		if (u != 'src'){
			edge[k,2] <- donornode
			node2time[[ donornode ]] <- t[i] 
			if (trcount==1){
				edge[k, 1]  <- 	recip2donorNode[u] 
			} else{
				edge[k,1] <- paste(u, trcount-1, sep='_')
			}
			k <- k + 1
		}	
	}
	
	# add terminals
	tip.label <- paste0( pids, '_0' ) # corresponding to death / removal
	for ( i in 1:length( pids )){
		pid <- pids[i] 
		tl <- tip.label[i]
		node2time[[tl]] <- tremoval[i] 
		if ( donor2trCounter[pid] == 0 ){
			lastnode <- recip2donorNode[pid]
		} else{
			trcount <- donor2trCounter[pid] 
			lastnode <- paste(pid, trcount,  sep='_')
		}
		edge[k,1] <- lastnode
		edge[k,2] <- tl
		k <- k + 1
	}
	
	internalNodes <- unique( as.vector( edge ) )
	internalNodes <- internalNodes[ !grepl('_0', internalNodes ) ]
	
	i_internalNodes <- setNames( N + 1:length(internalNodes), internalNodes )
	i_tips <- setNames( 1:N, tip.label)
	nodemap <- c( i_internalNodes, i_tips)
	
	edge <- edge[!is.na(edge[,1]), ]
	edge2 <- matrix(NA, nrow =nrow(edge), ncol =ncol(edge))
	edge2[,1] <- nodemap[edge[,1] ]
	edge2[,2] <- nodemap[edge[,2] ]
	
	edge.length <- rep(NA, nrow(edge2))
	edge.length <-  unname( unlist( node2time[edge[,2]] ) - unlist(  node2time[ edge[,1] ] ) )
	
	o <- list( edge = edge2, tip.label = tip.label, edge.length = edge.length, Nnode = Nnode )
	class(o) <- 'phylo'
	multi2di( read.tree(text= write.tree( o )) ) 
}


if (F)
{
	ebovtab <- read.table( 'ebola_pset_0_Events.out', header = T )
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
		if (any(!is.na(x))) return( max(x  , na.rm=T) ) #TODO some vals too small
		NA
	})
	tremoval[is.na(tremoval)] <- max( tremoval, na.rm=T)
	
	tr_tab <- data.frame( donor = donor, recip = recip, t = t)
	rem_tab <- data.frame( pids = pids, tremoval = tremoval)
	
	tre <- epi2tree( tr_tab, rem_tab )
	
	#downsample to x%
	tre <- drop.tip( tre, sample( tre$tip.label, floor( .9*length(tre$tip.label)), replace=F) )
	
	# plot
	plot( tre, show.tip=F)
	
	# quick phylodynamic
	require(phylodyn) # Karcher + suchard effective population size estimator 
	BNPR( tre ) -> o
	plot_BNPR( o , ylim = c(200, 2e5 ))
}
