RPdistributions <- function(Node_connection_within_household, 
                            Node_connection_between_household,
                            fname){
  
  #Get HH>1 locations from within HH network:
  W <- Node_connection_within_household
  W1 <- W[, c("V1", "V2", "V3")]
  W2 <- W[, c("V4", "V5", "V6")]
  colnames(W1) <- c("node", "long", "lat")
  colnames(W2) <- c("node", "long", "lat")
  net <- rbind(W1, W2)
  net <- net[!duplicated(net$node), ]
  longs <- matrix(NA, max(net$node), 1)
  longs[net$node, ] <- net$long
  
  #Workplace network
  B <- Node_connection_between_household
  colnames(B) <- c("nodei", "longi", "lati","nodej", "longj", "latj")
  
  #B1 <- B[,c(1,2,3)]
  #B2 <- B[,c(4,5,6)]
  #Replace longj/all with value taken from net
  B$longi <- net$long[B$nodei]
  B$lati <- net$lat[B$nodei]
  B$longj <- net$long[B$nodej]
  B$lati <- net$lat[B$nodei]
  
  dlat <- abs(B$latj-B$latj)*111
  dlong <- abs(B$longi-B$longj)*110
  distances <- sqrt(dlat^2+dlong^2)
  h <- hist(distances, breaks=seq(0,1,.02),xlim=c(0,1),ylim=c(0,6*10^5),main=fname)
}

################################################################

RPHHnet <- function(Node_connection_within_household, 
                    Node_connection_between_household,
                    fname){
  #Dissolve HH intp table of HHs and members:
  H <- Node_connection_within_household
  
  alpha <- 0
  alpha <-1#"Within" in wrong format - some missing due to no wp connections
  if (alpha==1){
    B <- Node_connection_between_household
    colnames(B) <- c("V1", "V2", "V3","V4", "V5", "V6")
    B1 <- B[, c("V1", "V2", "V3")]
    B2 <- B[, c("V4", "V5", "V6")]
    colnames(B1) <- c("node", "long", "lat")
    colnames(B2) <- c("node", "long", "lat")
    net <- rbind(B1, B2)
    net <- net[!duplicated(net$node), ]
    numnodes <- max(net$node)
    
    Hdata <- data.frame(seq(1,numnodes))
    Hdata$lat[net$node] <- net$long
    Hdata$lat[net$node] <- net$lat
    
    longs <- matrix(NA, max(net$node), 1)
    longs[net$node, ] <- net$long
    #Workplace network
    B <- Node_connection_between_household
    colnames(B) <- c("nodei", "longi", "lati","nodej", "longj", "latj")
    #Replace longj/all with value taken from net
    B$longi <- net$long[B$nodei]
    B$lati <- net$lat[B$nodei]
    B$longj <- net$long[B$nodej]
    B$lati <- net$lat[B$nodei]
  }
  
  H1 <- H[,1:3]
  H2 <- H[,4:6]
  colnames(H1) <- c('node', 'c1', 'c2')
  colnames(H2) <- c('node', 'c1', 'c2')
  Hfull <- rbind(H1,H2)
  Hlocs <- unique(Hfull[, c('c1', 'c2')])
  Hlocs[, 'hh'] <- 1:nrow(Hlocs)
  
  Hall <- merge(Hfull, Hlocs, by=c('c1','c2'))
  Hindex <- Hall[, c('node', 'hh')]
  Hindex <- unique(Hindex)#node and household index
  
  maxind <- max(Hindex$node)
  v1 <- seq(1, maxind, 1)
  v2 <- rep(0, maxind)
  v <- data.frame(v1, v2)
  colnames(v) <- c('node', 'household')
  v$household[Hindex$node] <- Hindex$hh
  maxhh <- max(Hindex$hh)
  #Take care of hh size 1:
  hh1 <- which(v$household==0)
  lh1 <- length(hh1)
  if (lh1>0) {
    v$household[hh1] <- seq(maxhh+1, maxhh+lh1, 1)
  }
  
  #Find links between HHs:
  W <- Node_connection_between_household
  W1 <- W[, 1:3]
  W2 <- W[, 4:6]
  colnames(W1) <- c('node', 'c1', 'c2')
  colnames(W2) <- c('node', 'c1', 'c2')
  
  X <- W[, c(1, 4)]
  X[, 1] <- v$household[X[, 1]]
  X[, 2] <- v$household[X[, 2]]
  
  write.csv(X, fname)
}

################################################################

RPcombine <- function(Node_connection_within_household, 
                      Node_connection_between_household){
  require(Matrix)
  #
  within <- Node_connection_within_household[, c("V1", "V4")]
  between <- Node_connection_between_household[, c("I", "J")]
  colnames(within) <- c("I", "J")
  net <- rbind(within, between)
  #
  #Option to write CSV with (i,j) pairs:
  #write.csv(net, "gos5p5.csv")
  #
  #Create matrix
  #net <- net[-which(net[, 1]<net[, 2]), ]
  numNodes <- max(net)
  
  G <- sparseMatrix(net$I, net$J, x=1, dims=c(numNodes, numNodes), symmetric = TRUE)
  
  #G <- cbind(c(net$I, net$J),c(net$J, net$i))
  #colnames(G) <- c("I", "J")
  #G %>% distinct()
  RPcombine <- G
}

################################################################

RPmcluster <- function(G,N){
  #Change/option to loop through N and output hist/stats for each value
  G@x[G@x>1] <- 1
  n <- nrow(G)
  #eps <- .01
  #bins <- seq(0, eps, 1)
  measureA <- matrix(list(), N, 1)
  measureB <- measureA
  measureA[[1, 1]] <- G
  measureB[[1, 1]] <- G
  A <- G
  B <- G
  X <- matrix(0, n, 1)
  #
  if(N==1){
    #Standard clustering:
    for (i in 1:n){
      vi <- A[i, ]
      vfind <- which(vi==1)
      lv <- length(vfind)
      Wi <- B[vfind, vfind]
      X[i, 1] <- sum(upper.tri(Wi, diag=FALSE))/(lv*(lv+1)/2)
    }
  }else{
    #N-clustering, N>=2
    for (i in 2:N){
      nextA <- A%*%G
      #nextA <- subset(nextA, select=-diag(nextA))
      diag(nextA) <- 0
      nextA@x[B@x==1] <- 0
      nextA@x[nextA@x>1] <- 1
      nextB <- B+nextA
      nextB@x[nextB@x>1] <- 1
      measureA[[i, 1]] <- nextA
      measureB[[i, 1]] <- nextB
      A <- nextA
      B <- nextB
    }
    Ax <- measureA[[N, 1]]
    #subset(Ax, select=-diag(Ax))
    diag(Ax) <- 0
    Bx <- measureB[[N, 1]]
    #subset(Bx, select=-diag(Bx))
    diag(Bx) <- 0
    B1 <- measureB[[1, 1]]
    BNm1 <- measureB[[N-1, 1]]
    for (i in 1:n){#N=1 - this loop with Bx=A - can't use BNm1
      vi <- Bx[i, ]
      vfind <- which(vi==1)
      lv <- length(vfind)
      Wi <- Bx[vfind, vfind]
      X[i, 1] <- sum(upper.tri(Wi, diag=FALSE))/(lv*(lv+1)/2)
    }
  }
  #m-clustering for whole network:
  RPmcluster <- X
  #Summary stats:
  RPmclusterStats <- c(mean(X, na.rm=TRUE), var(X, na.rm=TRUE), quantile(X, c(0, .05, .25, .5, .75, .95, 1), na.rm=TRUE))
}

################################################################

#MATLAB: RPmclusterSample:
RPmclusterSampleStats <- function(G,mmax,sampleNumber){
  #G must be simple: undirected with no self loops
  N <- mmax#Just in case - variable re-named
  #G@x[G@x>1] <- 1
  n <- nrow(G)
  v <- sample(n, sampleNumber, replace=FALSE)#v is indices
  lv <- length(v)
  clusterVec1 <- matrix(0, mmax, sampleNumber)#MATLAB: transpose
  
  function(i) which(G[i,] == 1)
  GG <- apply(G, 1, function(x) which (x == 1))
  
  imvw1 <- t(v)#node
  imvw2 <- rep(1, lv)#info order
  Cv <- vector("list", n)#List of data for each node
  lengthC <- lv
  
  #Standard clustering:
  for (i in 1:sampleNumber){
    #Find neighbours of node i:
    vx <- v[i]
    v1 <- mNbr(GG, vx, 1, vx)#Indices of neighbours
    #Update stored data:
    Cv[[i]] <- v1
    #
    Gred <- G[v1, ]#Adj mat of neighbourhood
    Gred <- Gred[, t(v1)]
    lv1 <- length(v1)#Number of neighbours
    #How many links exist/how many could exist:
    clusterVec1[1, i] <- sum(upper.tri(Gred, diag=FALSE))/lv1/(lv1-1)*2
  }
  #N-clustering, N>=2
  if (mmax>1){
    for (m in 2:mmax){
      for (i in 1:lv){
        #m-neighbours of v(i):
        vi <- v[i]
        vx <- Cv[[i]]#m-1 OR m (if already xcalculated in loop) nbrs of v(i)
        imvw2i <- imvw2[i]
        if (imvw2i<m){#i.e.==m-1
          vm <- mNbr(GG, vx, (m-imvw2i), vi)#Neighbours of i
          #vm=vm(vm~=vi) here?
          imvw2[i] <- m
          Cv[[i]] <- vm
        } else{
          vm <- vx
        }
        vmx <- vm#(vm~=vi); #Exclude node i
        lvm <- length(vmx)
        #m-clustering:
        #
        links1 <- 0
        #m-neighbours of each m-neighbour:
        for (j in 1:lvm){#parfor
          vmj <- vmx[j]#vmOrder(j);
          #if vmj~=vi
          #Faster to order or to find all links?
          if (vmj %in% imvw1==TRUE){#If have already calculated some nbhds
            index <- which(imvw1==vmj)#How many
            thisFar <- imvw2[index]
            vmjSoFar <- Cv[[index]]
            if (thisFar<m){
              mNj <- mNbr(GG, vmjSoFar, m-thisFar, vmj)
              imvw2[index] <- m
              Cv[[index]] <- mNj
            } else{
              mNj <- Cv[[index]]
            }
          } else{
            mNj <- mNbr(GG, vmj, m, vmj)
            lengthC <- lengthC+1
            index <- lengthC
            imvw1[index] <- vmj
            imvw2[index] <- m
            Cv[[index]] <- mNj
          }
          mNjx <- mNj#xlinks1 <- links1+length(intersect(vmx,mNjx))
        }
        clusterVec1[m, i] <- links1/lvm/(lvm-1)
      }
    }
  }
  #m-clustering for whole sample:
  RPmclusterSample <- X
  #Summary stats:
  RPmclusterSampleStats <- c(mean(X, na.rm=TRUE), var(X, na.rm=TRUE), quantile(X, c(0, .05, .25, .5, .75, .95, 1), na.rm=TRUE))
}

################################################################

#This is needed for RPclusterSample
mNbr <- function(GG, vx1, mMore, vi){
  vx <- as.vector(vx1)
  G1 <- GG[vx]
  G1 <- unlist(G1)
  #sumG1 <- colSums(G1)
  #v1 <- which(sumG1>0)
  v1 <- unique(G1)
  v1 <- v1[v1!=vi]
  vx <- as.vector(v1)
  if (mMore>1){
    k=2
    while (k<=mMore){
      Gm <- GG[vx]
      Gm <- unlist(Gm)
      #sumGm <- colSums(Gm)
      #vm <- which(sumGm>0)
      vm <- unique(Gm)
      vm <- union(vx, vm)
      vm <- setdiff(vm, vx)
      vx <- vm
      k=k+1
    }
  }
  mNbr <- vx
}

################################################################

RPmultiple <- function(G){
  n=2
  eps <- .01
  bins <- seq(0, eps, 1)
  #Hist here
  stats <- matrix(0, n, 9)
  for (i in 1:n){
    I <- 4+i
    g1 <- RPcluster(G, I)
    stats[I, ] <- g1
  }
  RPmultiple <- stats
}

################################################################

if (FALSE){
  install.packages("Matrix")
  require(Matrix)
  source('~/networkProcessing.R', echo=TRUE)
  sampleNumber <- 200
  adjMat <- RPcombine(Node_connection_within_household, 
                      Node_connection_between_household)
  clusterStats <- RPclusterSample(adjMat, 2, sampleNumber)
}