# ++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# Author: Leandro Ariza-Jimenez (larizaj@eafit.edu.co) #
# ++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# Core script of the Subtractive clustering algorithm.
# This particular function doesn't compute by itself the potential of each
# data point of being a cluster center, it takes that information as input
# vector.

subclustCore <- function(X,PointPotential,Ra=0.5,sqshFactor=1.25,
                         acceptRatio=0.5,rejectRatio=0.15) {
  
  # _____________________________________________________________________
  # Normalizing dataset into a unit hyperbox (for convenience)
  # _____________________________________________________________________
  
  # Looking for the minimum and maximum values per column
  minX <- apply(X,2,min)
  maxX <- apply(X,2,max)
  
  # Verifying wheter the dataset has constant colums 
  idx_constant_col <- which(maxX == minX)
  if (length(idx_constant_col) > 0) {
    
    # The dataset has columns with constant values (so minX and maxX are equal),
    # then calculate a small range relative to the data
    minX[idx_constant_col] <- minX[idx_constant_col] - 0.0001*(1 + abs(minX[idx_constant_col]))
    maxX[idx_constant_col] <- maxX[idx_constant_col] + 0.0001*(1 + abs(maxX[idx_constant_col]))
  
    # Normalizing dataset into a unit hyperbox using modified minX and maxX
    for (i in 1:numVars) {
      X[,i] = (X[,i] - minX[i]) / (maxX[i] - minX[i]);
    }
  }
  else {
    # Normalizing dataset into a unit hyperbox using minX and maxX
    X <- apply(X, 2, function(x)(x-min(x))/(max(x)-min(x)))
  }
  X[X < 0] <- 0 # Previnting undeflow values after normalization
  X[X > 1] <- 1 # Preventing overflow values after normalization

  # _____________________________________________________________________
  # Preliminary
  # _____________________________________________________________________
  
  # Dimensions of the dataset
  numPoints <- nrow(X)
  numVars <- ncol(X)

  # Re-defining neighborhood radii
  beta  <- -4/(sqshFactor * Ra)^2 # rb is usually defined as rb = sqshFactor*ra
  

  # _____________________________________________________________________
  # Finding the data point with the highest potential value
  # _____________________________________________________________________
  
  # This point will be the first cluster!
  # The highest potential value will be used as reference for 
  # accepting/rejecting other data points as cluster Centers.
  ref_max_PointPotential <- max(PointPotential)
  idx_max_PointPotential <- which.max(PointPotential)
  
  # _____________________________________________________________________
  # Preparing output variables
  # _____________________________________________________________________
  
  Centers <- matrix(,nrow=0,ncol=numVars) # Cluster Centers will be stored here
  numClusters <- 0 # number of detected potential centers
  CentersPotential <- c() # densities of the potential centers will be stored here
  CentersIdx <- c() # identity of the potential centers will be stored here

  # _____________________________________________________________________
  # Iterative procedure to find cluster centers and subtracting potential
  # from data points around found cluster centers
  # _____________________________________________________________________
  
  current_max_PointPotential <- ref_max_PointPotential
  
  # flag that controls the iterative procedure
  find_more_centers <- 1 
  
  # here begins the iterative procedure
  while (find_more_centers & current_max_PointPotential) {
    # If a data point has a potential that leads to a ratio
    # smaller than "acceptRatio" and "rejectRatio", then this flag will
    # not change its value and the iterative procedure will stop.
    find_more_centers <- 0
    
    # finding datapoint with the highest potential
    X_max_PointPotential <- X[idx_max_PointPotential,]
    
    # Ratio between the first highest potential and the current one.
    # We will accept a data point as a potential cluster when this ratio
    # is largest than the "acceptRatio".
    max_PointPotential_ratio <- current_max_PointPotential/ref_max_PointPotential
  
    if (max_PointPotential_ratio > acceptRatio) {
      # the ratio is high enough, this data point will be accepted as a
      # cluster center!
      find_more_centers <- 1
    }
    else if (max_PointPotential_ratio > rejectRatio) {
      # Accept the data point only if it has a sufficient potential and is
      # sufficiently far from the current set of cluster centers!

      # Finding the shortest distance between the datapoint with the highest
      # potential and all the previously found cluster centers
      d <- apply(Centers,1,function(x)sqrt(sum((x - X_max_PointPotential)^2)))
      dmin <- min(d)
      
      # Evaluation of the potential cluster center
      if ((dmin/Ra + max_PointPotential_ratio) >= 1) {
        find_more_centers <- 1	# This data point could be accepted as a cluster center
      }
      else {
        find_more_centers <- 2	# Remove this point from further consideration, and continue
      }
    }
    
    # Proceeding according the "find_more_centers" flag
    if (find_more_centers == 1) {
      # adding the data point to the list of cluster centers
      Centers <- unname(rbind(Centers,X_max_PointPotential))
      # updating the number of clusters
      numClusters <- numClusters + 1
      # storing the center identity (which data point is) and its density measure
      CentersIdx <- c(CentersIdx,idx_max_PointPotential)
      CentersPotential <- c(CentersPotential,current_max_PointPotential)
      
      # Subtracting potential from data points near the current cluster center
		  tmp <- matrix(X_max_PointPotential,nrow=numPoints,ncol=numVars,byrow=TRUE)
      PointPotential <- PointPotential - current_max_PointPotential*exp(beta*rowSums((tmp - X)^2))
		  
      # Annulling negative potentials
		  PointPotential[PointPotential < 0] <- 0
		  
		  # Finding the data point with the highest potential value.
		  # This point will be the next potential cluster!
		  current_max_PointPotential <- max(PointPotential)
		  idx_max_PointPotential <- which.max(PointPotential)
		}
		else if (find_more_centers == 2) {
		  # Removing data point from further consideration by setting its potential
		  # to zero
		  PointPotential[idx_max_PointPotential] <- 0

		  # Finding the data point with the highest potential value.
		  # This point will be the next potential cluster!
		  current_max_PointPotential <- max(PointPotential)
		  idx_max_PointPotential <- which.max(PointPotential)
		}
  }
  # here ends the iterative procedure
  
  # _____________________________________________________________________
  # Scaling cluster Centers from the unit hiperbox to the original
  # data space
  #_____________________________________________________________________
  
  for (i in 1:numVars) {
    Centers[,i] = (Centers[,i] * (maxX[i] - minX[i])) + minX[i];
  }
  
  # _____________________________________________________________________
  # Building ouput object
  # _____________________________________________________________________
  subcluster <- list(centers=Centers,centers.number=numClusters,
                     centers.identity=CentersIdx,centers.densities=CentersPotential,
                     input.parameters=list(ra=Ra,squash.factor=sqshFactor,
                                           accept.ratio=acceptRatio,reject.ratio=rejectRatio))

  return(subcluster)
}
