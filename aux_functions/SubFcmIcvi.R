# ++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# Author: Leandro Ariza-Jimenez (larizaj@eafit.edu.co) #
# ++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# An unsupervised clustering method
# parallel computation of the proposed method based on the combination of
# Subtractive + "Fuzzy C-means" + "internal clustering validition index"
#
# This script is a shorter version of the script used in:
# -------------------------------------------------------------------------------------
# Ariza-Jiménez, L., Quintero, O. L., & Pinel, N. (2018). Unsupervised fuzzy binning of
# metagenomic sequence fragments on three-dimensional Barnes-Hut t-Stochastic Neighbor
# Embeddings.
# In 2018 40th Annual International Conference of the IEEE Engineering in Medicine and
# Biology Society (EMBC). https://doi.org/10.1109/EMBC.2018.8512529
# -------------------------------------------------------------------------------------


SubFcmIcvi <- function(X,icvi){
  
  
  # parallel stuff
  numCores <- detectCores()
  cl <- makeCluster(numCores, type="FORK")
  registerDoParallel(cl)
  
  
  # values for Ra (a subtractive clustering input parameter)
  ra.vec <- seq(from=0.1, to=0.7, by=0.05)
  
  
  # computing a pair-wise distance matrix for the input embedding
  dist.mat.norm <- as.matrix(parDist(dataNorm(X), method = "euclidean"))
  
  
  ###################################################################
  # computing the potential of being a cluster center for every
  # object in the current embedding as a function of its distance to all
  # the other data points
  # a data point with many neighboring data points (i.e. with a high density
  # of surrounding data points) will have a high potential value.
  # BELOW IS THE FIRST STAGE OF THE SUBTRACTIVE CLUSTERING ALGORITHM!!!
  ###################################################################
  
  
  # parallel computation of the potential of every object in the current embedding for
  # different values of ra
  pointPotential <- foreach(X=isplitCols(dist.mat.norm, chunks=numCores), .combine='+') %dopar% 
  {
    foreach(ra=ra.vec, .combine='cbind') %do% {rowSums(exp(-4/ra^2*(X)^2))}
  }
  
  
  # unsupervised clustering
  # parallel computation of the proposed method based on the combination of
  # Subtractive + "Fuzzy C-means" + "internal clustering validition index"
  SubFcmIcvi <- foreach(j=1:length(ra.vec)) %dopar%
  {
    # running the subtractive clustering algorithm based on the embedding data,
    # the already-computed potentials of every embedding data point, and a
    # given ra value
    subclust <- subclustCore(X=X,
                             PointPotential=pointPotential[,j],
                             Ra=ra.vec[j])
    
    # Some errors arises when the clustering has only one group
    # and when FCM is computed based on this result
    if (nrow(subclust$centers) > 1) {
      # fcm should run, in the future, on the normalized data space!!!
      fcmclust <- cmeans(X,
                         centers=subclust$centers, iter.max=25, verbose=FALSE,
                         dist = "euclidean", method = "cmeans", m = 2)
      
      # computing cluster validation index
      icvi.val <- intCriteria(X, fcmclust$cluster, crit=icvi)
      
    }  else {
      fcmclust <- NULL # nothing to do
      icvi.val <- NaN # NaN are used for tracking purpouses
    }
    
    # output
    list(subclust=subclust, fcmclust=fcmclust, icvi.val=icvi.val)
  }
  
  # retrieving the Silhouette index for each Ra instance
  icvi.vec <- lapply(1:length(ra.vec), function(x) SubFcmIcvi[[x]]$icvi.val)
  icvi.vec <- unname(unlist(icvi.vec))
  
  # fixing NaNs introduced when the number of centers is equal to 1,
  # (replacing NaNs with the worst case value of Silhouette index)
  icvi.vec[is.nan(icvi.vec)] <- 0
  
  # determining the optimal clustering accoding to Silhouette index
  bestRa.pos <- bestCriterion(icvi.vec, crit=icvi)
  
  # optimal clustering according to Silhouette index
  SubFcmIcvi <- SubFcmIcvi[[bestRa.pos]]
  
  # stop cluster
  stopCluster(cl)
  
  # return
  return(SubFcmIcvi)
}