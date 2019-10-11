# ++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# Author: Leandro Ariza-Jimenez (larizaj@eafit.edu.co) #
# ++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# libraries
library(foreach) # for parallelization
library(parallel) # for parallelization
library(iterators) # for parallelization
library(itertools) # for parallelization
library(doParallel) # for parallelization
library(parallelDist) # for computing distance in parallel
library(e1071) # for FCM
library(clusterCrit) # to compute internal clustering validity indices (icvis)
#library(Rtsne) # to use the BH-SNE algorithm
library(clusterGeneration) # to generate the toy dataset
library(RColorBrewer) # color palettes
source("./aux_functions/subclustCore.R") # auxiliar function
source("./aux_functions/Fmeasure.R") # auxiliar function
source("./aux_functions/SubFcmIcvi.R") # auxiliar function (this is the main clustering algorithm)


# toy dataset
X <- genRandomClust(numClust=7, # number of clusters
                    sepVal=0.2, # cluster separation
                    numNonNoisy=2, 
                    numNoisy=0, 
                    numOutlier=0, 
                    numReplicate=1, 
                    fileName="test",  
                    clustszind=2, 
                    clustSizeEq=50, 
                    rangeN=c(50,1000), # cluster size range
                    rangeVar=c(10,100), # cluster variance range
                    clustSizes=NULL,
                    outputDatFlag=FALSE, 
                    outputLogFlag=FALSE, 
                    outputEmpirical=FALSE,
                    outputInfo=FALSE)
Y <- apply(X$datList$test_1, 2, function(x)(x-min(x))/(max(x)-min(x))) # data normalization
labels.true <- X$memList$test_1


# unsupervised clustering
# other icvi options: "Silhouette", "Calinski_Harabasz", "Dunn"
clustering.output <- SubFcmIcvi(Y, icvi="Calinski_Harabasz")

# computing Fmeasure
labels.estim <- clustering.output$fcmclust$cluster
F1 <- Fmeasure(labels.true,labels.estim)

# visualizing clustering output
oldpar <- par(no.readonly = TRUE) # back-up par
layout(matrix(1:2,1,2))
par(mai = c(0.1, 0.1, 0.5, 0.1))
col.palette <- brewer.pal(n = 12, name = "Paired")
plot(Y, col=labels.true, asp=1, pch=16, axes=FALSE, frame.plot=T,
     xlab="", ylab="", main="True labels")
plot(Y, col=col.palette[labels.estim], asp=1, pch=16, axes=FALSE, frame.plot=T,
     xlab="", ylab="", main="Estimate labels")
par(oldpar)