# ++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# Author: Leandro Ariza-Jimenez (larizaj@eafit.edu.co) #
# ++++++++++++++++++++++++++++++++++++++++++++++++++++ #

Fmeasure <- function(true_labels,estimated_labels){
  
    # table
  A <- table(estimated_labels,true_labels)
  
  # Precision
  den <- sum(A)
  num <- sum(apply(A,1,max)) # sum of maximum values per row
  P <- num/den
  
  # Recall
  num <- sum(apply(A,2,max)) # sum of maximum values per row
  R <- num/den
  
  # F-measure
  Fmeas <- 2*R*P/(P+R)
  
  return(Fmeas)
  
}