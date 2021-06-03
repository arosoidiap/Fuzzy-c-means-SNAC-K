cMeansMembershipSummary <- function(mca,num.clusters,repetitions,m)
{
  library(snow)
  library(foreach)
  library(doSNOW)
  library(e1071)
  
  #Empty data frame to store clusterization per repetition 
  #membership.current <- data.frame(array(0,dim=c(nrow(mca$ind$coord),num.clusters,repetitions))
  membership.current <- matrix(0,nrow=nrow(mca),ncol=num.clusters)
  set.seed(1234)
  
  #Setting of the progress bar
  pb <- txtProgressBar(max = repetitions, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  membership <- list()
  
  
  #Parallelized computation of the repetitions
  membership <- foreach(j=1:repetitions,combine=rbind,.options.snow=opts) %dopar%{

    #Computing the indexes per number of clusters
    library(e1071)
    #Proper place to set the seed and obtain different results per each iteration
    
    cmean <- cmeans(mca,centers=num.clusters,iter.max=1000,m=m)
    
    #changing the cluster numbers to have it ordered by increasing number of elements
    #All the clusters from the several realizations will have the same cluster order and name
    ord <- order(cmean$size)
    for(i in 1:length(cmean$size)){
      membership.current[,i] <- cmean$membership[,ord[i]]
    }
    #Upgrading the state of the progress bar
    setTxtProgressBar(pb, j)
    
    return(membership.current=membership.current)
  }
  
  #Compute mean of memberships across all realizations
  membership.current.ave <- array(0,dim=c(nrow(mca),num.clusters,repetitions))
  for(i in 1:repetitions){
    membership.current.ave[,,i] <- membership[[i]]
  }
  
  
  #Mean of memberships matrices across realizations
  membership$membership <- apply(membership.current.ave,1:2,mean)
  
  #Computing cluster with greater membership per individual
  membership$cluster <- apply(membership$membership,1,which.max)
  
  #Obtaining ids
  #membership$id <- id[as.numeric(row.names.data.frame(mca)),]
  
  return(list(membership=membership))
  
}