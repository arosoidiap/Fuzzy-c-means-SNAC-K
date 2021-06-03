indexes <- function(pca,num.clusters,repetitions,m)
{
  library(e1071)
  library(snow)
  library(foreach)
  library(doSNOW)
  library(fpc)
  
  #Empty data frames to store indexes per repetition
  fukuyama <- data.frame(matrix(0,nrow = length(num.clusters),ncol=repetitions))
  xie <- data.frame(matrix(0,nrow = length(num.clusters),ncol=repetitions))
  part.coeff <- data.frame(matrix(0,nrow = length(num.clusters),ncol=repetitions))
  part.ent <- data.frame(matrix(0,nrow = length(num.clusters),ncol=repetitions))
  calinski <- data.frame(matrix(0,nrow = length(num.clusters),ncol=repetitions))
  
  #Proper place to set the seed and obtain different results per each iteration
  set.seed(1234)
  
  #Setting of the progress bar
  pb <- txtProgressBar(max = repetitions, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  index <- list()
  
  #Parallelized computation of the repetitions
  index <- foreach(j=1:repetitions,combine=rbind,.options.snow=opts) %dopar%{
    
    library(e1071)
    library(fpc)
    
    #Upgrading the state of the progress bar
    setTxtProgressBar(pb, j)
    
    #Computing the indexes per number of clusters
    for(i in 1:length(num.clusters)){ 
      cmeans_prueba <- cmeans(pca,num.clusters[i],m=m)
      fukuyama[i,j] <- fclustIndex(cmeans_prueba,pca,index="fukuyama.sugeno")
      xie[i,j] <- fclustIndex(cmeans_prueba,pca,index="xie.beni")
      part.coeff[i,j] <- fclustIndex(cmeans_prueba,pca,index="partition.coefficient")
      part.ent[i,j] <- fclustIndex(cmeans_prueba,pca,index="partition.entropy")
      calinski[i,j] <- calinhara(pca,cmeans_prueba$cluster,num.clusters[i])
      print(paste0("Processing number of clusters", num.clusters[i], " realization ", j))
    }
    return(list(fukuyama=fukuyama[,j],xie=xie[,j],part.coeff=part.coeff[,j],part.ent=part.ent[,j],calinski=calinski[,j]))
  }
  
  
  #Empty data frames to store average indexes
  fukuyama.ave <- vector(mode="numeric", length=length(num.clusters))
  xie.ave <- vector(mode="numeric", length=length(num.clusters))
  part.coeff.ave <- vector(mode="numeric", length=length(num.clusters))
  part.ent.ave <- vector(mode="numeric", length=length(num.clusters))
  calinski.ave <- vector(mode="numeric", length=length(num.clusters))
  
  #Adding indexes from all the repetitions  
  for(k in 1:repetitions){
    fukuyama.ave <- fukuyama.ave + index[[k]]$fukuyama
    xie.ave <- xie.ave + index[[k]]$xie
    part.coeff.ave <- part.coeff.ave + index[[k]]$part.coeff
    part.ent.ave <- part.ent.ave + index[[k]]$part.ent
    calinski.ave <- calinski.ave + index[[k]]$calinski
  }
  
  #Dividing by number of repetitions to obtain the mean indexes per cluster
  fukuyama.ave <- fukuyama.ave / repetitions
  xie.ave <- xie.ave / repetitions
  part.coeff.ave <- part.coeff.ave / repetitions
  part.ent.ave <- part.ent.ave / repetitions
  calinski.ave <- calinski.ave / repetitions
  
  return(list(index,fukuyama.ave=fukuyama.ave, xie.ave=xie.ave, 
              part.coeff.ave=part.coeff.ave, part.ent.ave=part.ent.ave, calinski.ave=calinski.ave ))
}