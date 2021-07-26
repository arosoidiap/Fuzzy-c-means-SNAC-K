divergence_matrix_optimum_cluster_number <- function(centroids1, centroids2, fuzzy_cov1, fuzzy_cov2){
  
  #Compute the kullback-Leibler divergence for model obtained by different fuzzy c-means clusterings
  
  library(monomvn)
  
  nstates1 <- nrow(centroids1)
  nstates2 <- nrow(centroids2)
  
  
  divMat <- matrix(rep(0, nstates1*nstates2), nstates1, nstates2)
  
  for(i in 1:nstates1){
    for(j in 1:nstates2){
      divMat[i,j] <- kl.norm(mu1=centroids1[i,],S1=fuzzy_cov1[[i]],mu2=centroids2[j,],S2=fuzzy_cov2[[j]])
    }
  }
  
  return(divMat)
}