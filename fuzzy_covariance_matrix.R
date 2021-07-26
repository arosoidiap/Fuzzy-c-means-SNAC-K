fuzzy_covariance_matrix <- function(membership, pca, centroids, m = 1.1){
  
  #Compute the Fuzzy covariance matrix given the features of individuals and centroids  
  
  nstates <- nrow(centroids)
  dim <- ncol(centroids)
  
  individuals <- nrow(pca)
  
  fuzzy_covariance <- list()
  
  cov <- 0*diag(dim)
  
  for(i in 1:nstates){
    
    temp_cov <- cov
    den_cov <- 0
    
    for(j in 1:individuals){

      temp_cov <- temp_cov + (membership[j,i])^m * ((pca[j,] - centroids[i,]) %*% t(pca[j,] - centroids[i,]))
      den_cov <- den_cov + (membership[j,i])^m
    }
    fuzzy_covariance[[i]] <- temp_cov / den_cov
  }
  
  return(fuzzy_covariance)
}