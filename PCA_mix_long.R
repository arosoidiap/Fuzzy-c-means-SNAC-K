setwd("X:\\Analysis internship\\GitHub/")

load("X:\\Analysis internship\\Paper II\\Results\\pop_long_MM.RData")

head(pop_long_MM)
summary(pop_long_MM)

#PCAMix
library(PCAmixdata)

times<-c("BL","F2","F4")

years<-list(a=which(strtrim(names(pop_long_MM),2)=="BL"),
          b=which(strtrim(names(pop_long_MM),3)=="F2B"),
            c=which(strtrim(names(pop_long_MM),3)=="F4B"))

for(i in 1:length(times))
{
  y<-pop_long_MM[,years[[i]]]
  row.names(y)<-pop_long_MM$lopnr
  # Missing values are replaced by means for quantitative variables and by zeros in the indicator matrix for qualitative variables.
  y<-na.omit(y)
  names(y)
  
  y2<-data.frame(apply(y,2,as.factor))
  
  head(y2)
  
  summary(y2)
  
  quanti<-NULL
  quali<-1:length(years[[i]])
  
  X1<-NULL
  X2<-y2

  gc()
  # If X.quanti is NULL, only qualitative variables are available and standard MCA is performed
  system.time(obj <- PCAmix(X.quanti=X1,X.quali=X2,rename.level = TRUE,graph = FALSE))
  obj
  dim(obj$ind$coord)
  summary(obj)
  plot(obj)
  
  obj$eig
  sum(obj$eig[,"Eigenvalue"])
  barplot(obj$eig[,"Eigenvalue"])
  
  n<-nrow(X2)
  # p<-ncol(X1)+ncol(X2)
  p<-ncol(X2)
  
  
  lambda<-obj$eig[,"Eigenvalue"]
  lambda.p<-obj$eig[,"Cumulative"]
  
  #Kaiser - Guttman
  which(lambda>1)
  ndim1<-length(which(lambda>1))
  
  gc()
  obj1 <- PCAmix(X.quanti=X1,X.quali=X2,ndim = ndim1,rename.level = TRUE,graph = FALSE)
  Z1<-obj1$ind$coord
  dim(Z1)
  
  #Karlis - Saporta - Spinaki (2003)
  k<-1+2*sqrt((p-1)/(n-1))
  which(lambda>k)
  ndimk<-length(which(lambda>k))
  
  gc()
  objk <- PCAmix(X.quanti=X1,X.quali=X2,ndim = ndimk,rename.level = TRUE,graph = FALSE)
  Zk<-objk$ind$coord
  dim(Zk)
  
  #Variance prop.
  which(lambda.p<=50)
  ndimp<-length(which(lambda.p<=50))
  
  gc()
  objp <- PCAmix(X.quanti=X1,X.quali=X2,ndim = ndimp,rename.level = TRUE,graph = FALSE)
  Zp<-objp$ind$coord
  dim(Zp)
  
  ###Input Data
  # data<-obj$ind$coord
  data<-Zk
  dim(data)
  head(data)
  
  library(psych)
  par(mar=c(1,1,1,1))
  summary(data)
  multi.hist(data)
  apply(data,2,mean)
  apply(data,2,var)
  sum(apply(data,2,var))
  
  save(data,file = paste("PCAmix_",times[i],".RData",sep=""))
}