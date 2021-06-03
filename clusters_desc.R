load("X:\\Analysis internship\\GitHub/Indexes_data_BL.Rdata")

library(plyr)
library(dtplyr)
library(xlsx)
library(data.table)

library(doSNOW)
cl<-makeCluster(16)
registerDoSNOW(cl)

num.clusters <- 6
repetitions <- 100
m <- 1.1
mca <- pca

source('cMeansMembershipSummary.R')
# source('C:\\Users\\Albert.Roso\\Dropbox\\Albert_Roso_internship\\cMeansMembershipSummary.R')

### Memberships c-means
set.seed(123)
cmeans_prueba_6c <- cMeansMembershipSummary(mca,num.clusters,repetitions,m)
cmeans_prueba_6c_memb <- cmeans_prueba_6c$membership$membership
rm(cmeans_prueba_6c)
table(apply(cmeans_prueba_6c_memb,1,sum))
stopCluster(cl)

###Diseases
##original data
load("X:\\Analysis internship\\Paper I/Results/data_SNACK_new_rev.RData")

#Exclude no MM
pop_MM<-pop[pop$Chron_num>=2,]
rm(pop)

#Exclude disease less than <2%
BL<-which(strtrim(names(pop_MM),2)=="BL")
fp_BL<-apply(pop_MM[,BL],2,sum)
n_BL<-length(which(!is.na(pop_MM[,BL[1]])))
which(fp_BL/n_BL>=0.02)
BL_2pct<-names(which(fp_BL/n_BL>=0.02))

prob.paci.prevalent.num <- as.data.table(pop_MM[,BL_2pct])

#Storing as a matrix and transposing
prob.paci.prevalent.matrix <- t(as.matrix(prob.paci.prevalent.num)) 

#multiplying per membership 
memb_pond <- prob.paci.prevalent.matrix %*% cmeans_prueba_6c_memb

#sum of memberships per cluster
sum_memb_fuzzy <- apply(cmeans_prueba_6c_memb,2,sum)

#prevalence with memberships
prev_memberships <- memb_pond

for(i in 1:num.clusters){
  prev_memberships[,i] <- 100*memb_pond[,i] / sum_memb_fuzzy[i]
}

prev_memberships.df <- as.data.frame(prev_memberships)

for(j in 1:num.clusters){
  prevalence <- paste0("Prev c",j)
  names(prev_memberships.df)[j] <- prevalence
}

#Mean prevalence to compute difference observed vs expected
prev_memberships.num <- apply(prev_memberships.df,2,as.numeric)

prevalence.mean <- prev_memberships.num %*% sum_memb_fuzzy/sum(sum_memb_fuzzy) 

prev_memberships.df <- cbind(prev_memberships.num,prev.mean=prevalence.mean)

colnames(prev_memberships.df)[num.clusters+1] <- "prev.mean"

#Chronic categhories
prev_memberships.df <- cbind(prev_memberships.df,Problem=rownames(prev_memberships))


for(i in 1:num.clusters){
  oe <- paste0("OE c",i)
  prev_memberships.df <- cbind(prev_memberships.df,as.numeric(as.character(prev_memberships.df[,i]))/as.numeric(as.character(prev_memberships.df[,num.clusters+1])))
  colnames(prev_memberships.df)[i+2+num.clusters] <- oe
}

#Computing exclusivity
sum_malalt <- apply(memb_pond,1,sum)

exclusivity <- memb_pond

for(i in 1:nrow(exclusivity)){
  exclusivity[i,] <- 100*memb_pond[i,] / sum_malalt[i]
}

#Adding exclusivity to Data sheet

prev_memberships.df <- cbind(prev_memberships.df,exclusivity)

for(i in 1:num.clusters){
  exc <- paste0("Exc c",i)
  colnames(prev_memberships.df)[i+2+2*num.clusters] <- exc
}

for(i in 1:(num.clusters+1)){
  prev_memberships.df[,i] <- as.numeric(as.character(prev_memberships.df[,i]))
}

prev_memberships.fuzzy <- prev_memberships.df

write.xlsx(prev_memberships.fuzzy, "ProbsXClust_6c_100realiz_BL.xlsx",sheetName = "Clusters")

###Table with baseline variables

##numeric

data.comp.num <- pop_MM[,c("age","Chron_num","mmse","ADL_1","IADL1_ADL","dis_scoreADL1","balance","grip_strength",
                           "walking_speed","number_drugs","S_Albumin_grL","S_Creatinine_umolL","CRP_mmolL",
                           "prop_LSI","prop_worry","prop_orient","prop_outlookI","prop_resisI","prop_roleI",
                           "hoursmonthIADL.x", "hoursmonthADL.x", "hourmonthbaseform",
                           "hoursmonthIADL.y", "hoursmonthADL.y", "hourmonthbaseinform",
                           "n_spec0", "n_spec1", "n_spec2",
                           "tot_unplanned0", "tot_unplanned1", "tot_unplanned2",
                           "t_up_v0", "t_up_v1", "t_up_v2")]
library(DataExplorer)

plot_missing(data.comp.num)

summary.cont<-matrix(nrow=ncol(data.comp.num),ncol=num.clusters+2)
row.names(summary.cont)<-colnames(data.comp.num)

for(j in 1:ncol(data.comp.num))
{
data.comp.num.aux<-cbind(data.comp.num[,j],cmeans_prueba_6c_memb)
data.comp.num.aux<-na.omit(data.comp.num.aux)

data.comp.num.var <- t(as.matrix(data.comp.num.aux[,1])) 
data.comp.num.memb <-as.matrix(data.comp.num.aux[,-1])  

memb_pond_num <- data.comp.num.var %*% data.comp.num.memb
sum_memb_fuzzy_num<-apply(data.comp.num.memb,2,sum)

for(i in 1:num.clusters){
  memb_pond_num[,i] <- memb_pond_num[,i] / sum_memb_fuzzy_num[i]
}

summary.cont[j,1:num.clusters]<-memb_pond_num

summary.cont.mean <- summary.cont[j,1:num.clusters] %*% sum_memb_fuzzy_num/sum(sum_memb_fuzzy_num)

print(summary.cont.mean);print(mean(data.comp.num[,j],na.rm=TRUE))

summary.cont[j,(num.clusters+1):(num.clusters+2)]<-c(summary.cont.mean,sum(sum_memb_fuzzy_num))

}

summary.cont

colnames(summary.cont)<-c(paste("Cluster",1:num.clusters),"All","n")

summary.cont<-rbind(summary.cont,c(sum_memb_fuzzy,mean(sum_memb_fuzzy),sum(sum_memb_fuzzy)))

summary.cont<-rbind(summary.cont,summary.cont[nrow(summary.cont),]/sum(sum_memb_fuzzy))

rownames(summary.cont)[(nrow(summary.cont)-1):nrow(summary.cont)]<-c("No individuals","No individuals (%)")

write.xlsx(summary.cont, "ProbsXClust_6c_100realiz_BL.xlsx",sheetName = "Variables_cont", append = T)

##categorical
data.comp.cat <- pop_MM[,c("sex","agegr12","selfrated_health","chair_stand","education","occupation",
                           "economy","civilstatus_num","livingplace","smoking","alcohol","bmi",
                           "physical_activity","socialnetwork_tertiles","TExtrav_cat","TNeurot_cat",
                           "TOpenn_cat","death")]

data.comp.cat[] <-lapply(data.comp.cat,factor)

plot_missing(data.comp.cat)

summary.cat<-NULL

library(dummies)

for(j in 1:ncol(data.comp.cat))
{
  data.comp.cat.aux<-cbind(as.numeric(data.comp.cat[,j]),cmeans_prueba_6c_memb)
  data.comp.cat.aux<-na.omit(data.comp.cat.aux)
  
  data.comp.cat.var <- t(as.matrix(dummy(data.comp.cat.aux[,1]))) 
  data.comp.cat.memb <-as.matrix(data.comp.cat.aux[,-1])  
  
  memb_pond_cat <- data.comp.cat.var %*% data.comp.cat.memb
  sum_memb_fuzzy_cat<-apply(data.comp.cat.memb,2,sum)
  
  for(i in 1:num.clusters){
    memb_pond_cat[,i] <- memb_pond_cat[,i] / sum_memb_fuzzy_cat[i]
  }
  
  rownames(memb_pond_cat)<-levels(data.comp.cat[,j])
  
  summary.cat.mean <- memb_pond_cat %*% sum_memb_fuzzy_cat/sum(sum_memb_fuzzy_cat)
  
  summary.cat.aux<-cbind(memb_pond_cat,summary.cat.mean,sum(sum_memb_fuzzy_cat))
  
  summary.cat<-rbind(summary.cat,summary.cat.aux)
  
}

summary.cat

rownames(summary.cat)[1:2]<-c("Male","Female")
rownames(summary.cat)[(nrow(summary.cat)-1):nrow(summary.cat)]<-c("Alive","Dead")

colnames(summary.cat)<-c(paste("Cluster",1:num.clusters),"All","n")

summary.cat<-rbind(summary.cat,c(sum_memb_fuzzy,mean(sum_memb_fuzzy),sum(sum_memb_fuzzy)))

summary.cat<-rbind(summary.cat,summary.cat[nrow(summary.cat),]/sum(sum_memb_fuzzy))

rownames(summary.cat)[(nrow(summary.cat)-1):nrow(summary.cat)]<-c("No individuals","No individuals (%)")

write.xlsx(summary.cat, "ProbsXClust_6c_100realiz_BL.xlsx",sheetName = "Variables_cat", append = T)

###K-means proxy

summary(apply(cmeans_prueba_6c_memb,1,max))
apply(cmeans_prueba_6c_memb,1,which.max)
table(apply(cmeans_prueba_6c_memb,1,which.max))
 
library(compareGroups)
 
data.comp.num$cluster<-apply(cmeans_prueba_6c_memb,1,which.max)
dim(data.comp.num)

res1<-compareGroups(cluster~.,data=data.comp.num,max.ylev = num.clusters)
res1
createTable(res1,show.all = TRUE,show.n = TRUE,show.p.mul = TRUE)
export2xls(createTable(res1,show.all = TRUE,show.n = TRUE,show.p.mul = TRUE),file = "cont_var_6c_BL.xls")

data.comp.cat$cluster<-apply(cmeans_prueba_6c_memb,1,which.max)
dim(data.comp.cat)

res2<-compareGroups(cluster~.-agegr12,data=data.comp.cat,max.ylev = num.clusters,max.xlev=12,chisq.test.perm=T)
res2
createTable(res2,show.all = TRUE,show.n = TRUE,show.p.mul = TRUE)
export2xls(createTable(res2,show.all = TRUE,show.n = TRUE,show.p.mul = TRUE),file = "cat_var_6c_BL.xls")
