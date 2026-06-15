#Helper Functions for Spatially-Informed Iterative Block (SIB) ComBat

#Adjust data from batch i, given some standardized data from batch i
#(assume you previously calculate and save Ainv, gamma_i*)
adjust_batchi_iter <- function(Z,Ainv,gamistar){
  
  #modmean: (calculated in stand_data)
  #standmean: (calculated in stand_data)
  #sigma_v: (calculated in stand_data)
  #gamistar: nvox x 1 vector
  
  nids <- ncol(Z)
  gamistar_full <- tcrossprod(gamistar,rep(1,nids))
  Yshift1 <- Ainv %*% (Z - gamistar_full)
  
  return(Yshift1)
}

#Perform CAM ComBat on one block of features for a given iteration
spa_combat1_clust2 <- function(Z,Batch,bids,nbatch){
  
  #Calculate overarching and batch specific covariance matrices
  All_Sigmas <- cov_est(Z,Batch)
  Sigma <- cov(t(Z))
  
  #Iterating through batches, adjust data 
  Yadj <- matrix(NA,nrow=nrow(Z),ncol=ncol(Z))
  for(j in 1:nbatch){
    bid <- bids[j]
    bind <- which(Batch==bid)
    
    #Extract batch specific data
    Zb <- Z[,bind]
    
    #Extract batch specific covariance and save in prespecified list
    Sigma_i <- All_Sigmas[[j]]
    
    #Calculate gamma_i*
    gamistar <- batch_eb_gami(Zb,Sigma_i)
    
    #Adjust the data accordingly and save harmonized data
    Ainv <- Ainv_est(Sigma_i,Sigma)
    Yharm <- adjust_batchi_iter(Zb,Ainv,gamistar)
    Yadj[,bind] <- Yharm
  }
  
  ##Output harmonized data
  return(Yadj)
}

#Perform CAM ComBat across all blocks of features for one iteration
spa_combat1_diag2 <- function(Z,Batch,Clust,bids,nbatch,ncores,pres=FALSE){
  
  #Split the respective components into their clusters
  clust_ind <- ncol(Z) #number of individuals
  Z_clusts <- split(Z, Clust)
  Z_clusts <- lapply(Z_clusts,matrix,ncol=clust_ind)
  
  Yadj_clusts <- mclapply(Z_clusts,function(args)
  {spa_combat1_clust2(args,Batch,bids,nbatch)},
  mc.preschedule=pres,mc.cores=ncores)
  
  #Iterate through clusters to concatenate harmonized features together
  nclust <- length(unique(Clust))
  Yadj <- matrix(NA,nrow=nrow(Z),ncol=ncol(Z))
  for(j in 1:nclust){
    Yadj_clustj <- Yadj_clusts[[j]]
    cids <- which(Clust==j)
    Yadj[cids,] <- Yadj_clustj
  }
  
  #Output harmonized data
  return(Yadj)
}
