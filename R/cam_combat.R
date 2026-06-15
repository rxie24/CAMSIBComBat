#' Apply CAM ComBat
#' 
#' This function applies Covariance-Aware Multivariate (CAM) ComBat to harmonize a given dataset.
#' 
#' @param Y (p x n) data matrix for which the p rows are features, and the n columns are participants
#' @param Batch (n x 1) numeric or character vector of length n indicating the site/scanner/study id
#' @param X (n x b) matrix, with n as the number of individuals, and b as the number of covariates
#' @param res option to determine if covariates should be included or not (default is FALSE)
#' 
#' @return Returns list in which the last element is the final harmonized dataset and the other elements are the batch-specific covariance estimates
#' 
#' @examples 
#' data(Y)
#' data(X)
#' data(Batch)
#' Youtput3 <- cam_combat(Y,Batch,X)
#' 
#' @export
 

#Perform strategy 1 for spatial combat on all of the data
cam_combat <- function(Y,Batch,X,res=FALSE){
  
  #Make Batch into a factor variable if it is not already
  if(is.factor(Batch)==FALSE){
    Batch <- as.factor(Batch)
  }
  bids <- unique(Batch)
  nbatch <- length(bids)
  
  #Make X into a matrix if it is not
  if(is.matrix(X)==FALSE & !is.null(X)){
    X <- as.matrix(X)
  }
  
  #1st step: standardize the data and save estimates of sigma_v,alpha,beta
  Z_ests <- stand_data(Y,Batch,X)
  Z <- Z_ests[[1]]
  sigma_v <- Z_ests[[2]]
  standmean <- Z_ests[[3]]
  modmean <- Z_ests[[4]]
  
  #2nd: calculate overarching and batch specific covariance matrices
  All_Sigmas <- cov_est(Z,Batch)
  Sigma <- cov(t(Z))
  
  #3rd: iterating through batches, calculate gamma_i* and adjust data 
  #Save harmonized data and batch specific covariances into a list
  results_out <- rep_len(list(0),(1+nbatch))
  
  Yadj <- matrix(NA,nrow=nrow(Y),ncol=ncol(Y))
  for(j in 1:nbatch){
    bid <- bids[j]
    bind <- which(Batch==bid)
    
    #Extract batch specific data
    Zb <- Z[,bind]
    modb <- modmean[,bind]
    standb <- standmean[,bind]
    sigma_vb <- sigma_v[,bind]
    
    #Extract batch specific covariance and save as jth element in prespecified list
    Sigma_i <- All_Sigmas[[j]]
    results_out[[j]] <- Sigma_i
    
    #Calculate gamma_i*
    gamistar <- batch_eb_gami(Zb,Sigma_i)
    
    #Adjust the data accordingly and save harmonized data
    Ainv <- Ainv_est(Sigma_i,Sigma)
    Yharm <- adjust_batchi(Zb,modb,standb,sigma_vb,Ainv,gamistar,res)
    Yadj[,bind] <- Yharm
  }
  
  #Save harmonized data as last element in list
  results_out[[(nbatch+1)]] <- Yadj
  
  #Output harmonized data along with estimates of batch specific coariances
  return(results_out)
}

