#' Apply SIB ComBat
#' 
#' This function applies Spatially-Informed Iterative Block (SIB) ComBat to harmonize a given dataset.
#' 
#' @param Y (p x n) data matrix for which the p rows are features, and the n columns are participants
#' @param Batch (n x 1) numeric or character vector of length n indicating the site/scanner/study id
#' @param X (n x b) matrix, with n as the number of individuals, and b as the number of covariates
#' @param nclust number of clusters to divide features in if no assignments are inputted
#' @param Clust_list list of cluster assignments for features for each SIB ComBat iteration (if NULL, assignments are determined randomly)
#' @param niter number of SIB ComBat iterations to perform
#' @param ncores number of cores to parallelize SIB ComBat harmonization across
#' @param res option to determine if covariates should be included or not (default is FALSE)
#' @param pres option to determine if jobs should be assigned dynamically when parallelizing (default is FALSE)
#' @param cb option to determine if ComBat should be applied first (default is FALSE)
#' 
#' @return Returns final harmonized data after performing SIB ComBat
#' @import parallel
#' @importFrom neuroCombat neuroCombat
#' 
#' @examples 
#' data(Y)
#' data(X)
#' data(Batch)
#' Yharm4 <- sib_combat(Y,Batch,X,nclust=2,Clust_list=NULL,niter=2,ncores=1,res=FALSE,pres=FALSE,cb=FALSE)
#' 
#' @export
#' 

sib_combat <- function(Y,Batch,X,nclust=NULL,Clust_list=NULL,niter,ncores,res=FALSE,pres=FALSE,cb=FALSE){
  
  #Make Batch into a factor variable if it is not already
  if(is.factor(Batch)==FALSE){
    Batch <- as.factor(Batch)
  }
  bids <- unique(Batch)
  nbatch <- length(bids)
  
  #Make X into a matrix if it is not
  if(is.matrix(X)==FALSE){
    X <- as.matrix(X)
  }
  
  #If specified, harmonize with ComBat first
  if(isTRUE(cb)){
    Youtput1 <- neuroCombat(dat=Y,batch=Batch,mod=X,verbose=FALSE)
    Y <- Youtput1$dat.combat
  }
  
  #1st step: standardize the data and save estimates of sigma_v,alpha,beta
  Z_ests <- stand_data(Y,Batch,X)
  Z <- Z_ests[[1]]
  sigma_v <- Z_ests[[2]]
  standmean <- Z_ests[[3]]
  modmean <- Z_ests[[4]]
  
  #If you do not provide a list of cluster assignments, determine them randomly for each iteration
  if(is.null(Clust_list)){
    #For first iteration, split data evenly along nclust clusters
    nfea <- nrow(Z)
    a <- seq(1,nfea)
    clust_assign <- split(a, cut(seq_along(a), breaks = nclust, labels = FALSE))
    no_clust <- as.numeric(unlist(lapply(clust_assign,length)))
    Clust <- rep(NA,nfea)
    start <- 1
    for(k in 1:nclust){
      clust_length <- no_clust[k]
      end <- start + clust_length - 1
      Clust[start:end] <- k
      start <- end + 1
    }
    
    for(j in 1:niter){
      #Get shifted data for iteration 
      Z <- spa_combat1_diag2(Z,Batch,Clust,bids,nbatch,ncores,pres=pres)
      #Create new cluster assignment for next iteration
      Clust <- sample(Clust)
    }
    
  }else{
    for(j in 1:niter){
      #Extract cluster assignment from list
      Clust <- Clust_list[[j]]
      #Get shifted data for iteration 
      Z <- spa_combat1_diag2(Z,Batch,Clust,bids,nbatch,ncores,pres=pres)
    }
  }
  
  if(isTRUE(res)){
    Yadj <- sigma_v * Z
  }else{
    Yadj <- (sigma_v * Z) + modmean + standmean
  }
  
  #Output harmonized data
  return(Yadj)
}
