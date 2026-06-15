#Helper Functions for Covariance-Aware Multivariate (CAM) ComBat

#Function to calculate moore_penrose_inverse
moore_penrose_inverse <- function(A){
  #If a scalar, output scalar value as matrix
  if(length(A)==1){
    A_pinv <- as.matrix(1/A)
  }else{
    AtA <- t(A) %*% A
    eig <- eigen(AtA)
    
    # Eigen decomposition
    values <- eig$values
    vectors <- eig$vectors
    
    # Invert only eigenvalues that explain at least 99.99% of variability; otherwise set to 0
    cum_eigs <- cumsum(values)/sum(values)
    tol_ind <- which(cum_eigs > 0.999)[1]
    tol <- values[tol_ind]
    values_inv <- ifelse(values >= tol, 1 / values, 0)
    
    # Construct pseudoinverse of A^T A
    AtA_pinv <- vectors %*% diag(values_inv) %*% t(vectors)
    
    # Generalized inverse
    A_pinv <- AtA_pinv %*% t(A)
  }
  return(A_pinv)
}

#Create function to standardize data (used code for neuroCombat from github Jfortin1)
stand_data <- function(Y,Batch,X=NULL){
  
  #Make Batch into a factor variable if it is not already
  if(is.factor(Batch)==FALSE){
    Batch <- as.factor(Batch)
  }
  #Find the total number of batches
  nbatches <- nlevels(Batch)
  
  #Make matrix corresponding to which site each individual belongs to
  batch_mat <- model.matrix(~-1 + Batch)
  
  #Find the number of individuals from each site/batch and total number of individuals
  ni_array <- tabulate(Batch)
  N <- ncol(Y)
  
  #Create design matrix (without intercept)
  if(is.null(X)){
    Xdes <- batch_mat
    mod_mean <- matrix(0,nrow(Y),N)
  } else{
    #Make X into a matrix if it is not already
    if(is.matrix(X)==FALSE){
      X <- as.matrix(X)
    }
    Xdes <- cbind(batch_mat,X)
  }
  
  ##Calculate matrix of OLS estimates for all voxels/features: (XTX)^-1XTY
  XTX <- crossprod(Xdes)
  XTXinv <- solve(XTX)
  XTXinvXT <- tcrossprod(XTXinv,Xdes)
  ols_mat <- tcrossprod(XTXinvXT,Y)
  
  #Calculate grand and standardized mean (based on neuroCombat github)
  grand_mean <- crossprod(ni_array/N,ols_mat[1:nbatches,]) #1xV vector
  stand_mean <- crossprod(grand_mean, t(rep(1,N))) #VxN matrix
  
  #Calculate pooled variance (based on neuroCombat github)
  factors <- (N/(N-1))
  Yres <- Y - t(Xdes %*% ols_mat)
  pooled_var <- apply(Yres,1,var)/factors
  #replacing when pooled variance is 0 to median of calculated pooled variances
  pooled_var[pooled_var==0] <- median(pooled_var[pooled_var!=0],na.rm=TRUE) 
  
  if(!is.null(X)){
    Xdes2 <- Xdes
    Xdes2[,1:nbatches] <- 0
    mod_mean <- t(Xdes2 %*% ols_mat)
  }
  
  Sigma_v <- tcrossprod(sqrt(pooled_var),rep(1,N)) #Vx1 vector of sigma_v estimates
  Z <- (Y - stand_mean - mod_mean)/Sigma_v
  
  #Output standardized data with other variables used in standardization
  return(list(Z,Sigma_v,stand_mean,mod_mean))
}

#Create functions to calculate batch-specific method of moments estimators 
#(given some standardized data Z from batch i, where Z is a ni x V matrix
#where each column corresponds to one feature/voxel)

#MM for batch gamma_i
batch_mm_gami <- function(Z){
  
  #Find average Zij to get method of moments estimator
  gam_i <- rowMeans(Z)
  return(gam_i)
}

#MM for batch Sigma_i
batch_mm_Sigmai <- function(Z){
  
  #Calculate method of moments estimator for gamma_i
  gam_i <- batch_mm_gami(Z)
  
  #Calculate number of unique voxels/features and subject IDs
  nvox <- nrow(Z)
  nids <- ncol(Z)
  
  #Calculate Sigma_i
  Sigma_i <- matrix(0,nrow=nvox,ncol=nvox)
  for(j in 1:nids){
    Zvox <- Z[,j]
    Zdiff2 <- tcrossprod(Zvox - gam_i)
    Sigma_i <- Sigma_i + Zdiff2
  }
  Sigma_i <- Sigma_i/nids #or nids - 1
  
  return(Sigma_i)
}

#MM for batch gamma_0i
batch_mm_gam0i <- function(Z){
  
  #Calculate method of moments estimator for gamma_i
  gam_i <- batch_mm_gami(Z)
  gam_0i_mean <- mean(gam_i)
  
  return(gam_0i_mean)
}

#MM for batch tau_0i^2
batch_mm_tau0i <- function(Z){
  
  #Calculate method of moments estimator for gamma_i
  gam_i <- batch_mm_gami(Z)
  
  #Calculate estimate for tau_0i^2 by calculating variance
  tau_0i <- var(gam_i)
  
  return(tau_0i)
}

#Create function to get EB Batch Parameter effect estimate for gamma_i (batch specific)
#(given Sigma_i calculated previously)
batch_eb_gami <- function(Z,Sigma_i){
  
  #Calculate method of moments estimator for gamma_0i, tau_0i, and Sigma_i
  gam_0i <- batch_mm_gam0i(Z)
  tau_0i <- batch_mm_tau0i(Z)
  
  ni <- ncol(Z) #assuming each column corresponds to a different individual
  nvox <- nrow(Z) #assuming each row corresponds to a different feature
  
  Sigmai_inv <- moore_penrose_inverse(Sigma_i)
  
  gamistar11 <- ni*tau_0i*Sigmai_inv + diag(1,nvox)
  gamistar1 <- moore_penrose_inverse(gamistar11)
  
  sumZij <- rowSums(Z,na.rm=TRUE) #assuming each column corresponds to a different individual
  
  gamistar2 <- gam_0i + tau_0i*(Sigmai_inv %*% sumZij)
  gamistar <- gamistar1 %*% gamistar2
  
  return(gamistar)
}


#Create function to calculate estimate of underlying Covariance matrix, as well as batch specific covariance
#(given standardized data)
cov_est <- function(Z,Batch){
  
  #Assume Batch is already categorized as a factor/categorical variable
  batch_ids <- unique(Batch)
  nvox <- nrow(Z)
  ni_array <- tabulate(Batch)
  nbatch <- length(batch_ids)

  Sigma_list <- rep_len(list(0),nbatch)
  
  #Iterate through batches to calculate Sigma
  for(p in 1:nbatch){
    bid <- batch_ids[p]
    bind <- which(Batch==bid)
    ni <- ni_array[p]
    Zb <- Z[,bind]
    
    #Calculate covariance for batch i and add to list of covariances
    Sigma_i <- batch_mm_Sigmai(Zb)
    Sigma_list[[p]] <- Sigma_i
  }
  
  return(Sigma_list)
}

#Create function to calculate estimate of A_i inverse from batch i for adjustment 
#(given batch specific covariance and overarching covariance matrix)
Ainv_est <- function(Sigma_i,Sigma){
  
  #Calculate matrix from Spectral decomposition for Sigmai
  Leigs_i <- eigen(Sigma_i)
  Lvecs_i <- Leigs_i$vectors
  
  #Replace values where eigenvalues are less than 1e-8 to be 0
  Lvals_i <- Leigs_i$values
  cum_eigs_i <- cumsum(Lvals_i)/sum(Lvals_i)
  zero_ind_i <- which(cum_eigs_i > 0.99) #0.95
  Lvals_i[zero_ind_i] <- 0
  Lvals_sr_i <- sqrt(Lvals_i)
  
  #Take square root of eigenvalues of Sigmai
  diag_eigs_i <- diag(Lvals_sr_i)
  
  #Calculate matrix from Spectral decomposition for Sigma
  Leigs <- eigen(Sigma)
  Lvecs <- Leigs$vectors
  
  #Replace values where eigenvalues are less than 1e-8 to be 0
  Lvals <- Leigs$values
  cum_eigs <- cumsum(Lvals)/sum(Lvals)
  zero_ind <- which(cum_eigs > 0.99) #0.95
  Lvals[zero_ind] <- 0
  Lvals_sr <- sqrt(Lvals)
  
  #Take square root of eigenvalues of Sigma
  diag_eigs <- diag(Lvals_sr)
  
  Li <- Lvecs_i %*% diag_eigs_i
  P <- Lvecs %*% diag_eigs
  Pinv <- moore_penrose_inverse(P)
  A_i <- Li %*% Pinv
  Ai_inv <- moore_penrose_inverse(A_i)
  
  return(Ai_inv)
}

#Adjust data from batch i, given some standardized data from batch i
#(assume you previously calculate and save alpha, beta, sigma, Ainv, gamma_i*)
adjust_batchi <- function(Z,modmean,standmean,sigma_v,Ainv,gamistar,res=FALSE){
  
  #modmean: (calculated in stand_data)
  #standmean: (calculated in stand_data)
  #sigma_v: (calculated in stand_data)
  #gamistar: nvox x 1 vector
  
  nids <- ncol(Z)
  gamistar_full <- tcrossprod(gamistar,rep(1,nids))
  Yshift1 <- Ainv %*% (Z - gamistar_full)
  
  if(isTRUE(res)){
    Yadj <- sigma_v * Yshift1
  }else{
    Yadj <- (sigma_v * Yshift1) + modmean + standmean
  }
  
  return(Yadj)
}
