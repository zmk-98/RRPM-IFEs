# Load estimation functions for RRPM with IFEs
source("RPC.R")
# -> RPC(X, y, factor1, factor2, index1)

#' @param X the predictors, a matrix of size NT * p 
#' @param y the response, a vector of length NT
#' @param index1 section of each observation
#' @param K.max the maximum possible number of common components, X * \theta_k
#' @param R.max the maximum possible number of latent factors, F
#' @param cc.K tuning mutiplicative constant in the penalty for K
#' @param cc.R tuning mutiplicative constant in the penalty for R
#' 
#' @return A list containing: selected (K,R), estimated (Theta, B, F, L) for different criteria
choose.KR = function(X, y, v,
                     index1 = c(),
                     K.max = K.max,
                     R.max = R.max,
                     cc.K = 1, cc.R = 2)
{ 
  XX = X
  YY = y
  
  p = ncol(X)
  n = length(y)
  
  ind = unique(index1)
  N = length(ind)
  Ti = n/N
  NTi = (N + Ti) / (N * Ti)
  
  if  (N >0)
  {
    I = list()
    for (i in 1:N)
    {
      I[[i]] = which(index1==ind[i])
    }
  }
  
  fit = list()
  for (K in 0:K.max)
  {
    for (R in 0:R.max)
    { 
      # fit[[K*(R.max+1)+R+1]] = RRPM1(X=X, y=y, factor1=K, factor2=R, index1=index1)
      fit[[K*(R.max+1)+R+1]] = RPC(X=X, y=y, v=v, factor1=K, factor2=R, index1=index1)
      cat("(K=", K, ", R=", R, ") done", "\n", sep="")
    }
  }
  
  sigma2 = fit[[(K.max+1)*(R.max+1)]]$RSS / (N*Ti-K.max*(p+N-K.max)-(Ti+N-R.max)*R.max)
  # sigma2 = fit[[(K.max+1)*(R.max+1)]]$RSS/n
  
  RSS = rep(0, (K.max+1)*(R.max+1))
  MSS = rep(0, (K.max+1)*(R.max+1))
  
  PCp1 = rep(0, (K.max+1)*(R.max+1))
  ICp1 = rep(0, (K.max+1)*(R.max+1))
  PCp2 = rep(0, (K.max+1)*(R.max+1))
  ICp2 = rep(0, (K.max+1)*(R.max+1))
  PCp3 = rep(0, (K.max+1)*(R.max+1))
  ICp3 = rep(0, (K.max+1)*(R.max+1))
  

  
  for (K in 0:K.max)
  {
    for (R in 0:R.max)
    {
      RSS[K*(R.max+1)+R+1] = fit[[K*(R.max+1)+R+1]]$RSS
      
      if (N + Ti < 200) {
        MSS[K*(R.max+1)+R+1] = RSS[K*(R.max+1)+R+1]/(N*Ti-K*(p+N-K)-R*(N+Ti-R))
      } else {
        MSS[K*(R.max+1)+R+1] = RSS[K*(R.max+1)+R+1]/(N*Ti)
      }
      
      PCp1[K*(R.max+1)+R+1] = MSS[K*(R.max+1)+R+1] + sigma2*(cc.K*K+cc.R*R) * NTi*log(1/NTi)
      ICp1[K*(R.max+1)+R+1] = log(MSS[K*(R.max+1)+R+1]) + (cc.K*K+cc.R*R) * NTi*log(1/NTi)
      PCp2[K*(R.max+1)+R+1] = MSS[K*(R.max+1)+R+1] + sigma2*(cc.K*K+cc.R*R) * NTi*log(min(N,Ti))
      ICp2[K*(R.max+1)+R+1] = log(MSS[K*(R.max+1)+R+1]) + (cc.K*K+cc.R*R) * NTi*log(min(N,Ti))
      PCp3[K*(R.max+1)+R+1] = MSS[K*(R.max+1)+R+1] + sigma2*(cc.K*K+cc.R*R) * log(min(N,Ti))/min(N,Ti)
      ICp3[K*(R.max+1)+R+1] = log(MSS[K*(R.max+1)+R+1]) + (cc.K*K+cc.R*R) * log(min(N,Ti))/min(N,Ti)
      # penalty of the third criterion is usually smaller than the first and second criteria
      
    }
  }
  
  
  KR.PCp1 = which.min(PCp1)
  K.PCp1 = (KR.PCp1-1) %/% (R.max+1); R.PCp1 = KR.PCp1 - K.PCp1*(R.max+1) - 1
  KR.ICp1 = which.min(ICp1)
  K.ICp1 = (KR.ICp1-1) %/% (R.max+1); R.ICp1 = KR.ICp1 - K.ICp1*(R.max+1) - 1
  
  KR.PCp2 = which.min(PCp2)
  K.PCp2 = (KR.PCp2-1) %/% (R.max+1); R.PCp2 = KR.PCp2 - K.PCp2*(R.max+1) - 1
  KR.ICp2 = which.min(ICp2)
  K.ICp2 = (KR.ICp2-1) %/% (R.max+1); R.ICp2 = KR.ICp2 - K.ICp2*(R.max+1) - 1
  
  KR.PCp3 = which.min(PCp3)
  K.PCp3 = (KR.PCp3-1) %/% (R.max+1); R.PCp3 = KR.PCp3 - K.PCp3*(R.max+1) - 1
  KR.ICp3 = which.min(ICp3)
  K.ICp3 = (KR.ICp3-1) %/% (R.max+1); R.ICp3 = KR.ICp3 - K.ICp3*(R.max+1) - 1
  
  MSS.matrix = matrix(MSS, K.max+1, R.max+1, byrow = T)
  PCp1.matrix = matrix(PCp1, K.max+1, R.max+1, byrow = T)
  ICp1.matrix = matrix(ICp1, K.max+1, R.max+1, byrow = T)
  PCp2.matrix = matrix(PCp2, K.max+1, R.max+1, byrow = T)
  ICp2.matrix = matrix(ICp2, K.max+1, R.max+1, byrow = T)
  PCp3.matrix = matrix(PCp3, K.max+1, R.max+1, byrow = T)
  ICp3.matrix = matrix(ICp3, K.max+1, R.max+1, byrow = T)
  
  
  if (1==1)
  {
    rownames = c()
    colnames = c()
    for (K in 0:K.max)
    {
      rownames = c(rownames, paste("K=", K, sep=""))
    }
    for (R in 0:R.max)
    {
      colnames = c(colnames, paste("R=", R, sep=""))
    }
    
    rownames(MSS.matrix) = rownames
    rownames(PCp1.matrix) = rownames
    rownames(ICp1.matrix) = rownames
    rownames(PCp2.matrix) = rownames
    rownames(ICp2.matrix) = rownames
    rownames(PCp3.matrix) = rownames
    rownames(ICp3.matrix) = rownames

    colnames(MSS.matrix) = colnames
    colnames(PCp1.matrix) = colnames
    colnames(ICp1.matrix) = colnames
    colnames(PCp2.matrix) = colnames
    colnames(ICp2.matrix) = colnames
    colnames(PCp3.matrix) = colnames
    colnames(ICp3.matrix) = colnames
  }
  
  
  best.PCp1 = fit[[KR.PCp1]]
  best.ICp1 = fit[[KR.ICp1]]
  best.PCp2 = fit[[KR.PCp2]]
  best.ICp2 = fit[[KR.ICp2]]
  best.PCp3 = fit[[KR.PCp3]]
  best.ICp3 = fit[[KR.ICp3]]
   
  
  cat("the selected K from PCp1 is:", K.PCp1, "\n")
  cat("the selected R from PCp1 is:", R.PCp1, "\n")
  cat("the selected K from ICp1 is:", K.ICp1, "\n")
  cat("the selected R from ICp1 is:", R.ICp1, "\n")
  cat("the selected K from PCp2 is:", K.PCp2, "\n")
  cat("the selected R from PCp2 is:", R.PCp2, "\n")
  cat("the selected K from ICp2 is:", K.ICp2, "\n")
  cat("the selected R from ICp2 is:", R.ICp2, "\n")
  cat("the selected K from PCp3 is:", K.PCp3, "\n")
  cat("the selected R from PCp3 is:", R.PCp3, "\n")
  cat("the selected K from ICp3 is:", K.ICp3, "\n")
  cat("the selected R from ICp3 is:", R.ICp3, "\n")
  
  
  
  return(list(fit.collection=fit, 
    
              MSS=MSS.matrix, 
              ICp1=ICp1.matrix, PCp1=PCp1.matrix, 
              ICp2=ICp2.matrix, PCp2=PCp2.matrix, 
              ICp3=ICp3.matrix, PCp3=PCp3.matrix, 
              
              K.ICp1=K.ICp1, R.ICp1=R.ICp1, K.PCp1=K.PCp1, R.PCp1=R.PCp1, 
              K.ICp2=K.ICp2, R.ICp2=R.ICp2, K.PCp2=K.PCp2, R.PCp2=R.PCp2, 
              K.ICp3=K.ICp3, R.ICp3=R.ICp3, K.PCp3=K.PCp3, R.PCp3=R.PCp3)
         )
}