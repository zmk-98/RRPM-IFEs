#' @example loss of parsimony
#' 
#' @model y_{it} = (X_{it,1}+X_{it,2}+X_{it,3})*b_i + F_t*L_i + \sigma*\epsilon_{it},
#' where all predictors are cross-sectional invariant (of rank 1), or X_{i,1}=v_1, X_{i,2}=v_2, X_{i,3}=v_3, and F=v_3+v_4
#' 
#' @case For the case of sigma = 0:
#' @description The original model is excluded from the model with the orthogonality restriction, if we keep K=1 unchanged
#' that is, fit1$MSE>0 when K=1 and R=1
#' @description Even if R is a fixed large number, fit2$MSE>0 when K=1
#' @description The original model is included in the model with the orthogonality restriction again, if we try a larger K, 
#' that is, fit3$MSE=0 when K=2 and R=1
#' 
#' @case For the case of sigma = 0.2:
#' @description fit1$MSE and fit2$MSE obviously larger than fit3$MSE
#' 
#' @conclusion A choice of larger K under the additional orthogonality restriction shows the loss of parsimony 

setwd("codes/simu/verify examples in Section 2")

library("foreach"); library("doParallel"); cl <- makeCluster(36); registerDoParallel(cl) 

source("RPC.R")


N.simu = 2
OUTALL = c()

sigma = 0.2
# sigma = 0.2
NTi.list = list(c(30,30), c(50,50), c(70,70), c(100,100), c(150,150), c(200,200))
p = 3; K = 1; R = 2
p.a = 3; K.a = 2; R.a = 2

for (i in 1:6) {
  NTi = NTi.list[[i]]
  N = NTi[1]
  Ti = NTi[2]
  
  OUT <- foreach(iter=1:N.simu, .combine='rbind') %dopar%
  {
    M = sqrt(Ti) * svd(matrix(rnorm(Ti*5), Ti, 5))$u
    X1 = matrix(rep(M[,1], N), Ti, N)
    X2 = matrix(rep(M[,2], N), Ti, N)
    X3 = matrix(rep(M[,3], N), Ti, N)
    F0 = cbind(M[,3] + M[,4], M[,5])
    
    Theta = c(1, 1, 1)
    # Theta.a = cbind(c(1, 1, 1), c(0, 0, 1))
    
    B = runif(N) + 0.5
    L = matrix(runif(N*2) + 0.5, N, 2)
    
    index = c()
    y = c()
    X = c()
    for (i in 1:N) {
      index = c(index, rep(i, Ti))
      Xi = cbind(X1[,i], X2[,i], X3[,i])
      yi = Xi %*% Theta * B[i] + F0 %*% L[i,] +  sigma * rnorm(Ti)
      
      X = rbind(X, Xi)
      y = c(y, yi)
    }
    
    v = cbind(X1[,1], X2[,1], X3[,1])
    v.proj = v %*% solve(t(v)%*%v) %*% t(v)
    v.proj.perp = diag(1,Ti) - v.proj
    
    # fit1 = RPC(X=X, y=y, v=NULL, index1=index, factor1=1, factor2=1, max_iter=200) 
    # print(round(fit1$RSS/(N*Ti), 4)) 
    # if sigma = 0, RSS1 = 0
    
    fit1 = RPC(X=X, y=y, v=v, index1=index, factor1=1, factor2=2, max_iter=200)
    print(round(fit1$RSS/(N*Ti), 4)) 
    # if sigma = 0, RSS1 is positive
    
    fit2 = RPC(X=X, y=y, v=v, index1=index, factor1=1, factor2=10, max_iter=200)
    print(round(fit2$RSS/(N*Ti), 4)) 
    # if sigma = 0, RSS2 is positive
    
    fit3 = RPC(X=X, y=y, v=v, index1=index, factor1=2, factor2=2, max_iter=200)
    print(round(fit3$RSS/(N*Ti), 4)) 
    # if sigma = 0, RSS3 = 0
    
    out = c(fit1$RSS/(N*Ti), fit2$RSS/(N*Ti), fit3$RSS/(N*Ti)) 
    out
  }
  
  out = c(sigma, N, Ti, round(colMeans(OUT), 4))
  
  print(out)
  
  OUTALL = rbind(OUTALL, out)
  colnames(OUTALL) = c("sigma", "N", "Ti", "fit1", "fit2", "fit3")
  print(OUTALL)
  
  write.csv(OUTALL, file = "Verify_Example3_sigma0.csv")
}


