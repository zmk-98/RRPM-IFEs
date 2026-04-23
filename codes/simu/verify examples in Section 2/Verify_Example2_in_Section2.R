#' @example the need for rank-one decomposition
#' 
#' @model y_{it} = (X_{it,1}/\sqrt{2} + X_{it,2}/\sqrt{2})*b_{1i} + X_{it,3}*b_{2i} + F_{t1}*L_{i1} + F_{t2}*L_{i2} + \sigma*\epsilon_{it},
#' where X_1 and X_2 are both high-rank, X_3 is low-rank and Rank(X_3)=2, and X_3=v_1*w_1^\top+v_2*w_2^\top, F_1=v_2, F_2=v_3
#' 
#' @case For the case of sigma = 0:
#' @description We can't explain the model with the orthogonality restriction, when we don't do the rank-one decomposition, 
#' i.e., fit1$RSS>0 when K=2 and R=2
#' @description Even if K equals p and R is a fixed large number, fit2$RSS>0 
#' @description We successfully fit the original model again with the orthogonality restriction, when we do the rank-one decomposition, 
#' i.e., fit3$RSS=0 when K=3 and R=1
#' 
#' @case For the case of sigma = 0.2:
#' @description fit1$MSE and fit2$MSE obviously larger than fit3$MSE
#' 
#' @conclusion This shows the need for rank-one decomposition with the orthogonality restriction

setwd("codes/simu/verify examples in Section 2")

library("foreach"); library("doParallel"); cl <- makeCluster(36); registerDoParallel(cl) 

source("RPC.R")


N.simu = 100
OUTALL = c()

sigma = 0
NTi.list = list(c(30,30), c(50,50), c(70,70), c(100,100), c(150,150), c(200,200))
p = 3; K = 2; R = 2
p.a = 4; K.a = 3; R.a = 1

for (i in 1:6) {
  NTi = NTi.list[[i]]
  N = NTi[1]
  Ti = NTi[2]
  
  OUT <- foreach(iter=1:N.simu, .combine='rbind') %dopar%
  {
    M = sqrt(Ti) * svd(matrix(rnorm(Ti*4), Ti, 4))$u
    X1 = matrix(rnorm(Ti*N), Ti, N)
    X2 = matrix(rnorm(Ti*N), Ti, N)
    w1 = runif(N) + 0.5
    w2 = runif(N) + 0.5
    X3 = M[,1] %*% t(w1) + M[,2] %*% t(w2)
    F0 = cbind(M[,2]+M[,3], M[,4])
    
    Theta = cbind(c(1/sqrt(2), 1/sqrt(2), 0), c(0, 0, 1))
    B = cbind(runif(N)+0.5, runif(N)+0.5)
    L = cbind(runif(N)+0.5, runif(N)+0.5)
  
    index = c()
    y = c()
    X = c()
    for (i in 1:N) {
      index = c(index, rep(i, Ti))
      Xi = cbind(X1[,i], X2[,i], X3[,i])
      yi = Xi %*% Theta %*% B[i,] + F0 %*% L[i,] +  sigma * rnorm(Ti)
      
      X = rbind(X, Xi)
      y = c(y, yi)
    }
    
    X3.1 = M[,1] %*% t(w1)
    X3.2 = M[,2] %*% t(w2)
    # svd.X3 = svd(X3) # alternatively, one can use svd to do rank-one decomposition
    # X3.1 = svd.X3$d[1] * svd.X3$u[,1] %*% t(svd.X3$v[,1])
    # X3.2 = svd.X3$d[2] * svd.X3$u[,2] %*% t(svd.X3$v[,2])
    Xnew = cbind(as.vector(X1), as.vector(X2), as.vector(X3.1), as.vector(X3.2))
    
    v = cbind(M[,1], M[,2])
    v.proj = v %*% solve(t(v)%*%v) %*% t(v)
    v.proj.perp = diag(1, Ti) - v.proj
    
    # fit1 = RPC(X=X, y=y, v=NULL, index1=index, factor1=2, factor2=2, max_iter=200)
    # print(round(fit1$RSS/(N*Ti), 4)) # if sigma = 0, RSS = 0
    
    fit1 = RPC(X=X, y=y, v=v, index1=index, factor1=2, factor2=2, max_iter=200)
    print(round(fit1$RSS)/(N*Ti), 4) # if sigma = 0, RSS is positive
    
    fit2 = RPC(X=X, y=y, v=v, index1=index, factor1=3, factor2=10, max_iter=200)
    print(round(fit2$RSS/(N*Ti), 4)) # if sigma = 0, RSS is positive
    
    fit3 = RPC(X=Xnew, y=y, v=v, index1=index, factor1=3, factor2=2, max_iter=200)
    print(round(fit3$RSS/(N*Ti), 4)) # if sigma = 0, RSS = 0
    
    out = c(fit1$RSS/(N*Ti), fit2$RSS/(N*Ti), fit3$RSS/(N*Ti)) 
    out
  }
  
  out = c(sigma, N, Ti, round(colMeans(OUT), 4))
  
  print(out)
  
  OUTALL = rbind(OUTALL, out)
  colnames(OUTALL) = c("sigma", "N", "Ti", "fit1", "fit2", "fit3")
  print(OUTALL)
  
  write.csv(OUTALL, file = "Verify_Example2.csv")
}

