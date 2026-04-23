#' Bias Correction for Restricted Principal Component (RPC) Estimator
#' 
#' @description We only consider the heteroskedasticity in variance, and not consider the serial or cross-sectional correlation
#' @details Special case - If K=0, no need for bias correction of hat Theta
#' @details Special case - If R=0,bias xi is 0, the second term of bias zeta is 0 in theory
#' @details Special case - Serial homogeneity in variance: bias zeta is 0 in theory
#' @details Special case - Cross-sectional homogeneity in variance: bias xi is 0 in theory
#' @details Special case - X and epsilon are independent: bias eta is 0 in theory
#' 
#' @param X: an NT * p matrix, predictors
#' @param v: a T * p2 matrix, low-rank predictors
#' @param y: a vector of length NT, responses  
#' @param index1: a vector of length NT, show the section of each observation
#' @param factor1: K, the number of linear combinations, X*theta_k (K > 0)
#' @param factor2: R, the number of latent factors
#' @param m: the bandwidth for estimated Nickell bias
#' @param (Theta, B, F0, L0): true parameters
#' @param E: a T * N matrix, realization of errors
#' @param Omega: a T * N matrix, variance of errors
#' @param RPC.return: return from function RPC(...) when (K,R) are correctly specified, (Theta, B, F, L)
#' 
#' @return bias corrected RPC estimator
#' 

setwd("codes/simu/bias correction")

library("foreach"); library("doParallel"); cl <- makeCluster(36); registerDoParallel(cl) 

# source("choose_KR.R")
source("RPC.R")
source("BiasCorrectedRPC.R")

p = 3
K0 = 1 # reduced-rank
R0 = 1 # latent factors
NTi.list = list(c(30,30), c(50,50), c(70,70), c(100,100), c(150,150), c(200,200))

error.type = 1 
# 1: IID normal; 2: heavy tailed (t-distribution); 3: heteroskedastic; 

OUTALL = c()

N.simu = 100

for(i.NTi in 1:6)
{
  NTi = NTi.list[[i.NTi]]; N = NTi[1]; Ti = NTi[2]
  
  OUT <- foreach(iter=1:N.simu, .combine='rbind') %dopar%
  {
    set.seed(iter)
    
    index = c()
    X = c()
    y = c()
    e = c()

    # generating the data (Xi,yi)
    if (K0 > 0) {
      Theta = svd(matrix(rnorm((p-1) * K0), p-1, K0))$u
    } else {
      stop("K=0: No need for bias correction.")
    }
    
    rho.F = 0.5
    sigma.F = 0.5
    if (R0 > 0) {
      L = matrix(runif(N*R0)+0.5, N, R0)
      F1 = matrix(NA, Ti+100, R0)
      F1[1,] = rnorm(R0)
      for (t in 1:(Ti+99)) {
        F1[t+1,] = rho.F * F1[t,] + rnorm(R0, mean=0, sd=sqrt(1-rho.F^2)*sigma.F)
      }
      F0 = as.matrix(F1[101:(Ti+100),])
    } else {
      F1 = matrix(0, Ti+100, 1)
      F0 = matrix(0, Ti, 1)
      L = matrix(0, N, 1)
    }
    
    if (R0 > 0) {
      F0.proj.perp = diag(1, Ti) - F0 %*% solve(t(F0)%*%F0) %*% t(F0)
    } else {
      F0.proj.perp = diag(1, Ti)
    }
    
    v1 = rnorm(Ti+100)
    Z = F0.proj.perp %*% v1[101:(Ti+100)]
    v1[101:(Ti+100)] = Z * norm(v1[101:(Ti+100)],"2") / norm(Z,"2")
    v = matrix(v1[101:(Ti+100)], Ti, 1)
    
    B = matrix(0, N, K0)
    Err = matrix(0, Ti+100, N)
    E = matrix(0, Ti, N)
    Omega = matrix(0, Ti, N)
    rho.y = 0.5
    
    for (i in 1:N) {
      bk = 0.5 + runif(K0)
      B[i,] = bk
      index = c(index, rep(i, Ti))
      if (error.type == 1) {
        Err[,i] = 2 * rnorm(Ti+100)
      } else if (error.type == 2) {
        Err[,i] = 2 * rt(Ti+100, df=10)
      } else if (error.type == 3) {
        Err[,i] = 2 * rnorm(Ti+100)
        Odd = ((1:(Ti+100))%%2==1)
        Even = ((1:(Ti+100))%%2==0)
        Err[Odd, i] = sqrt(2/3) * Err[Odd, i]
        Err[Even, i] = sqrt(4/3) * Err[Even, i]
      }
      
      y1 = c()
      X1 = c()
      y1[1] = 0
      X1 = rbind(X1, c(y1[1], rnorm(p-2), v1[1]))
      for (t in 2:(Ti+100)) {
        X1 = rbind(X1, c(y1[t-1], rnorm(p-2), v1[t]))
        y1[t] = y1[t-1] * rho.y * bk[1] + X1[t,2:p] %*% Theta %*% bk + t(F1[t,])%*%L[i,] + Err[t,i]
        y1 = c(y1, y1[t])
      }
      X = rbind(X, X1[101:(Ti+100),])
      y = c(y, y1[101:(Ti+100)])
    }
    if (K0 > 1) {
      Theta0 = rbind(c(rho.y, 0), Theta)
    } else {
      Theta0 = matrix(c(rho.y, Theta), p, 1)
    }
    
    fit.oracle = RPC(X=X, y=y, v=v, index1 = index, 
                     factor1=K0, factor2=R0, max_iter=500)
    
    err.Theta = 0 # measure the error of hat Theta
    err.Theta.BC = 0 # measure the error of bias-corrected hat Theta
    
    if (K0 > 0) {
      Theta0.proj = Theta0 %*% solve(t(Theta0)%*%Theta0) %*% t(Theta0)
      
      hat.Theta = fit.oracle$hat.Theta
      hat.Theta.proj = hat.Theta %*% solve(t(hat.Theta)%*%hat.Theta) %*% t(hat.Theta)
      err.Theta = norm(hat.Theta.proj-Theta0.proj, "f")
      
      BC.Theta.proj.vec = bias_corrected_RPC(X=X, y=y, v=v, factor1=K0, factor2=R0,
                                             m = ceiling(1.5*Ti^(3/8)), index1 = index, 
                                             RPC.return = fit.oracle)
      err.Theta.BC = norm(BC.Theta.proj.vec-as.vector(Theta0.proj), "2")
    }
    
    out = c(err.Theta, err.Theta.BC)
    out
  }
  
  
  # estimation error results: err.Theta, err.Theta.BC, err.F
  out = c(K0, R0, N, Ti, round(colMeans(OUT),4))
  
  print(out)
  
  OUTALL = rbind(OUTALL, out)
  colnames(OUTALL) = c("K0", "R0", "N", "Ti", "err.Theta", "err.Theta.BC")
  print(OUTALL)
  
  write.csv(OUTALL, file = "RPCerr1_OUT.csv")
}

# K0 R0   N  Ti err.Theta err.Theta.BC
# out  1  1  30  30    0.1213       0.1160
# out  1  1  50  50    0.0646       0.0615
# out  1  1  70  70    0.0479       0.0399
# out  1  1 100 100    0.0283       0.0265
# out  1  1 150 150    0.0146       0.0129

