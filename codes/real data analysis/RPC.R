#' Estimate a Reduced-Rank Panel Model (RRPM) with Interactive Fixed Effects (IFEs) by Restricted Principal Component (RPC) approach
#' 
#' @param X: an NT * p matrix, predictors
#' @param v: an T * p2 matrix, low-rank predictors
#' @param y: a vector of length NT, responses
#' @param index1: a vector of length NT, show the section of each observation
#' @param factor1: K, the number of linear combinations, X*theta_k
#' @param factor2: R, the number of latent factors
#' 
#' @return A list containing estimated (Theta, B, F, L)

RPC = function(X, y, v = NULL,
                 factor1 = 1, factor2 = 0,
                 index1 = c(),
                 lambda = ncol(X)/length(y), max_iter = 200)
{

  XX = X
  YY = y
  
  p = ncol(X)
  n = length(y)
  K = factor1
  R = factor2
  
  ind = unique(index1)
  N = length(ind)
  # assume balanced panel data: each observation has a time series of the same length
  Ti = n/N
  
  lambda = lambda
  
  # p1 - NO. of high-rank predictors
  # p2 - NO. of low-rank predictors
  if (is.null(v)) {
    p2 = 0
  } else {
    p2 = dim(v)[2]
  }
  p1 = p - p2
  
  if (p2 >= 1) {
    v.proj = v %*% solve(t(v) %*% v) %*% t(v)
    v.proj.perp = diag(1, Ti) - v.proj
  } else {
    v.proj.perp = diag(1, Ti)
  }
  
  # YY = YY - mean(YY)
  if (N>0) {
    X = list()
    y = list()
    I = list()
    
    for (i in 1:N) { 
      I[[i]] = which(index1==ind[i])
      X[[i]] = XX[I[[i]],]
      y[[i]] = YY[I[[i]]]

#      X[[i]] = X[[i]] - matrix(colMeans(X[[i]]), Ti, ncol(X[[i]]), byrow=TRUE)
#      y[[i]] = y[[i]]-mean(y[[i]])
      
#      YY[I[[i]]] = y[[i]]
#      XX[I[[i]], ] = X[[i]]
      
    }
  } 
  
  
  Theta = NULL
  B = NULL
  F0 = NULL
  L0 = NULL
  
  # null model
  if (K==0 & R==0) {
    resi = YY
    RSS = sum(YY^2)
  }
  # print(RSS)
  
  # pure factor model
  if (K==0 & R>0) { 
    OP = matrix(0, Ti, Ti)
    for (i in 1:N) {
      OP = OP + y[[i]] %*% t(y[[i]])
    }
    
    F0 = sqrt(Ti) * as.matrix(eigen(v.proj.perp %*% OP %*% v.proj.perp)$vectors[,1:R])
    
    L0 = matrix(0, N, R)
    resi = rep(0, n)
    for (i in 1:N) {
      L0[i,] = t(F0) %*% y[[i]]/Ti
      resi[I[[i]]] = y[[i]] - F0 %*% L0[i,]
    }
    RSS = sum(resi^2)
  }

  # RRPM (K>0) with interactive effects (R>0) or without interactive effects (R=0)
  if (K>0) { 
    # starting value for (Theta, B)
    Theta = matrix(0, p, K)
    Beta = c()
    S = eigen(cov(XX))
    S = S$vectors %*% diag(1/sqrt(S$values+1/n)) %*% t(S$vectors)
    
    Xinv = XX%*%S
    
    for (i in 1:N) {
      Xi = Xinv[I[[i]],]
      beta_i = solve(t(Xi)%*%Xi+ diag(lambda,p))%*%(t(Xi)%*%y[[i]])
      Beta = cbind(Beta, beta_i) 
    }
    Theta = as.matrix(S %*% eigen(Beta%*%t(Beta))$vectors[,1:K])

    for (k in 1:K) {
      Theta[,k] = Theta[,k]/sqrt(sum(Theta[,k]^2))
    }
    
    B = matrix(0, N, K)
    r1 = rep(0, n)
    for (i in 1:N) {
      Xi.Theta = as.matrix(X[[i]]%*%Theta)
      B[i,] = solve(t(Xi.Theta)%*%Xi.Theta + diag(lambda,K)) %*% (t(Xi.Theta)%*%y[[i]])
      r1[I[[i]]] = y[[i]] - Xi.Theta%*%B[i,]
    }
    RSS = sum(r1^2)
    

    # iteration between (Theta, B) and (F, L)
    for (iter1 in 1:max_iter) { 
      # update (F, L)
      Z = rep(0, n)
      if (R>0) {
        OP = matrix(0, Ti, Ti)
        for (i in 1:N) {
          OP = OP + r1[I[[i]]]%*%t(r1[I[[i]]])
        }
        F0 = sqrt(Ti) * as.matrix(eigen(v.proj.perp %*% OP %*% v.proj.perp)$vectors[,1:R])
        L0 = matrix(0, N, R)
        for (i in 1:N) {
          L0[i,] = t(F0) %*% r1[I[[i]]]/Ti
          Z[I[[i]]] = y[[i]] - F0 %*% L0[i,]          
        }
      } else {
        Z = YY
      }
      
      # update (Theta, B)
      for (k in 1:K) { 
        theta.k.0 = Theta[,k]
        r1 = rep(0, n)
        B = matrix(0, N, K)
        for (i in 1:N) {
          # update B
          Xi.Theta = X[[i]]%*%Theta
          B[i,] = solve(t(Xi.Theta)%*%Xi.Theta+diag(lambda,K)) %*% t(Xi.Theta)%*%Z[I[[i]]]
          r1[I[[i]]] = Z[I[[i]]] - Xi.Theta%*%B[i,] + Xi.Theta[,k]*B[i,k]
        }
        
        # update Theta
        XX.b = XX
        lambda1 = 0
        for (i in 1:N) {
          XX.b[I[[i]],] = XX[I[[i]],]*B[i,k]
          lambda1 = lambda1 + B[i,k]^2 
        }
        theta.k = solve(t(XX.b)%*%XX.b+diag(lambda1*lambda/N, p)) %*% t(XX.b) %*% r1
#        theta.k = solve(t(XX.b)%*%XX.b+diag(lambda, p)) %*% t(XX.b)%*% r1
        B[,k] = B[,k]*sqrt(sum(theta.k^2))
        theta.k = theta.k/sqrt(sum(theta.k^2))
        
        Theta[,k] = theta.k
      }
      
      r1 = rep(0, n)
      resi = rep(0, n)
      B = matrix(0, N, K)
      for (i in 1:N) { 
        Xi.Theta = X[[i]]%*%Theta
        B[i,] = solve(t(Xi.Theta)%*%Xi.Theta+diag(lambda,K)) %*% t(Xi.Theta)%*%Z[I[[i]]]
        r1[I[[i]]] = y[[i]] - Xi.Theta%*%B[i,]
        resi[I[[i]]] = Z[I[[i]]] - Xi.Theta%*%B[i,]
      }
      
      if(abs(RSS-sum(resi^2)) < 1.0e-6)
        break
      
      RSS = sum(resi^2)
      # cat(iter1, "\n")
      # print(RSS)
    }
  }
  
  resi = rep(0, n)
  if (K>0 & R>0)
  {
    OP = matrix(0, Ti, Ti)
    for (i in 1:N)
    {
      OP = OP + r1[I[[i]]] %*% t(r1[I[[i]]])
    }
    F0 = sqrt(Ti) * as.matrix(eigen(v.proj.perp %*% OP %*% v.proj.perp)$vectors[,1:R])
    L0 = matrix(0, N, R)
    for (i in 1:N)
    {
      L0[i,] = t(F0) %*% r1[I[[i]]]/Ti
      resi[I[[i]]] = r1[I[[i]]] - F0 %*% L0[i,]
    }
    RSS = sum(resi^2)
  }
  
  return(list(resi=resi, RSS=RSS, hat.Theta=Theta, hat.B=B, hat.F=F0, hat.L=L0, K=K, R=R))
}

