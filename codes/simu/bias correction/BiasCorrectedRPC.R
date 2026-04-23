#' Bias Correction for Restricted Principal Component (RPC) Estimator
#' 
#' @description We only consider the heteroskedasticity in variance, and not consider the serial or cross-sectional correlation
#' @details Special case - If K=0, no need for bias correction of hat Theta
#' @details Special case - If R=0, bias xi is 0, the second term of bias zeta is 0
#' @details Special case - Serial homogeneity in variance: bias zeta is 0
#' @details Special case - Cross-sectional homogeneity in variance: bias xi is 0
#' @details Special case - X and epsilon are independent: bias eta is 0
#' 
#' @param X: an NT * p matrix, predictors
#' @param v: an T * p2 matrix, low-rank predictors
#' @param y: a vector of length NT, responses  
#' @param index1: a vector of length NT, show the section of each observation
#' @param factor1: K, the number of linear combinations, X*theta_k (K > 0)
#' @param factor2: R, the number of latent factors
#' @param m: the bandwidth for estimated Nickell bias
#' @param RPC.return: return from function RPC(...) when (K,R) are correctly specified, (Theta, B, F, L)
#' 
#' @return bias corrected RPC estimator
#' 

source("RPC.R")

bias_corrected_RPC <- function(X, y, v = NULL, 
                               factor1 = 1, factor2 = 0, m = 5,
                               index1 = c(),
                               lambda = ncol(X)/length(y), RPC.return = RPC.return)
{ 
  XX = X
  YY = y
  
  p = ncol(X)
  n = length(y)
  K = factor1
  R = factor2
  
  if (K==0)
  {
    stop("K=0: No need for bias correction.")
  }
  
  ind = unique(index1)
  N = length(ind)
  # assume balanced panel data: each observation has a time series of the same length
  Ti = n/N
  
  lambda = lambda
  
  # p1 - NO. of high-rank predictors
  # p2 - NO. of low-rank predictors
  if (is.null(v)) {
    p2 = 0
    v = matrix(0, Ti, 1)
    Sigma.v.inv = 0
  } else {
    p2 = dim(v)[2]
    Sigma.v.inv = solve(t(v) %*% v / Ti)
  }
  p1 = p - p2
  
  v.proj = v %*% Sigma.v.inv %*% t(v) / Ti
  v.proj.perp = diag(1, Ti) - v.proj
  
  
  # YY = YY - mean(YY)
  if (N > 0) {
    X = list()
    y = list()
    I = list()
    
    for (i in 1:N) { 
      I[[i]] = which(index1==ind[i])
      X[[i]] = XX[I[[i]],]
      y[[i]] = YY[I[[i]]]
    }
  }
  
  Theta = RPC.return$hat.Theta
  B = RPC.return$hat.B
  hat.F = RPC.return$hat.F
  L = RPC.return$hat.L
  
  if (R > 0) {
    Sigma.L.inv = solve(t(L) %*% L / N)
  } else { 
    L = matrix(0, N, 1)
    Sigma.L.inv = 0
    hat.F = matrix(0, Ti, 1)
  }
  W = L %*% Sigma.L.inv %*% t(L)
  F.proj = hat.F %*% t(hat.F) / Ti
  F.proj.perp = diag(1, Ti) - F.proj
  Theta.proj = Theta %*% solve(t(Theta)%*%Theta) %*% t(Theta)
  
  J = cbind(matrix(0, p-K, K), diag(1, p-K)) # (p-K)*p
  beth = matrix(0, Ti*max(R,1), Ti*max(R,1)) # TR*TR
  A2 = matrix(0, Ti*max(R,1), (p-K)*K)
  E = matrix(0, Ti, N) # residual matrix
  
  A1 = list()  
  A1.sqr = list()
  H = list()
  Q = list()
  S = list()
  U = list()
  V = list()
  Z = list() # transformed predictors
  
  for (i in 1:N) {
    A1[[i]] = F.proj.perp %*% X[[i]] %*% Theta
    A1.sqr[[i]] = t(A1[[i]]) %*% A1[[i]]
    H[[i]] = A1[[i]] %*% solve(A1.sqr[[i]]) %*% t(A1[[i]])
    beth = beth + (1/N) * kronecker(L[i,] %*% t(L[i,]), diag(1,Ti) - v.proj.perp %*% H[[i]] %*% v.proj.perp) 
    Q[[i]] = (diag(1,Ti) - H[[i]]) %*% F.proj.perp %*% X[[i]] %*% t(J)
    S[[i]] = (diag(1,p) - Theta %*% solve(A1.sqr[[i]]) %*% t(A1[[i]])%*%X[[i]]) %*% t(J)
    V[[i]] = X[[i]] %*% Theta %*% solve(A1.sqr[[i]]) %*% t(A1[[i]])
    U[[i]] = diag(1,Ti) - V[[i]]
    A2 = A2 + (1/N) * kronecker(L[i,] %*% t(B[i,]), v.proj.perp %*% Q[[i]])
    E[,i] = y[[i]] - X[[i]] %*% Theta %*% B[i,] - hat.F %*% L[i,]
  }
  if (R > 0) {
    beth.inv = solve(beth)
  } else {
    beth.inv = matrix(0, Ti, Ti)
  }
  
  D0 = matrix(0, (p-K)*K, (p-K)*K)
  DZ = matrix(0, (p-K)*K, (p-K)*K)
  for (i in 1:N) {
    Z[[i]] = kronecker(t(B[i,]), Q[[i]]) - kronecker(t(L[i,]), diag(1,Ti)-H[[i]]) %*% beth.inv %*% A2
    D0 = D0 + t(Z[[i]]) %*% Z[[i]]
    for (t in 1:Ti) {
      DZ = DZ + E[t,i]^2 * Z[[i]][t,] %*% t(Z[[i]][t,])
    }
  }
  D0 = D0 / (N*Ti)
  DZ = DZ / (N*Ti)
  
  A0 = (diag(1,p^2) + commutation_matrix(p,p)) %*% kronecker(Theta %*% solve(t(Theta) %*% Theta), (diag(1,p) - Theta.proj) %*% t(J)) # p^2 * K(p-K) matrix 
  
  
  # incidental parameter bias due to (B, L)
  zeta1 = rep(0, (p-K)*K)
  zeta2 = matrix(0, Ti, max(R,1))
  for (i in 1:N) {
    zeta1 = zeta1 + (1/N) * as.vector(t(Q[[i]]) %*% diag(E[,i]^2) %*% A1[[i]] %*% solve(A1.sqr[[i]]))
    zeta2 = zeta2 + (1/N) * F.proj.perp %*% U[[i]] %*% diag(E[,i]^2) %*% t(U[[i]]) %*% hat.F
  }
  bias.zeta1 = A0 %*% solve(D0) %*% zeta1
  bias.zeta2 = -A0 %*% solve(D0) %*% ((1/Ti) * t(A2) %*% beth.inv %*% as.vector(zeta2))
  
  
  # incidental parameter bias due to F
  A3 = list()
  A4 = list()
  A5 = list()
  A6 = list()
  for (i in 1:N) {
    A3[[i]] = (1/Ti) * kronecker(B[i,], t(S[[i]]) %*% t(X[[i]]) %*% hat.F) # (p-K)K*R
    A4[[i]] = Ti^(-1/2) * kronecker(L[i,], t(V[[i]]) %*% hat.F) # TR*R
  }
  for (i in 1:N) { 
    A5[[i]] = matrix(0, (p-K)*K, max(R,1))
    A6[[i]] = matrix(0, Ti*max(R,1), max(R,1))
    for (k in 1:N) {
      A5[[i]] = A5[[i]] + W[i,k] * A3[[k]]
      A6[[i]] = A6[[i]] + W[i,k] * A4[[k]]
    }
    A5[[i]] = A5[[i]]/N; A6[[i]] = A6[[i]]/N
  }
  xi1 = rep(0, (p-K)*K)
  xi2 = rep(0, (p-K)*K)
  xi3 = rep(0, Ti*max(R,1))
  xi4 = rep(0, Ti*max(R,1))
  for (i in 1:N) { 
    xi1 = xi1 + mean(E[,i]^2) * A3[[i]] %*% Sigma.L.inv %*% L[i,] 
    xi2 = xi2 + mean(E[,i]^2) * A5[[i]] %*% Sigma.L.inv %*% L[i,]
    xi3 = xi3 + mean(E[,i]^2) * A4[[i]] %*% Sigma.L.inv %*% L[i,]
    xi4 = xi4 + mean(E[,i]^2) * A6[[i]] %*% Sigma.L.inv %*% L[i,]
  }
  xi1 = xi1/N; xi2 = xi2/N; xi3 = xi3/N; xi4 = xi4/N
  bias.xi1 = -A0 %*% solve(D0) %*% xi1
  bias.xi2 = A0 %*% solve(D0) %*% xi2
  bias.xi3 = -Ti^(-1/2) * A0 %*% solve(D0) %*% t(A2) %*% beth.inv %*% xi3
  bias.xi4 = Ti^(-1/2) * A0 %*% solve(D0) %*% t(A2) %*% beth.inv %*% xi4
  
  
  # Nickell bias
  eta1 = rep(0, K*(p-K))
  eta2 = rep(0, K*(p-K))
  
  A7 = list() # T^{-1} X_i^\top \hat F
  A8 = list() # T^{-1} X_i^\top v
  A9 = list() # \sum_{s=t+1}^{t+M} X_{is} X_{is}^\top
  A10 = list() # \sum_{s=t+1}^{t+M} X_{is} \hat F_s^\top
  A11 = list() # \sum_{s=t+1}^{t+M} X_{is} v_s^\top
  
  for(i in 1:N) {
    A7[[i]] = (1/Ti) * t(X[[i]]) %*% hat.F
    A8[[i]] = (1/Ti) * t(X[[i]]) %*% v
    A9[[i]] = list()
    A10[[i]] = list()
    A11[[i]] = list()
    for (t in 1:(Ti-1)) { 
      if (Ti-t >= 2) {
        A9[[i]][[t]] = t(X[[i]][(t+1):min(t+m,Ti),]) %*% X[[i]][(t+1):min(t+m,Ti),] # p * p
        A10[[i]][[t]] = t(X[[i]][(t+1):min(t+m,Ti),]) %*% hat.F[(t+1):min(t+m,Ti),] # p * max(R,1)
        A11[[i]][[t]] = t(X[[i]][(t+1):min(t+m,Ti),]) %*% v[(t+1):min(t+m,Ti),] # p * max(p2,1)
      } else {
        A9[[i]][[t]] =  X[[i]][Ti,] %*% t(X[[i]][Ti,])
        A10[[i]][[t]] = X[[i]][Ti,] %*% t(hat.F[Ti,])
        A11[[i]][[t]] = X[[i]][Ti,] %*% t(v[Ti,])
      }
    }
  }
  
  for (i in 1:N) {
    for (t in 1:(Ti-1)) {
      # Temp is an amount associated with (i,t)
      Temp = (
        -J %*% (A9[[i]][[t]] - A10[[i]][[t]] %*% t(A7[[i]]) - A7[[i]] %*% t(A10[[i]][[t]])) +
        J %*% (t(X[[i]]) %*% X[[i]]/Ti - A7[[i]] %*% t(A7[[i]])) %*% 
          Theta %*% solve(A1.sqr[[i]]/Ti) %*% t(Theta) %*% (A9[[i]][[t]] - A10[[i]][[t]] %*% t(A7[[i]]) - A7[[i]] %*% t(A10[[i]][[t]]))
        ) %*% Theta %*% solve(A1.sqr[[i]]/Ti) %*% t(Theta) 
      ###
      eta1 = eta1 + 1/(N*Ti) * kronecker(
        B[i,], 
        Temp %*% X[[i]][t,] * E[t,i]
        )
      eta2 = eta2 - 1/(N*Ti) * kronecker(
        B[i,], 
        (Temp %*% A7[[i]] + t(S[[i]]) %*% A10[[i]][[t]]) %*% hat.F[t,] * E[t,i]
        )
    }
  }
  
  
  eta3 = rep(0, K*(p-K))
  eta4 = rep(0, K*(p-K))
  
  A12 = list() # list of length Ti \sum_{j=1}^N ...  
  A13 = list() # list (of length N) of list (of length Ti) 
  
  for (s in 1:Ti) {
    A12[[s]] = matrix(0, K*(p-K), max(R,1))
    for (j in 1:N) {
      A12[[s]] = A12[[s]] + (1/N) * kronecker(
        B[j,] %*% t(L[j,]) %*% Sigma.L.inv, 
        t(S[[j]]) %*% X[[j]][s,]
        ) # dim: K(p-K)*R
    }
  }
  
  for (i in 1:N) {
    A13[[i]] = list()
    for (s in 1:Ti) {
      A13[[i]][[s]] = kronecker(
        L[i,], 
        (t(X[[i]][s,]) - t(hat.F[s,]) %*% t(A7[[i]]) - t(v[s,]) %*% Sigma.v.inv %*% t(A8[[i]])) %*% Theta %*% solve(A1.sqr[[i]]/Ti) %*% t(Theta)
        ) # dim: R*p
    }
  }
    
  A14 = matrix(0, K*(p-K), Ti*max(R,1))
  for (j in 1:N) {
    A14 = A14 + (1/N) * kronecker(
      B[j,] %*% t(L[j,]) %*% Sigma.L.inv, 
      t(S[[j]]) %*% t(X[[j]])
      ) # dim: K(p-K)*TR  
  }
  A15 = list() # list (of length N) of list (of length Ti)
  for (i in 1:N) {
    A15[[i]] = list()
    for (t in 1:Ti) {
      A15[[i]][[t]] = (- (hat.F %*% hat.F[t,] + v %*% Sigma.v.inv %*% v[t,]) %*% t(X[[i]][t,]) 
                       - (X[[i]] - hat.F %*% t(A7[[i]]) - v %*% Sigma.v.inv %*% t(A8[[i]])) %*% Theta %*% 
                         solve(A1.sqr[[i]]/ Ti) %*% t(Theta) %*% (X[[i]][t,] %*% t(X[[i]][t,]) - X[[i]][t,] %*% t(hat.F[t,]) %*% t(A7[[i]]) - A7[[i]] %*% hat.F[t,] %*% t(X[[i]][t,])) 
      ) %*% Theta %*% solve(A1.sqr[[i]]/Ti) %*% t(Theta) 
      A15[[i]][[t]] = kronecker(
        L[i,], 
        A15[[i]][[t]]
        ) # dim: TR*p
    }
  }
  
  for (i in 1:N) {
    for (t in 1:(Ti-1)) {
      for (s in (t+1):min(t+m, Ti)) {
        eta3 = eta3 + 
          1/(N*Ti) * A12[[s]] %*% A13[[i]][[s]] %*% X[[i]][t,] * E[t,i] + 
          1/(N*Ti^2) * A14 %*% A15[[i]][[s]] %*% X[[i]][t,] * E[t,i]
        eta4 = eta4 - 
          1/(N*Ti) * A12[[s]] %*% A13[[i]][[s]] %*% A7[[i]] %*% hat.F[t,] * E[t,i] -
          1/(N*Ti^2) * A14 %*% A15[[i]][[s]] %*% A7[[i]] %*% hat.F[t,] * E[t,i] - 
          1/(N*Ti) * A14 %*% kronecker(
            L[i,], 
            v.proj.perp %*% A1[[i]] %*% solve(A1.sqr[[i]]) %*% t(Theta) %*% (X[[i]][s,] %*% t(hat.F[s,])) %*% hat.F[t,] * E[t,i]
            )
      }
    }
  }
  
  
  eta5 = rep(0, Ti*max(R,1))
  eta6 = rep(0, Ti*max(R,1))
  
  A16 = matrix(0, (p-K)*K, Ti*max(R,1))
  for (j in 1:N) {
    A16 = A16 + 1/(N*sqrt(Ti)) * kronecker(
      B[j,]%*%t(L[j,]), 
      t(S[[j]])%*%t(X[[j]])
      )
  }
  A16 = A16 %*% beth.inv %*% kronecker(
    diag(1,max(R,1)), 
    v.proj.perp %*% F.proj.perp
    )
  
  A17 = list() # list of length Ti \sum_{j=1}^N ...
  A18 = list() # list (of length N) of list (of length Ti) 
  for (s in 1:Ti) {
    A17[[s]] = matrix(0, Ti*max(R,1), max(R,1))
    for (k in 1:N) {
      A17[[s]] = A17[[s]] + 
        (sqrt(Ti)/N) * kronecker(
          L[k,]%*%t(L[k,])%*%Sigma.L.inv, 
          A1[[k]]%*%solve(A1.sqr[[k]])%*%t(Theta)%*%X[[k]][s,]
          )
    }
  } # A13[[i]][[s]] can still be used here, combined with A17
  
  A18 = matrix(0, Ti*max(R,1), Ti*max(R,1))
  for (k in 1:N) {
    A18 = A18 + (sqrt(Ti)/N) * kronecker(
      L[k,]%*%t(L[k,])%*%Sigma.L.inv, 
      A1[[k]]%*%solve(A1.sqr[[k]])%*%t(Theta)%*%t(X[[k]])
      ) # dim: TR*TR
  } # A15[[i]][[s]] can still be used here, combined with A18
  
  for (i in 1:N) {
    for (t in 1:(Ti-1)) {
      for (s in (t+1):min(t+m,Ti)) {
        eta5 = eta5 + 
          1/(N*Ti) * A17[[s]] %*% A13[[i]][[s]] %*% X[[i]][t,] * E[t,i] + 
          1/(N*Ti^2) * A18 %*% A15[[i]][[s]] %*% X[[i]][t,] * E[t,i]
        eta6 = eta6 - 
          1/(N*Ti) * A17[[s]] %*% A13[[i]][[s]] %*% A7[[i]] %*% hat.F[t,] * E[t,i] - 
          1/(N*Ti^2) * A18 %*% A15[[i]][[s]] %*% A7[[i]] %*% hat.F[t,] * E[t,i] - 
          1/(N*Ti) * A18 %*% kronecker(L[i,], v.proj.perp %*% A1[[i]] %*% solve(A1.sqr[[i]]) %*% 
                                             t(Theta) %*% (X[[i]][s,] %*% t(hat.F[s,])) %*% hat.F[t,] * E[t,i])
      }
    }
  }
  eta5 = A16 %*% eta5
  eta6 = A16 %*% eta6
  
  bias.eta1 = A0 %*% solve(D0) %*% eta1
  bias.eta2 = A0 %*% solve(D0) %*% eta2
  bias.eta3 = A0 %*% solve(D0) %*% eta3
  bias.eta4 = A0 %*% solve(D0) %*% eta4
  bias.eta5 = A0 %*% solve(D0) %*% eta5
  bias.eta6 = A0 %*% solve(D0) %*% eta6
  
  cat("estimated zeta: ", round(bias.zeta1 + bias.zeta2, 5), "\n")
  cat("estimated xi:", round(bias.xi1 + bias.xi2 + bias.xi3 + bias.xi4, 5), "\n")
  cat("estimated eta:", round(bias.eta1 + bias.eta2 + bias.eta3 + bias.eta4 + bias.eta5 + bias.eta6, 5), "\n")
  ###
  BC.Theta.proj.vec = as.vector(Theta.proj) - 
    Ti^(-1)*(bias.zeta1 + bias.zeta2) - 
    N^(-1)*(bias.xi1 + bias.xi2 + bias.xi3 + bias.xi4) - 
    Ti^(-1)*(bias.eta1 + bias.eta2 + bias.eta3 + bias.eta4 + bias.eta5 + bias.eta6)
  
  return(BC.Theta.proj.vec)
}




commutation_matrix <- function(m, n) {
  K <- matrix(0, m*n, m*n)
  for (i in 1:m) {
    for (j in 1:n) {
      row <- (j - 1) * m + i      # vec(A^T)
      col <- (i - 1) * n + j      # vec(A)
      K[row, col] <- 1
    }
  }
  K
}


