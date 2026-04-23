setwd("codes/simu/estimation error and MSC in Section 4")

library("foreach"); library("doParallel"); cl <- makeCluster(36); registerDoParallel(cl) 

source("choose_KR.R")
source("RPC.R")

p = 10
# N = 100
# Ti = 100
NTi.list = list(c(20,200), c(50,200), c(100,200), c(200,200), c(200,20), c(200,50), c(200,100))
# NTi.list = list(c(20,200), c(50,200), c(200,20), c(200,50))
# K = 2  # reduced-rank
# R = 3  # latent factors
KR.list = list(c(2,2), c(0,2), c(2,0), c(0,0))

error.type = 1
# error.type = 2
# error.type = 3
# 1: IID normal; 2: heavy tailed (t-distribution); 3: heteroskedastic; 
# 4: serial correlation; 5: cross-sectional correlation 
# no need for 4 and 5 under the dynamic setting as we allow for the dependence between X and E

N.simu = 100
N.simu = 3

K.max = 4
R.max = 4
# K.max = 3
# R.max = 3

OUTALL = c()

for (i.KR in 1:4) 
{ 

  KR = KR.list[[i.KR]]; K0 = KR[1]; R0 = KR[2]
  # for (i.NTi in 1:4)
  for (i.NTi in 1:7)
  {
    NTi = NTi.list[[i.NTi]]; N = NTi[1]; Ti = NTi[2]
    
    OUT <- foreach(iter=1:N.simu, .combine='rbind') %do%
    {
      set.seed(iter)
      
      index = c()
      X = c()
      y = c()
      e = c()
      
      noise = 0.5
#     noise =  sqrt( (K0+1)/(R0+1)) / sqrt(2)
      
      # generating the data (Xi, yi)
      if (error.type == 1) {
        E = matrix(rnorm((Ti+50)*N), Ti+50, N) 
      } else if (error.type == 2) {
        E = matrix(rt((Ti+50)*N, df=10), Ti+50, N) * sqrt((10-2)/10)
      } else if (error.type == 3) {
        E = matrix(rnorm((Ti+50)*N), Ti+50, N)
        E[(row(E) + col(E)) %% 2 == 0] = E[(row(E) + col(E)) %% 2 == 0] * sqrt(4/3)
        E[(row(E) + col(E)) %% 2 == 1] = E[(row(E) + col(E)) %% 2 == 1] * sqrt(2/3)
      }
      
      rho.F = 0.5
      gamma.F = sqrt(1-rho.F^2)
      if (R0 > 0) {
        sn = sqrt(R0)
        F1 = c()
        F1 = rbind(F1, rnorm(R0))
        for (t in 2:(Ti+50)) {
          F1 = rbind(F1, rho.F*F1[t-1,]+ gamma.F*rnorm(R0))
        }
        # F1 = matrix(rnorm((Ti+50)*R0), Ti+50, R0)
        F0 = as.matrix(F1[51:(Ti+50),])
        F0.proj.perp = diag(1, Ti) - F0 %*% solve(t(F0)%*%F0) %*% t(F0)
        L = sqrt(12) * matrix(runif(N*R0)-0.5, N, R0)
        # L = matrix(rnorm(N*R0), N, R0)
      } else {
        sn = 1
        F1 = matrix(0, Ti+50, 1)
        F0 = matrix(0, Ti, 1)
        F0.proj.perp = diag(1, Ti)
        L = matrix(0, N, 1)
      }
      
      rho.v = 0.6
      gamma.v = sqrt(1-rho.v^2)
      v1 = rnorm(1)
      for (t in 2:(Ti+50)) {
        v1 = c(v1, rho.v*v1[t-1] + gamma.v*rnorm(1))
      }
      Z = F0.proj.perp %*% v1[51:(Ti+50)]
      v = Z * norm(v1[51:(Ti+50)], "2") / norm(Z, "2")
      v1[51:(Ti+50)] = v
      v = as.matrix(v)
      w = runif(N) + 0.5
      
      rho.y = 0.5
      rho.X = 0.6
      gamma.X = sqrt(1-rho.X^2)
      if (K0 > 0) {
        Theta1 = svd(matrix(rnorm((p-1) * K0), p-1, K0))$u
        Theta = rbind(c(rho.y, rep(0, K0-1)), Theta1)
      } else {
        Theta = rep(0, p)
      }
      B = c()
      for (i in 1:N) { 
        index = c(index, rep(i, Ti))
        Xi = c()
        yi = c()
        yi[1] = 0
        h = rnorm(1)
        for (j in 2:(p-2)) {
          h = c(h, rho.X*h[j-1] + gamma.X*rnorm(1)) 
        }
        Xi = rbind(Xi, c(0, h, v1[1]*w[i]))
        if (K0 > 0) {
          bi = runif(K0) + 0.5
        } else {
          bi = 0
        }
        B = rbind(B, bi)
        for (t in 2:(Ti+50)) {
          # Xi = rbind(Xi, c(yi[t-1], rho.X*Xi[t-1,2:(p-1)] + gamma.X*rnorm(p-2), v1[t]*w[i]))
          h = rnorm(1)
          for (j in 2:(p-2)) {
            h = c(h, rho.X*h[j-1] + gamma.X*rnorm(1)) 
          }
          Xi = rbind(Xi, c(yi[t-1], h, v1[t]*w[i]))
          yi = c(yi, t(Xi[t,]) %*% Theta %*% bi + noise * (t(F1[t,]) %*% L[i,] + sn * E[t,i]))
        }
        
        y = c(y, yi[51:(Ti+50)])
        X = rbind(X, Xi[51:(Ti+50),])
      }
      
      # fit = choose.KR(X, y, v, index1 = index, 
      #                 K.max = K.max, R.max = R.max, 
      #                 cc.K = 1, cc.R = 1.6)
      fit = choose.KR(X, y, v, index1 = index,
                      K.max = K.max, R.max = R.max,
                      cc.K = 1, cc.R = 2)
      
      cancor.Theta.RPC = 0
      cancor.F.RPC = 0
      err.Theta.RPC = 0
      err.F.RPC = 0
      
      fit.oracle = fit$fit.collection[[K0*(K.max+1)+R0+1]]
      if (K0>0) {
        hat.Theta = fit.oracle$hat.Theta
        cancor.Theta.RPC = min(cancor(X%*%Theta, X%*%hat.Theta)$cor)
        
        Theta.proj = Theta %*% solve(t(Theta)%*%Theta) %*% t(Theta)
        hat.Theta.proj = hat.Theta %*% solve(t(hat.Theta)%*%hat.Theta) %*% t(hat.Theta)
        err.Theta.RPC = norm(hat.Theta.proj-Theta.proj, "f")
      }
      if (R0>0) {
        hat.F = fit.oracle$hat.F
        cancor.F.RPC = min(cancor(F0, hat.F)$cor)
        
        F.proj = F0 %*% solve(t(F0)%*%F0) %*% t(F0)
        hat.F.proj = hat.F %*% solve(t(hat.F)%*%hat.F) %*% t(hat.F)
        err.F.RPC = norm(hat.F.proj-F.proj, "f")
      }
      
      K.ICp1 = fit$K.ICp1; R.ICp1 = fit$R.ICp1
      K.PCp1 = fit$K.PCp1; R.PCp1 = fit$R.PCp1
      
      K.ICp2 = fit$K.ICp2; R.ICp2 = fit$R.ICp2
      K.PCp2 = fit$K.PCp2; R.PCp2 = fit$R.PCp2
      
      K.ICp3 = fit$K.ICp3; R.ICp3 = fit$R.ICp3
      K.PCp3 = fit$K.PCp3; R.PCp3 = fit$R.PCp3
      
      # K.BIC1 = fit$K.BIC1; R.BIC1 = fit$R.BIC1
      # K.AIC1 = fit$K.AIC1; R.AIC1 = fit$R.AIC1
      # 
      # K.BIC2 = fit$K.BIC2; R.BIC2 = fit$R.BIC2
      # K.AIC2 = fit$K.AIC2; R.AIC2 = fit$R.AIC2
      # 
      # K.BIC3 = fit$K.BIC3; R.BIC3 = fit$R.BIC3
      # K.AIC3 = fit$K.AIC3; R.AIC3 = fit$R.AIC3
      
      
      out = c(cancor.Theta.RPC, cancor.F.RPC, 
              err.Theta.RPC, err.F.RPC,
              K.PCp1, R.PCp1, K.ICp1, R.ICp1,  
              K.PCp2, R.PCp2, K.ICp2, R.ICp2, 
              K.PCp3, R.PCp3, K.ICp3, R.ICp3 
              # K.AIC1, R.AIC1, K.BIC1, R.BIC1,
              # K.AIC2, R.AIC2, K.BIC2, R.BIC2,
              # K.AIC3, R.AIC3, K.BIC3, R.BIC3
              )
      out
      
    }
    
    names = c("cancor.Theta.RPC", "cancor.F.RPC", 
              "err.Theta.RPC", "err.F.RPC",
              "K.PCp1", "R.PCp1", "K.ICp1", "R.ICp1", 
              "K.PCp2", "R.PCp2", "K.ICp2", "R.ICp2", 
              "K.PCp3", "R.PCp3", "K.ICp3", "R.ICp3" 
              # "K.AIC1", "R.AIC1", "K.BIC1", "R.BIC1",
              # "K.AIC2", "R.AIC2", "K.BIC2", "R.BIC2",
              # "K.AIC3", "R.AIC3", "K.BIC3", "R.BIC3"
              )
    colnames(OUT) = names
   
    
#    cat("result for (N=", N, ", T=", Ti, ", K=", K0, ", R=", R0, "):", "\n", sep="")
    # canonical correlation between: hat.Theta and Theta0; hat.F and F0
    # estimation error: distance between the projection of hat.Theta and Theta0; hat.F and F0, 
    # model selection: frequency of K (from 0 to K.max) and R (from 0 to R.max)
    
    out = matrix(0, 6, K.max + R.max + 10)
    # out = matrix(0, 12, K.max + R.max + 10)
    # rows: 6 different criteria
    # columns: (K, R, N, Ti, cancor.Theta, cancor.F, err.Theta, err.F, 0, 1, ..., K.max, 0, 1, ..., R.max)
    cnames = c("K0", "R0", "N", "Ti", 
               "cancor.Theta.RPC", "cancor.F.RPC", "err.Theta.RPC", "err.F.RPC",
               paste0("K=", 0:K.max), paste0("R=", 0:R.max))
    colnames(out) = cnames
      
    rnames = c()

    for (i in 1:6) 
    { 
      j = (i-1)*2 + 5

      A = as.numeric(table(factor(OUT[,j], levels=0:K.max)))
      B = as.numeric(table(factor(OUT[,j+1], levels=0:R.max)))
      kr = c(A, B)     
     
      out[i, ] = c(K0, R0, N, Ti, colMeans(OUT[,1:4]), kr)
      rnames = c(rnames, substring(colnames(OUT)[j],3)) 
    }
    
    rownames(out) = rnames
    
    print(out)
    
    OUTALL = rbind(OUTALL, out)
    
    write.csv(OUTALL, file = "RPCmsc_OUT.csv")
  }
}


