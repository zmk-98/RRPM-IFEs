setwd("codes/real data analysis")

library(pwt10) 

source("RPC.R")
source("choose_KR.R")
# source("NNR.R")
#############

# MSC (model selection criterion) approach

B0 = as.matrix(pwt10.01)

cols = c("country", "isocode", "rgdpe", 
         "pop", "year", 
         "rnna", "emp", "hc", "pl_c", "pl_i", "csh_i", "csh_g")
B = B0[,cols]
B = B0[as.numeric(B0[,"year"])>=1968, ]

nac = c()
for (i in 1:ncol(B))
  nac = c(nac, sum(is.na(B[,i])))

B = B[,nac<=2411]
isocode = unique(B[,2])
nac1 = c()
C = c()
for (i in 1:length(isocode))
{
  if (sum(is.na(B[which(B[,2]==isocode[i]),]))==0)
    C = rbind(C, B[which(B[,2]==isocode[i]),])
}

C = as.data.frame(C)             
y1 = log(as.numeric(C$rgdpe)/as.numeric(C$pop))

isocode = unique(C$isocode)

nyrs = nrow(C)


# consider N = 91, T = 50, p = 11
# 
# response:
#   log(rgdpe/pop)
# 
# 11 predictors:
#   lag 1 of log(rgdpe/pop), lag 2 of log(rgdpe/pop)
#   log(rnna/emp), lag 1 of log(rnna/emp)
#   log(hc), log(emp/pop)
#   log(pl_i/pl_c), lag 1 of log(pl_i/pl_c)
#   csh_i, csh_g, year
lag0 = (as.numeric(C$year)>1969)
lag1 = (as.numeric(C$year)>1968 & as.numeric(C$year)<2019)
lag2 = (as.numeric(C$year)>1967 & as.numeric(C$year)<2018)

Z1 = cbind(y1, C$rnna, C$emp, C$hc, C$pop, C$pl_i, C$pl_c, C$csh_i, C$csh_g, C$year)
Z1 = matrix(as.numeric(Z1), nrow=length(y1))

X1 = cbind(y1, # lag 0, lag 1, lag 2
           log(Z1[,2]/Z1[,3]), # log(rnna/emp) & lag 1
           log(Z1[,4]), # log(hc) 
           log(Z1[,3]/Z1[,5]), # log(emp/pop) 
           log(Z1[,6]/Z1[,7]), # log(pl_i/pl_c) & lag 1 
           Z1[,8], # csh_i
           Z1[,9], # csh_g
           Z1[,10] # year
)

library(bestNormalize)
for (i in 1:ncol(X1))
{
  bn_boxcox_obj <- bestNormalize::boxcox(X1[,i]-min(X1[,i])+sd(X1[,i])/10)
  X1[,i] <- bn_boxcox_obj$x.t
}
index1 = C$isocode[lag0]

X1 = scale(X1)

y = X1[lag0, 1] # log(rgdpe/pop)
X = cbind(X1[lag1,1], X1[lag2,1], # lag 1 & lag 2 of log(rgdpe/pop)
          X1[lag0,2], X1[lag1,2], # log(rnna/emp) & lag 1
          X1[lag0,3], # log(hc) 
          X1[lag0,4], # log(emp/pop) 
          X1[lag0,5], X1[lag1,5], # log(pl_i/pl_c) & lag 1
          X1[lag0,6], # csh_i
          X1[lag0,7], # csh_g
          X1[lag0,8]  # year
) 
# no pl_g, no csh_x+csh_m, p=10

n = dim(X)[1]
ind1 = C$year[lag0]



###

XX = X
YY = y

ind = unique(index1)
N = length(ind)      # N = 91
Ti = length(y)/N     # Ti = 50
p = ncol(X)          # p = 11
N; Ti; p


# The codes below shows that only the first predictor is low-rank (of rank 1)
# pred = list()
# for (j in 1:p) {
#   pred[[j]] = X[1:Ti,j]
#   for (i in 2:N) {
#     pred[[j]] = cbind(pred[[j]], X[((i-1)*Ti+1):(i*Ti),j])
#   }
#   print(qr(pred[[j]])$rank)
# } 
v = as.matrix(X[1:Ti,11] / norm(X[1:Ti,11], "2")) # if p = 11


# X is a matrix of size NT x p
# y is a vector of length NT

set.seed(123)

# --- Model selection criterion ---
MSC_result2 = choose.KR(X = X, y = y, v = v, index1 = index1,
                        K.max = 8, R.max = 8, cc.K = 1, cc.R = 2)

c(MSC_result2$K.ICp1, MSC_result2$R.ICp1, 
  MSC_result2$K.ICp2, MSC_result2$R.ICp2, 
  MSC_result2$K.ICp3, MSC_result2$R.ICp3) 

c(MSC_result2$K.PCp1, MSC_result2$R.PCp1, 
  MSC_result2$K.PCp2, MSC_result2$R.PCp2, 
  MSC_result2$K.PCp3, MSC_result2$R.PCp3) 



# Interpretation using RPC with (K=2, R=1)
Theta21 = MSC_result2$fit.collection[[2*(8+1)+1+1]]$hat.Theta
B21 = MSC_result2$fit.collection[[2*(8+1)+1+1]]$hat.B
UTheta21 = svd(Theta21%*%t(B21))$u[,1:2]
round(UTheta21, 3)
# [,1]   [,2]
# [1,] -0.989  0.003
# [2,]  0.143 -0.003
# [3,] -0.019 -0.479
# [4,] -0.025  0.637
# [5,] -0.024 -0.519
# [6,] -0.016  0.067
# [7,] -0.003  0.027
# [8,] -0.001  0.004
# [9,] -0.010  0.047
# [10,]  0.008  0.018
# [11,] -0.002  0.296
round(svd(Theta21%*%t(B21))$d[1:2], 4)
# [1] 10.2824  3.8701
round(svd(Theta21%*%t(B21))$d[1:2]/sqrt(N), 4)
# [1] 1.0779 0.4057

F21 = MSC_result2$fit.collection[[2*(8+1)+1+1]]$hat.F
L21 = MSC_result2$fit.collection[[2*(8+1)+1+1]]$hat.L
UF21 = svd(F21%*%t(L21))$u[,1]
round(svd(F21%*%t(L21))$d[1]/sqrt(N*Ti), 4)
# [1] 0.0618

SST = t(y) %*% y
SSR.theta1 = 0
SSR.theta2 = 0
SSR.XF1 = 0
SSR.all = 0

for (i in 1:N) {
  Xi = X[((i-1)*Ti+1):(i*Ti), ]
  yi = y[((i-1)*Ti+1):(i*Ti)]
  Xi.Theta = Xi %*% UTheta21
  
  SSR.theta1 = SSR.theta1 + t(yi) %*% (diag(1,Ti) - Xi.Theta[,1] %*% solve(t(Xi.Theta[,1])%*%Xi.Theta[,1]) %*% t(Xi.Theta[,1])) %*% yi
  SSR.theta2 = SSR.theta2 + t(yi) %*% (diag(1,Ti) - Xi.Theta %*% solve(t(Xi.Theta)%*%Xi.Theta) %*% t(Xi.Theta)) %*% yi
  
  XF.bind = cbind(Xi.Theta, UF21)
  SSR.XF1 = SSR.XF1 + t(yi) %*% (diag(1,Ti) - XF.bind %*% solve(t(XF.bind)%*%XF.bind) %*% t(XF.bind)) %*% yi
  
  SSR.all = SSR.all + t(yi) %*% (diag(1,Ti) - Xi %*% solve(t(Xi)%*%Xi) %*% t(Xi)) %*% yi
}

round(c(SST, SSR.theta1, SSR.theta2, SSR.XF1, SSR.all), 2)

R.sqr.theta1 = 1 - SSR.theta1/SST
R.sqr.theta2 = 1 - SSR.theta2/SST
R.sqr.XF1 = 1 - SSR.XF1/SST
R.sqr.all = 1 - SSR.all/SST

round(c(R.sqr.theta1, R.sqr.theta2, R.sqr.XF1, R.sqr.all), 4)
# [1] 0.9966 0.9970 0.9980 0.9984
# coefficients of hat Theta and R square under other settings can be obtained by the same method


iso_keep = unique(C$isocode)
years = 1970:2019
year_num = as.numeric(B0[,"year"])
iso_num = B0[,"isocode"]

D = B0[year_num >= 1970 & year_num <= 2019 & iso_num %in% iso_keep, ]


country_cor = function(B0, iso_keep, varname, f_hat, 
                       year_start = 1970, year_end = 2019, min_obs = 10) {
  year_num = as.numeric(B0[,"year"])
  iso_all = B0[,"isocode"]
  
  D = B0[year_num >= year_start & year_num <= year_end & iso_all %in% iso_keep, ]
  years = year_start:year_end
  
  out = NULL
  
  for (iso in iso_keep) {
    idx_iso = D[,"isocode"] == iso
    Di = D[idx_iso, ]
    if (nrow(Di) == 0) next
    
    yi = as.numeric(Di[,"year"])
    zi = suppressWarnings(as.numeric(Di[, varname]))
    
    ord = order(yi)
    yi = yi[ord]
    zi = zi[ord]
    
    idx = !is.na(zi)
    if (sum(idx) < min_obs) next
    
    f_i = f_hat[match(yi[idx], years)]
    z_i = zi[idx]
    
    ok = !is.na(f_i) & !is.na(z_i)
    Ti = sum(ok)
    if (Ti < min_obs) next
    
    ri = cor(f_i[ok], z_i[ok])
    
    out = rbind(out, data.frame(isocode = iso, T_i = Ti, cor = ri))
  }
  
  out
}

res_avh = country_cor(B0, iso_keep, "avh", F21, min_obs = 30)
mean(abs(res_avh$cor))
# [1] 0.184168
weighted.mean(abs(res_avh$cor), w = res_avh$T_i)
# [1] 0.1803419
nrow(res_avh)
# [1] 49
sum(abs(res_avh$cor)>0.3)
# [1] 11

res_labsh = country_cor(B0, iso_keep, "labsh", F21, min_obs = 30)
mean(abs(res_labsh$cor))
# [1] 0.1820747
weighted.mean(abs(res_labsh$cor), w = res_labsh$T_i)
# [1] 0.1820747
nrow(res_labsh)
# [1] 78
sum(abs(res_labsh$cor)>0.3)
# [1] 12

res_ctfp = country_cor(B0, iso_keep, "ctfp", F21, min_obs = 30)[-74,]
# USA is always set as 1, exclude this NA value
mean(abs(res_ctfp$cor))
# [1] 0.2364581
weighted.mean(abs(res_ctfp$cor), w = res_ctfp$T_i)
# [1] 0.2364581
nrow(res_ctfp)
# [1] 77
sum(abs(res_ctfp$cor)>0.3)
# [1] 23


#####
country_joint_reg = function(B0, iso_keep, f_hat,
                             year_start = 1970, year_end = 2019,
                             min_obs = 30) {
  
  year_num = as.numeric(B0[, "year"])
  iso_all  = B0[, "isocode"]
  
  D = B0[year_num >= year_start & year_num <= year_end & iso_all %in% iso_keep, ]
  years = year_start:year_end
  
  out = NULL
  for (iso in iso_keep[-87]) {
    # exclude "USA", because ctfp for USA = 1 is the benchmark
    idx_iso = D[, "isocode"] == iso
    Di = D[idx_iso, ]
    if (nrow(Di) == 0) next
    
    yi    = as.numeric(Di[, "year"])
    avh   = suppressWarnings(as.numeric(Di[, "avh"]))
    labsh = suppressWarnings(as.numeric(Di[, "labsh"]))
    ctfp  = suppressWarnings(as.numeric(Di[, "ctfp"]))
    
    ord = order(yi)
    yi    = yi[ord]
    avh   = avh[ord]
    labsh = labsh[ord]
    ctfp  = ctfp[ord]
    
    f_i = f_hat[match(yi, years)]
    
    ok = !is.na(f_i) & !is.na(avh) & !is.na(labsh) & !is.na(ctfp)
    Ti = sum(ok)
    if (Ti < min_obs) next
    
    dat_i = data.frame(
      F = f_i[ok],
      avh = avh[ok],
      labsh = labsh[ok],
      ctfp = ctfp[ok]
    )
    
    fit = lm(F ~ avh + labsh + ctfp, data = dat_i)
    sm = summary(fit)
    cf = sm$coefficients
    
    out = rbind(out, data.frame(
      isocode = iso,
      T_i = Ti,
      beta_avh = cf["avh", "Estimate"],
      t_avh    = cf["avh", "t value"],
      beta_labsh = cf["labsh", "Estimate"],
      t_labsh    = cf["labsh", "t value"],
      beta_ctfp = cf["ctfp", "Estimate"],
      t_ctfp    = cf["ctfp", "t value"],
      r2 = sm$r.squared
    ))
  }
  
  out
}

res_joint = country_joint_reg(B0, iso_keep, F21, min_obs = 30)
nrow(res_joint)
# [1] 45
summary(res_joint$r2)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01737 0.11136 0.17998 0.21447 0.31181 0.62884 

mean(abs(res_joint$t_avh))
# [1] 1.69494
mean(abs(res_joint$t_labsh))
# [1] 1.508174
mean(abs(res_joint$t_ctfp))
# [1] 2.003057

mean(abs(res_joint$t_avh) > 2, na.rm = TRUE)
# [1] 0.3555556
mean(abs(res_joint$t_labsh) > 2, na.rm = TRUE)
# [1] 0.3555556
mean(abs(res_joint$t_ctfp) > 2, na.rm = TRUE)
# [1] 0.4222222
