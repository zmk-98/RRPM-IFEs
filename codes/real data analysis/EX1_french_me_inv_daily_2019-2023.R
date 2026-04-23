#' @case 2019.1.1-2023.12.31
#' 
#' 
#' @description Five factors:
#' MKT-RF: market excess return (market premium), 
#' SMB: small-minus-big (size factor), 
#' HML: high-minus-low (value factor),
#' RMW: robust-minus-weak (profitability factor),
#' CMA: conservative-minus-aggressive (investment factor)
#' @description use c.K = 1 and c.R = 2 or 1.6
#' 

setwd("codes/real data analysis")

# RRPM-IFEs
source("RPC.R")
source("choose_KR.R")

Factors = read.csv("5factors_daily.csv")
Portfolios = read.csv("100portfolios_me_inv_daily.csv") # we use value-weighted return

X0 = as.matrix(Factors[, 2:6])
Y = Portfolios[,2:101]
Y = Y*(Y> -99) + 0

dates = as.Date(as.character(Factors[,1]), "%Y%m%d")
idx_2019_first = which(dates >= as.Date("2019-01-01"))[1]
idx_2023_last = tail(which(dates <= as.Date("2023-12-31")), 1)


X0 = X0[idx_2019_first:idx_2023_last, ]
Y = Y[idx_2019_first:idx_2023_last, ]
X0 = scale(X0)

ind = c()
X = c()
for (i in 1:100)
{
  ind = c(ind, rep(i, length(X0[,1])))
  X = rbind(X, X0)
}
y = c(as.matrix(Y))
index1 = ind


XX = X
YY = y

ind = unique(index1)
N = length(ind)     # N = 100
Ti = length(y)/N    # Ti = 1258
p = ncol(X)         # p = 5

Xall = X
yall = y
X = list()
y = list()
I = list()
for (i in 1:N)
{
  I = which(index1==ind[i])
  X[[i]] = XX[I,]
  y[[i]] = YY[I]
  
  X[[i]] = X[[i]] - matrix(colMeans(X[[i]]), Ti, ncol(X[[i]]), byrow=TRUE)
  
  
  y[[i]] = y[[i]]-mean(y[[i]])
  
  Xall[I,] = X[[i]]
  yall[I] = y[[i]]
}
X_list = X
y_list = y
X = Xall
y = yall

v = X_list[[1]]
for (j in 1:5) {
  v[,j] = v[,j] / norm(v[,j], "2")
}


# X is a matrix of size NT x p
# y is a vector of length NT
# v is a matrix of size T x p (all predictors are low-rank)
# X_list is a list of length N, each matrix of size T x p 
# y_list is a list of length N, each vector of length T

set.seed(123)

# N = 100, T = 1258 from 2019 to 2023
# --- Model selection criterion ---
MSC_result1 <- choose.KR(X = X, y = y, v = v, index1 = index1, 
                         K.max = 5, R.max = 5, cc.K = 1, cc.R = 2)

c(MSC_result1$K.ICp1, MSC_result1$R.ICp1, # (4, 0)
  MSC_result1$K.ICp2, MSC_result1$R.ICp2, # (4, 0)
  MSC_result1$K.ICp3, MSC_result1$R.ICp3) # (4, 0)
# [1] 4 0 4 0 4 0

c(MSC_result1$K.PCp1, MSC_result1$R.PCp1, # (4, 1)
  MSC_result1$K.PCp2, MSC_result1$R.PCp2, # (4, 1)
  MSC_result1$K.PCp3, MSC_result1$R.PCp3) # (4, 1)
# [1] 4 1 4 1 4 1


Theta41 = MSC_result1$fit.collection[[4*(5+1)+1+1]]$hat.Theta
B41 = MSC_result1$fit.collection[[4*(5+1)+1+1]]$hat.B
UTheta41 = svd(Theta41%*%t(B41))$u[,1:4]
round(UTheta41, 3)
# [,1]   [,2]   [,3]   [,4]
# [1,] -0.931  0.346  0.081 -0.062
# [2,] -0.342 -0.928 -0.074 -0.013
# [3,] -0.119  0.018 -0.696  0.458
# [4,]  0.048  0.139 -0.624 -0.037
# [5,] -0.007  0.007  0.338  0.886
round(svd(Theta41%*%t(B41))$d[1:4], 3)
# [1] 14.565  3.075  2.743  2.023
round(svd(Theta41%*%t(B41))$d[1:4]/sqrt(N), 3)
# [1] 1.457 0.307 0.274 0.202

F41 = MSC_result1$fit.collection[[4*(5+1)+1+1]]$hat.F
L41 = MSC_result1$fit.collection[[4*(5+1)+1+1]]$hat.L
round(svd(F41%*%t(L41))$d[1]/sqrt(N*Ti), 3)
# [1] 0.19


SST = t(y-mean(y)) %*% (y-mean(y))
SSR1 = 0
SSR2 = 0
SSR3 = 0
SSR4 = 0

SSR.XF = 0
SSR.3f = 0
SSR.5f = 0

for (i in 1:N) {
  Xi = X[((i-1)*Ti+1):(i*Ti), ]
  yi = y[((i-1)*Ti+1):(i*Ti)]
  Xi.Theta = Xi %*% UTheta41
  
  SSR1 = SSR1 + t(yi) %*% (diag(1,Ti) - Xi.Theta[,1] %*% solve(t(Xi.Theta[,1])%*%Xi.Theta[,1]) %*% t(Xi.Theta[,1])) %*% yi
  SSR2 = SSR2 + t(yi) %*% (diag(1,Ti) - Xi.Theta[,1:2] %*% solve(t(Xi.Theta[,1:2])%*%Xi.Theta[,1:2]) %*% t(Xi.Theta[,1:2])) %*% yi
  SSR3 = SSR3 + t(yi) %*% (diag(1,Ti) - Xi.Theta[,1:3] %*% solve(t(Xi.Theta[,1:3])%*%Xi.Theta[,1:3]) %*% t(Xi.Theta[,1:3])) %*% yi
  SSR4 = SSR4 + t(yi) %*% (diag(1,Ti) - Xi.Theta %*% solve(t(Xi.Theta)%*%Xi.Theta) %*% t(Xi.Theta)) %*% yi
  
  XF = cbind(Xi.Theta, F41)
  SSR.XF = SSR.XF + t(yi) %*% (diag(1,Ti) - XF %*% solve(t(XF)%*%XF) %*% t(XF)) %*% yi
  
  SSR.3f = SSR.3f + t(yi) %*% (diag(1,Ti) - Xi[,1:3] %*% solve(t(Xi[,1:3])%*%Xi[,1:3]) %*% t(Xi[,1:3])) %*% yi
  SSR.5f = SSR.5f + t(yi) %*% (diag(1,Ti) - Xi %*% solve(t(Xi)%*%Xi) %*% t(Xi)) %*% yi
}

round(c(SST, SSR1, SSR2, SSR3, SSR4, SSR.XF, SSR.3f, SSR.5f), 2)
# [1] 391595.79  80083.99  69666.91  58604.51  53743.61  49179.31  59013.94  53231.55

R.sqr1 = 1 - SSR1/SST
R.sqr2 = 1 - SSR2/SST
R.sqr3 = 1 - SSR3/SST
R.sqr4 = 1 - SSR4/SST

R.sqr.XF = 1 - SSR.XF/SST

R.sqr.3f = 1 - SSR.3f/SST
R.sqr.5f = 1 - SSR.5f/SST

round(c(R.sqr1, R.sqr2, R.sqr3, R.sqr4, R.sqr.XF, R.sqr.3f, R.sqr.5f), 4)
# [1] 0.7955 0.8221 0.8503 0.8628 0.8744 0.8493 0.8641
# R.sqr1 - R.sqr4 for adding 4 linear combinations of 5 factors sequentially
# R.sqr.XF means using 4 linear combinations + latent factor
# R.sqr.3f and R.sqr.5f for using 3 factors and 5 factors, respectively

# calculate the contribution
cov.X = cov(Xi[1:Ti,])
numerator = sum(diag(cov.X %*% Theta41 %*% solve(t(Theta41)%*%cov.X%*%Theta41) %*% t(Theta41) %*% cov.X))
round(numerator/sum(diag(cov.X)), 4)
# [1] 0.937






## =========================================================
## Figure: latent factor time series + loadings-size relation
## =========================================================

## ---------- 1. Prepare data ----------
f <- F41[, 1]

size_curvature <- rowMeans(Y[,1:10]) + rowMeans(Y[,91:100]) - 2*rowMeans(Y[,41:60])
cor(f, size_curvature)
# [1] 0.5537905
summary(lm(f~size_curvature))
# Call:
#   lm(formula = f ~ size_curvature)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -5.4579 -0.4574 -0.0153  0.3849  6.5166 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    0.002837   0.023495   0.121    0.904    
# size_curvature 0.480949   0.020404  23.571   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.8333 on 1256 degrees of freedom
# Multiple R-squared:  0.3067,	Adjusted R-squared:  0.3061 
# F-statistic: 555.6 on 1 and 1256 DF,  p-value: < 2.2e-16



## quadratic fit for loadings vs size
loadings = L41[,1]
size = rep(1:10, each = 10)
fit_quad <- lm(loadings ~ size + I(size^2))
summary(fit_quad)
# Call:
#   lm(formula = loadings ~ size + I(size^2))
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.17039 -0.04811 -0.01404  0.03095  0.78854 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.502533   0.038559   13.03   <2e-16 ***
#   size        -0.211486   0.016104  -13.13   <2e-16 ***
#   I(size^2)    0.016005   0.001427   11.22   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1037 on 97 degrees of freedom
# Multiple R-squared:  0.6961,	Adjusted R-squared:  0.6898 
# F-statistic: 111.1 on 2 and 97 DF,  p-value: < 2.2e-16



pdf("size_loading_plot.pdf", width = 8, height = 3.5)

par(mar = c(4, 4, 1, 1))  # 减少上下空白

size_grid <- seq(min(size), max(size), length.out = 300)
pred_quad <- predict(fit_quad, newdata = data.frame(size = size_grid))

ylim_bot <- range(c(loadings, pred_quad))
ylim_bot <- ylim_bot + c(-0.04, 0.04) * diff(ylim_bot)

plot(size, loadings,
     pch = 1,
     cex = 0.7,
     ylim = ylim_bot,
     xlab = "Size", ylab = "Latent factor loadings")

lines(size_grid, pred_quad, lwd = 1.8)

dev.off()
























# approximate factor model
source("bai_ng_select.R")
source("eigen_ratio_select.R")
source("growth_ratio_select.R")
source("Bootstrap factors.R")

Y_demean <- scale(Y, center = TRUE, scale = FALSE)
Y_std <- scale(Y, center = TRUE, scale = TRUE)

### Kaiser's criterion
# scaling or not by section has no effect on Kaiser's selection
eigvals = eigen(cor(Y))$values
sum(eigvals > 1 + sqrt(N/(Ti-1)))
# [1] 4
round(eigvals, 3)[1:15]
# [1] 74.700  3.942  2.743  1.616  0.923  0.856  0.559  0.539  0.479  0.444  0.408  0.362  0.347  0.344
# [15]  0.327

eigvals_demean = eigen(cor(Y_demean))$values
sum(eigvals_demean > 1 + sqrt(N/(Ti-1)))
round(eigvals_demean, 3)[1:15]

eigvals_std = eigen(cor(Y_std))$values
sum(eigvals_std > 1 + sqrt(N/(Ti-1)))
round(eigvals_std, 3)[1:15]


### Bai-Ng's model selection criteria
bai_ng_1 <- bai_ng_select(Y, rmax = 15, cc.R = 1)
c(bai_ng_1$IC1, bai_ng_1$IC2, bai_ng_1$IC3)
# [1] 8 8 8
c(bai_ng_1$PC1, bai_ng_1$PC2, bai_ng_1$PC3)
# [1]  9  9 10

bai_ng_demean_1 <- bai_ng_select(Y_demean, rmax = 15, cc.R = 1)
c(bai_ng_demean_1$IC1, bai_ng_demean_1$IC2, bai_ng_demean_1$IC3)
# [1] 8 8 8
c(bai_ng_demean_1$PC1, bai_ng_demean_1$PC2, bai_ng_demean_1$PC3)
# [1]  9  9 10

bai_ng_std_1 <- bai_ng_select(Y_std, rmax = 15, cc.R = 1)
c(bai_ng_std_1$IC1, bai_ng_std_1$IC2, bai_ng_std_1$IC3)
# [1] 6 6 6
c(bai_ng_std_1$PC1, bai_ng_std_1$PC2, bai_ng_std_1$PC3)
# [1] 6 6 6

bai_ng_1.6 <- bai_ng_select(Y, rmax = 15, cc.R = 1.6)
c(bai_ng_1.6$IC1, bai_ng_1.6$IC2, bai_ng_1.6$IC3)
# [1] 4 4 4
c(bai_ng_1.6$PC1, bai_ng_1.6$PC2, bai_ng_1.6$PC3)
# [1] 4 4 5

bai_ng_demean_1.6 <- bai_ng_select(Y_demean, rmax = 15, cc.R = 1.6)
c(bai_ng_demean_1.6$IC1, bai_ng_demean_1.6$IC2, bai_ng_demean_1.6$IC3)
# [1] 4 4 4
c(bai_ng_demean_1.6$PC1, bai_ng_demean_1.6$PC2, bai_ng_demean_1.6$PC3)
# [1] 4 4 5

bai_ng_std_1.6 <- bai_ng_select(Y_std, rmax = 15, cc.R = 1.6)
c(bai_ng_std_1.6$IC1, bai_ng_std_1.6$IC2, bai_ng_std_1.6$IC3)
# [1] 4 4 4
c(bai_ng_std_1.6$PC1, bai_ng_std_1.6$PC2, bai_ng_std_1.6$PC3)
# [1] 4 4 4

bai_ng_2 <- bai_ng_select(Y, rmax = 15, cc.R = 2)
c(bai_ng_2$IC1, bai_ng_2$IC2, bai_ng_2$IC3)
# [1] 4 4 4
c(bai_ng_2$PC1, bai_ng_2$PC2, bai_ng_2$PC3)
# [1] 4 4 4

bai_ng_demean_2 <- bai_ng_select(Y_demean, rmax = 15, cc.R = 2)
c(bai_ng_demean_2$IC1, bai_ng_demean_2$IC2, bai_ng_demean_2$IC3)
# [1] 4 4 4
c(bai_ng_demean_2$PC1, bai_ng_demean_2$PC2, bai_ng_demean_2$PC3)
# [1] 4 4 4

bai_ng_std_2 <- bai_ng_select(Y_std, rmax = 15, cc.R = 2)
c(bai_ng_std_2$IC1, bai_ng_std_2$IC2, bai_ng_std_2$IC3)
# [1] 3 3 3
c(bai_ng_std_2$PC1, bai_ng_std_2$PC2, bai_ng_std_2$PC3)
# [1] 4 4 4


### eigen ratio
er <- eigen_ratio_select(Y, rmax = 15)
er$r_hat
# [1] 1
round(er$table, 3)
# r lambda_r lambda_r1  ratio
# 1   1    1.708     0.088 19.489
# 2   2    0.088     0.077  1.145
# 3   3    0.077     0.042  1.807
# 4   4    0.042     0.025  1.728
# 5   5    0.025     0.021  1.184
# 6   6    0.021     0.018  1.166
# 7   7    0.018     0.016  1.095
# 8   8    0.016     0.014  1.187
# 9   9    0.014     0.011  1.253
# 10 10    0.011     0.009  1.151
# 11 11    0.009     0.009  1.068
# 12 12    0.009     0.008  1.051
# 13 13    0.008     0.008  1.115
# 14 14    0.008     0.007  1.061
# 15 15    0.007     0.007  1.034

# same result as original Y
er_demean <- eigen_ratio_select(Y_demean, rmax = 15)
er_demean$r_hat
round(er_demean$table, 3)

er_std <- eigen_ratio_select(Y_std, rmax = 15)
er_std$r_hat
# [1] 1
round(er_std$table, 3)
# r lambda_r lambda_r1  ratio
# 1   1    0.746     0.039 18.952
# 2   2    0.039     0.027  1.437
# 3   3    0.027     0.016  1.697
# 4   4    0.016     0.009  1.751
# 5   5    0.009     0.009  1.078
# 6   6    0.009     0.006  1.533
# 7   7    0.006     0.005  1.037
# 8   8    0.005     0.005  1.126
# 9   9    0.005     0.004  1.077
# 10 10    0.004     0.004  1.090
# 11 11    0.004     0.004  1.125
# 12 12    0.004     0.003  1.046
# 13 13    0.003     0.003  1.007
# 14 14    0.003     0.003  1.051
# 15 15    0.003     0.003  1.095


### growth ratio
gr <- growth_ratio_select(Y, rmax = 15)
gr$r_hat
# [1] 1
round(gr$table, 3)
# r numerator denominator growth_ratio
# 1   1     1.602       0.202        7.915
# 2   2     0.202       0.204        0.992
# 3   3     0.204       0.106        1.927
# 4   4     0.106       0.062        1.708
# 5   5     0.062       0.063        0.981
# 6   6     0.063       0.042        1.513
# 7   7     0.042       0.041        1.007
# 8   8     0.041       0.040        1.027
# 9   9     0.040       0.037        1.098
# 10 10     0.037       0.034        1.084
# 11 11     0.034       0.032        1.065
# 12 12     0.032       0.030        1.063
# 13 13     0.030       0.030        1.003
# 14 14     0.030       0.029        1.027
# 15 15     0.029       0.029        1.013

gr_demean <- growth_ratio_select(Y_demean, rmax = 15)
gr_demean$r_hat
round(gr_demean$table, 3)

gr_std <- growth_ratio_select(Y_std, rmax = 15)
gr_std$r_hat
round(gr_std$table, 3)
# r numerator denominator growth_ratio
# 1   1     1.618       0.197        8.220
# 2   2     0.197       0.205        0.961
# 3   3     0.205       0.099        2.075
# 4   4     0.099       0.067        1.474
# 5   5     0.067       0.044        1.530
# 6   6     0.044       0.040        1.083
# 7   7     0.040       0.037        1.084
# 8   8     0.037       0.034        1.098
# 9   9     0.034       0.034        1.013
# 10 10     0.034       0.032        1.035
# 11 11     0.032       0.031        1.048
# 12 12     0.031       0.031        0.990
# 13 13     0.031       0.029        1.061
# 14 14     0.029       0.028        1.050
# 15 15     0.028       0.027        1.042


### bootstrap method
out <- select_num_factors_bootstrap(
  X = Y,
  Kmax = 8,
  B = 1000,          # increase to 1000+ for serious use
  alpha = 0.05,
  multiplier = "normal",
  center = TRUE,
  scale. = FALSE,
  verbose = TRUE
)
# Testing H0: r <= 0 vs H1: r >= 1 
# p-value 0 
# 
# Testing H0: r <= 1 vs H1: r >= 2 
# p-value 0 
# 
# Testing H0: r <= 2 vs H1: r >= 3 
# p-value 0 
# 
# Testing H0: r <= 3 vs H1: r >= 4 
# p-value 0.016 
# 
# Testing H0: r <= 4 vs H1: r >= 5 
# p-value 0.061 










## ---------- latent factor time series ----------
threshold <- quantile(abs(f), (Ti-20)/Ti)
idx <- which(abs(f) >= threshold)
Factors[idx_2019_first + idx - 1, 1]
# [1] 20190102 20200317 20200324 20200413 20210106 20210126 20210127 20210128 20210129 20210203 20210208 20210216 20210223 20210224
# [15] 20210309 20210324 20210602 20210610 20220616 20231211
round(f[idx], 3)
# [1]  3.801  4.604 -2.908  3.335 -3.517  5.738  5.934 -5.323  5.271  3.930  3.158  3.463 -2.960  3.062  4.953 -3.801  2.998 -3.058  3.291
# [20] -3.124

dates0 <- as.Date(as.character(Factors[,1]), format="%Y%m%d")
dates1 <- dates0[idx_2019_first:idx_2023_last]

## Key event dates and labels
event_dates <- as.Date(c(
  "2019-01-02",
  "2020-03-17",
  "2020-03-24",
  "2021-01-06",
  "2021-01-27",
  "2021-03-09",
  "2021-06-02",
  "2022-06-16",
  "2023-12-11"
))

event_labels <- c(
  "Apple warning",
  "COVID crash",
  "Stimulus rebound",
  "Blue wave",
  "Meme frenzy",
  "Rates/reflation",
  "Meme resurgence",
  "Fed 75bp",
  "Year-end\n rally"
)

event_idx <- match(event_dates, dates1)
keep <- !is.na(event_idx)
event_idx <- event_idx[keep]
event_labels <- event_labels[keep]

xlim_ext <- range(dates1) + c(-10, 20)
ylim_top <- range(c(f, threshold, -threshold)) * 1.05

plot(dates1, f, type = "l",
     xlim = xlim_ext, ylim = ylim_top,
     lwd = 1.6,
     xlab = "Date", ylab = "Latent factor values")

abline(h = 0, lty = 1, lwd = 0.6, col = "gray80")
abline(h = threshold,  lty = 2, lwd = 1, col = "gray50")
abline(h = -threshold, lty = 2, lwd = 1, col = "gray50")

## extreme points
points(dates1[idx], f[idx], pch = 16, cex = 0.5)

## key events
points(dates1[event_idx], f[event_idx], pch = 21, bg = "white", cex = 1.1)

## label offsets
x_shift <- c(6, -6, 6, 6, 5, -6, 6, -6, -10)

y_shift <- ifelse(f[event_idx] > 0, 0.35, -0.35)

## 手动微调几个关键点（核心）
y_shift[1] <- 0.4
y_shift[5] <- 0.55    # Meme frenzy ↑
y_shift[6] <- 0.65    # Rates/reflation ↑↑（避免挤）
y_shift[4] <- -0.55   # Blue wave ↓
y_shift[9] <- -0.5    # Year-end rally ↓

## last label a bit lower to avoid crowding
y_shift[length(y_shift)] <- -0.45

par(xpd = FALSE)

for (j in seq_along(event_idx)) {
  segments(dates1[event_idx[j]], f[event_idx[j]],
           dates1[event_idx[j]] + x_shift[j],
           f[event_idx[j]] + y_shift[j],
           lwd = 0.8, col = "gray50")
}

adj_vec <- ifelse(x_shift > 0, 0, 1)
adj_vec[9] <- 1   # 强制最后一个靠右对齐（往左贴）

## first 8 labels
text(dates1[event_idx[-length(event_idx)]] + x_shift[-length(x_shift)],
     f[event_idx[-length(event_idx)]] + y_shift[-length(y_shift)],
     labels = event_labels[-length(event_labels)],
     cex = 0.75,
     adj = adj_vec[-length(adj_vec)])

## last label
text(dates1[event_idx[length(event_idx)]] - 8,
     f[event_idx[length(event_idx)]] - 0.45,
     labels = event_labels[length(event_labels)],
     cex = 0.75,
     adj = 1)

mtext("(b)", side = 3, adj = 0, line = 0.2, cex = 0.95)