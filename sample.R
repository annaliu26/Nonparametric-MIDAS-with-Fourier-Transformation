rm(list=ls())
list.of.packages <- c("igraph", "Rcpp", "RcppArmadillo", "foreach")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(igraph); library(Rcpp); library(RcppArmadillo); library(foreach)

set.seed(100)

FUNCcode_dir = 'I:/My Drive/Research/Nonparametric estimation/functions.R'
CPPcode_dir = 'I:/My Drive/Research/Nonparametric estimation/admmmcp.cpp'
source(FUNCcode_dir)  ## load Nonpara MIDAS code
sourceCpp(CPPcode_dir)  ## load CPP code for updating iteration

############################################################################################

## Set the frequency ratio, parameters of Fourier transformation, all tuning parameters
## and final iteration

L = 2; K = 3; v.len = L+2*K+1;
lambda2 = 1; theta = 2
theta.gamma = 0
lambda1s = seq(2, 3, 0.1)  ## choice of lambda1
iterations = 20000  ## final iteration num

############################################################################################

## load dataset, Y, X
## All response should be a column vector
## Predictors X should be a block diagonal matrix, each subject is one block

## Example of dataset below
## This is for simulation sample. In actual, we only know N = total num of subjects
T = 100  ## sample size in each subject. 
m = 40  ## frequency ratio
R = 15
N = 2*R  ## total number of subjects

SS = as.matrix(Fourier(m, K, L))

data = read.csv('data.csv', header = TRUE)  ## Example data = y|X
df = rep(NA,N)  ## degree of freedom
X = matrix(rep(0, T*v.len*N^2), nrow = N*T)  ## diagonal matrix
for (r in 1:N)
{
  X.sub = data[1:T+(r-1)*T,-1]
  X.trans = as.matrix(X.sub)%*%SS
  X[1:T+(r-1)*T,1:v.len+(r-1)*v.len] = X.trans
  df[r] = tr(X.trans%*%solve(t(X.trans)%*%X.trans)%*%t(X.trans))
}
y = as.vector(data[,1])  ## y


############################################################################################

## Set the penalty as 0 matrix in Zhu&Qu(2015)
D = matrix(rep(0,v.len*v.len), nrow=v.len);
diagD <- kronecker(diag(1, N), D);

mem_final <- matrix(rep(NA, length(lambda1s)*N), nrow=N)
MSE = k_final = BIC = final_iter <- rep(NA, length(lambda1s))
beta_s = matrix(rep(NA, N*v.len*length(lambda1s)), nrow=N*v.len)
beta.est = list()

## all results for all settings of lambda1
matrix = foreach(iii = 1:length(lambda1s),.combine=cbind) %dopar% 
{
  lambda1 = lambda1s[iii]; index = t(combn(N,2));
  
  ## if LS estimator cannot be calculated, try to use penalty matrix
  B_ini01 = update_B_ini(X, diagD, y, N, theta.gamma, index, lambda0 = 1) #lambda0 = theta in initial beta (Zhu&Qu)
  
  sol_final = prclust_admm(X, y, diagD, B_ini01, index,
                            gamma1 = theta.gamma, gamma2 = lambda1, 
                            theta = lambda2, tau = theta, N, v.len,
                            max_iter = iterations,
                            eps_abs = 1e-4, eps_rel = 1e-4)

  final_iter[iii] = sol_final$'iterations';  ## record num of iters
  beta.est[[iii]] = sol_final$'B_hat';  ## record estimated coeffs in each iter
  beta_s[,iii] = sol_final$'B_hat'[, final_iter[iii]]  ## final estimated coeffs
  
  Ad_final <- create_adjacency(sol_final$V, N);
  G_final <- graph.adjacency(Ad_final, mode = 'upper')
  cls_final <- components(G_final);  
  mem_final[,iii]<- cls_final$membership;  ## clustering panels
  k_final[iii]<- cls_final$no  ## total panel num
  
  BIC[iii] = log(t(y-X%*%beta_s[,iii])%*%(y-X%*%beta_s[,iii])/N)+log(N)*(k_final[iii]+v.len)/N  ## BIC for lambda1
}

BIC.select = which(BIC == min(BIC))
result = list(Iteration = final_iter[BIC.select],
              Beta.est = beta_s[,BIC.select],
              Panels = mem_final[,BIC.select],
              Total.panel = k_final[BIC.select],
              BIC.lambda1 = BIC[BIC.select])

