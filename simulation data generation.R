library(MASS)

# coefficients for transition healthy -> illness
b_dis = c(2,-1.5,0.1,-0.5,1,-2.5,-1,0) # Settings A and B
# b_dis = c(2,-1,0.1,-0.5,1,-1,-1,0) # Setting C

# coefficients for transition healthy -> death, settings A and B
b_dth1 = c(0.3,-0.2,0.4,0.7) 
b_dth1_ind = c(1,5,6,8) # 

# coefficients for transition illness -> death, settings A and B
b_dth2 = c(-0.3,0.9,0.05)
b_dth2_ind = c(7,8)

# coefficients for censoring
b_cens = c(1.5,0.5) #setting B
b_cens_ind = c(2,5) #settings B and C

# baseline hazard function constants
H_const_dis = 0.02 # settings A and B
# H_const_dis = 0.01 # setting C
H_const_dth1 = 0.02
H_const_dth2 = 0.05
H_const_cens = 0.05

UnitRescale = function(x) (x-min(x))/(max(x) - min(x))

N = 200 #number of simulations
n = 1500 #sample size
seedBase = 1
i = 1
while(i <= N)
{
  set.seed(seedBase + i)
  #nn is the number of observations generated, out of which we sample n observations that satisfy death > recruitment (T2>R).   
  nn = 5000 # if n = 1500
  # nn = 50000 # if n = 10000

  ## generating covariates:
  
  ## setting A:
  X = matrix(NA,nn,8)
  X[,1] = rgamma(nn,2,6) 
  X[,2] = rgeom(nn,1/10)
  X[,3] = rexp(nn,0.25)
  X[,4] =  rbeta(nn,2,8)
  X[,5] = rnorm(nn,0,2)
  X[,6] = rweibull(nn,3,4)
  X[,7] = rpois(nn,5)
  X[,8] = runif(nn)
  X = apply(X,2,UnitRescale)
  
  # ## settings B and C
  # rho = 0.8
  # S = matrix(rho,8,8)
  # diag(S) = 1
  # X = pnorm(mvrnorm(nn,mu=rep(0,8),Sigma = S))

  ## sampling disease and death-before-disease times:
  
  T1 = rexp(nn,H_const_dis * exp(X %*% b_dis))
  
  T2 = rexp(nn,H_const_dth1 * exp(X[,b_dth1_ind] %*% b_dth1)) #settings A and B
  
  #setting C:
  # mu = sin(X[,b_dth1_ind[1]] * pi) + 2 * abs(X[,b_dth1_ind[2]] - 0.5) + X[,b_dth1_ind[3]]^3
  # T2 = rexp(nn,0.04 / mu)

  ## sampling censoring times:
  C = rexp(nn,H_const_cens) # setting A
  # C = rexp(nn,H_const_cens * exp(X[,b_cens_ind] %*% b_cens)) #setting B
  # C = rlnorm(nn,meanlog = 3*abs(X[,b_cens_ind[1]]-0.5) + 2*X[b_cens_ind[2]], sdlog = 1.5) #setting C
  
  ## sampling recruitment age:
  
  #setting A:
  u = runif(nn)
  a = 4; b = 22; c = 13
  R = ifelse(u>0.5,b - sqrt((1-u)*(b-a)*(b-c)),a + sqrt(u * (b-a)*(c-a)))
  
  # R = pmax(0,1 + 5 * X[,1] + 7*X[,2]  + 10*X[,6] + rnorm(nn)) #setting B
  # R = pmax(0,1 + 5 * X[,1] + 6*X[,2]  + 4*X[,6] + rnorm(nn)) #setting C
  
  
  V = pmin(T1,T2,R + C)
  D1 = V == T1
  D2 = V == T2
  n_dis = sum(D1)
  
  ## sampling death after disease:

  #settings A and B:
  rates = as.vector(H_const_dth2 * exp(cbind(X[,b_dth2_ind],V)[D1,] %*% b_dth2))
  T2[D1] = rexp(n_dis,rates) + V[D1]  
  
  #setting C:
  # mu = 0.5 + (cos(pi*X[,b_dth2_ind[1]])^2 + 2*abs(X[,b_dth2_ind[2]]-0.5))[D1] + sqrt(V[D1])/3
  # T2[D1] = rgamma(n_dis,scale = 3, shape = mu) + V[D1]  

  D3 = (T2 < R + C) & D1
  W = pmin(T2, R + C)
  obs = which(T2 > R)
  samp = sample(obs,n)  # sample n observations out of nn, satisfying T2>R.

  ## keep the n sampled observations
  V = V[samp]
  R = R[samp]
  X = X[samp,]
  D1 = D1[samp]
  D2 = D2[samp]
  D3 = D3[samp]
  W = W[samp]
  
  ### here comes the estimation procedure, provided in a separate file
  
  i = i + 1
}
