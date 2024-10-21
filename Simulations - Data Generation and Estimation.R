library(MASS)
library(Rcpp)
library(survival)
library(xtable)

sourceCpp("Prevalents.cpp")

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
MADVar = function(x) #estimating the variance using MAD
{
  K = 1.4826
  (median(abs(x-median(x))) * K)^2
}

N = 200 #number of simulations
n = 1500 #sample size
numPairs = 50 # user-defined number of sampled pairs per observation. 

resPairwise = resPL =  resL2 = matrix(NA,N,length(b_dis))
SDBoot1 = SDBoot2 = SDBoot3 = SDBoot3MAD = matrix(NA,N,length(b_dis))
CoverBoot1 = CoverBoot2 = CoverBoot3 = CoverBoot3MAD = matrix(NA,N,length(b_dis))
seedBase = 1
i = 1
while(i <= N)
{
  ##### DATA GENERATION #### 
  
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
  
  ### ESTIMATION PROCEDURE  ###
  
  #fitting marginal Cox models using partial likelihood. Warnings are produced because the data include prevalent observations which standard partial likelihood estimation cannot use. This is expected and therefore suppressed. 
  suppressWarnings({  
    fit1 = coxph(Surv(R, V, D1, type = "counting") ~ X)
    fit2 = coxph(Surv(R, V, D2, type="counting") ~ X)
    fit3 = coxph(Surv(pmax(R,V), W, D3, type = "counting") ~ X + V, subset = D1)
    fitC = coxph(Surv(R,V,1-D1-D2,type = "counting") ~ X)
  })
  
  resPL[i,] = coef(fit1)  # the partial likelihood estimator for transition 1->2
  
  #Getting the estimated baseline hazard functions:  
  suppressWarnings({  
    H012_estimated = stepfun(basehaz(fit1,centered = F)[,2],c(0,basehaz(fit1,centered = F)[,1]))
    H013_estimated = stepfun(basehaz(fit2,centered = F)[,2],c(0,basehaz(fit2,centered = F)[,1]))
    H023_estimated = stepfun(basehaz(fit3,centered = F)[,2],c(0,basehaz(fit3,centered = F)[,1]))
    H0C_estimated = stepfun(basehaz(fitC,centered = F)[,2],c(0,basehaz(fitC,centered = F)[,1]))
  })

  elp13 = exp(X %*% coef(fit2))
  elp23 = exp(cbind(X,V) %*% coef(fit3))
  elpC = exp(X %*% coef(fitC))
  
  weights = rep(1,n)
 
  ZetaTerms = getZetaTermsC(R = R, V = V, D1 = D1, D2 = D2, X23 = as.matrix(X), 
                            elp13 = elp13, elp23 = elp23, elpC = elpC,
                            H013V = H013_estimated(V), H023V = H023_estimated(V),
                            H023R = H023_estimated(R), H0CV = H0C_estimated(V),
                            H0CR = H0C_estimated(R), b23 = coef(fit3),numPairsPerElement = numPairs)
  
  resPairwiseOpt = try(optim(coef(fit1),getPairwiseLogLike,X12=X,D1=D1,
                         ZetaTerms=ZetaTerms,H012V = H012_estimated(V),
                         method = "BFGS", gr = getPairwiseDeriv, 
                         numPairsPerElement = numPairs, weights = weights),T)
  resPairwise[i,] = resPairwiseOpt$par
  
 
  # estimation using Likelihood 2 while plugging-in nuisance parameters:
  # (this is not the proposed estimation procedure, but we provide simulations
  # demonstrating its bias when recruitment age does not start from zero)
  
  # Events12sorted = sort(V[D1])
  # jumps12 = diff(c(0,H012_estimated(Events12sorted)))
  # elp23_noV = exp(X %*% coef(fit3)[names(coef(fit3)) != "V"])  ## Likelihood 2 requires that V is not added as a covariate here
  # 
  # resL2Opt = try(optim(coef(fit1),getLogLike2,X12 = X,D1=D1,H012V = H012_estimated(V),
  #                      R = R, H012R = H012_estimated(R), H013R = H013_estimated(R),
  #                      elp13 = elp13, elp23 = elp23_noV, jumps12 = jumps12, H012T = H012_estimated(Events12sorted),
  #                      H013T = H013_estimated(Events12sorted), T12sorted = Events12sorted,
  #                      H023R = H023_estimated(R), H023T = H023_estimated(Events12sorted),
  #                      b23V = coef(fit3)["V"],gr = getLogLike2Deriv, method = "BFGS"),T)
  # resL2[i,] = resL2Opt$par
  
########## Variance Estimation #########
  ## Required for bootstrap 2 and 3:
  prev = which(V < R)
  ninc = n - length(prev)
  ord1 = order(c(V[-prev],R[-prev]))
  VRord1 = c(V[-prev],R[-prev])[ord1]
  ord3 = order(c(W[D1],pmax(R,V)[D1]))
  WRord = c(W[D1],pmax(R,V)[D1])[ord3]
  D1ord = c(D1[-prev],rep(0,ninc))[ord1]
  D2ord = c(D2[-prev],rep(0,ninc))[ord1]
  D3ord = c(D3[D1],rep(0,sum(D1)))[ord3]
  DCord = c((1-D1-D2)[-prev],rep(0,ninc))[ord1]
  VorR = c(rep(T,ninc),rep(F,ninc))[ord1]
  WorR = c(rep(T,sum(D1)),rep(F,sum(D1)))[ord3]
  ##
  
  b = 1
  BootBaseSeed = 1
  N_Boot = 100 # number of bootstrap repetitions
  boot_hess_scores = matrix(NA,N_Boot,ncol(X)) # For Bootstrap 3
  boot_betas = matrix(NA,N_Boot,ncol(X)) # For Bootstrap 2
  while(b <= N_Boot)
  {
    set.seed(BootBaseSeed + b)
    ###### sampling from beta and H0  - required for Bootstrap 2 and 3 ######
    weights = rexp(n)
    b_dis_samp = mvrnorm(1,mu = coef(fit1),Sigma = vcov(fit1))
    b_dth1_samp = mvrnorm(1,mu = coef(fit2),Sigma = vcov(fit2))
    b_dth2_samp = mvrnorm(1,mu = coef(fit3),Sigma = vcov(fit3))
    b_cens_samp = mvrnorm(1,mu = coef(fitC),Sigma = vcov(fitC))
    
    elp12_boot = exp(X %*% b_dis_samp)
    elp13_boot = exp(X %*% b_dth1_samp)
    elp23_boot = exp(cbind(X,V) %*% b_dth2_samp)
    elpC_boot = exp(X %*% b_cens_samp)
    
    elp12_boot_ord = c(elp12_boot[-prev],elp12_boot[-prev])[ord1]
    elp13_boot_ord = c(elp13_boot[-prev],elp13_boot[-prev])[ord1]
    elp23_boot_ord = c(elp23_boot[D1],elp23_boot[D1])[ord3]
    elpC_boot_ord = c(elpC_boot[-prev],elpC_boot[-prev])[ord1]
    
    weightsOrd1 = c(weights[-prev],weights[-prev])[ord1]
    weightsOrd3 = c(weights[D1],weights[D1])[ord3]
    
    jumps12_boot = getBreslow(elp12_boot_ord,VorR,D1ord,weightsOrd1)
    jumps13_boot = getBreslow(elp13_boot_ord,VorR,D2ord,weightsOrd1)
    jumps23_boot = getBreslow(elp23_boot_ord,WorR,D3ord,weightsOrd3)
    jumpsC_boot = getBreslow(elpC_boot_ord,VorR,DCord,weightsOrd1)
    
    H012_boot = stepfun(sort(V[-prev][D1[-prev]]),c(0,cumsum(jumps12_boot)))
    H013_boot = stepfun(sort(V[D2]),c(0,cumsum(jumps13_boot)))
    H023_boot = stepfun(sort(W[D3]),c(0,cumsum(jumps23_boot)))
    H0C_boot = stepfun(sort(V[(1-D1-D2)==1]),c(0,cumsum(jumpsC_boot)))
    
    ZetaTerms_boot = getZetaTermsC(R = R, V = V, D1 = D1, D2 = D2, X23 =  as.matrix(X),
                                   elp13 = elp13_boot, elp23 = elp23_boot, elpC = elpC_boot,
                                   H013V = H013_boot(V), H023V = H023_boot(V),
                                   H023R = H023_boot(R), H0CV = H0C_boot(V), 
                                   H0CR = H0C_boot(R), b23 = b_dth2_samp, numPairsPerElement = numPairs)
    
    # Bootstrap 2: #
    OptFitBoot = try(optim(resPairwise[i,],getPairwiseLogLike,X12 = X,D1 = D1,
                           ZetaTerms = ZetaTerms_boot, H012V = H012_boot(V),
                           method = "BFGS", gr = getPairwiseDeriv,
                           numPairsPerElement = numPairs,weights=weights),T)
    if(class(OptFitBoot) == "try-error") #If optimization failed, skip this bootstrap sample
    {
      BootBaseSeed = BootBaseSeed + 1
      next
    }
    
    boot_betas[b,] = OptFitBoot$par
    ####
    
    # Bootstrap 3: #
    weights = rep(1,n)
    score = getPairwiseDeriv(resPairwise[i,], X12 = X, D1 = D1, ZetaTerms = ZetaTerms_boot,
                             H012V = H012_boot(V),numPairsPerElement = numPairs,
                             weights=weights)
    hess = getPairwiseHessian(resPairwise[i,], X12 = X, D1 = D1, ZetaTerms = ZetaTerms_boot,
                              H012V = H012_boot(V),numPairsPerElement = numPairs,weights = weights)
    boot_hess_scores[b,] = solve(hess) %*% score
    ####
    
    b = b + 1
  }
  
  VarBoot2 = cov(boot_betas)  # estimated variance matrix based on Bootstrap 2
  SDBoot2[i,] = sqrt(diag(VarBoot2))
  ## for Bootstrap 3: ##
  #### Deriving the empirical cross-covariance matrix psi_ij and psi_il.  
  #numPairsPerElement = number of sampled pairs per observations used to derive the point estimates
  #numPairsPerElementVar = number of sampled pairs per observations for estimating the variance, can be the same, or smaller than numPairsPerElement.
  CrossCov = getCrossCov(resPairwise[i,], 1*D1, X12 = X, ZetaTerms = ZetaTerms, 
                         H012V = H012_estimated(V), numPairsPerElement = numPairs,
                         numPairsPerElementVar =  numPairs)
  CrossCov = CrossCov + t(CrossCov) - diag(diag(CrossCov))
  
  hess = getPairwiseHessian(b12 = resPairwise[i,], D1 = D1, X12 = X, ZetaTerms = ZetaTerms,
                            H012V = H012_estimated(V),numPairsPerElement = numPairs, weights = weights)
  VarBoot3 = solve(hess) %*% CrossCov %*% solve(hess) + cov(boot_hess_scores) # estimated variance matrix based on Bootstrap 3
  VarBoot3MAD = diag(solve(hess) %*% CrossCov %*% solve(hess)) + apply(boot_hess_scores,2,MADVar) # estimated variance matrix based on Bootstrap 3 - using MAD. This is for setting A with n = 1500.
  SDBoot3[i,] = sqrt(diag(VarBoot3))
  SDBoot3MAD[i,] = sqrt(VarBoot3MAD)
  
  ############# Bootstrap 1 - full weighted bootstrap: ###########
  
  boot_betas = matrix(NA,N_Boot,ncol(X)) 
  b = 1
  
  while(b <= N_Boot)
  {
    set.seed(BootBaseSeed + b)
    weights = rexp(n)
    suppressWarnings({  
    fit1_boot = try(coxph(Surv(R, V, D1, type = "counting") ~ X,weights = weights),T)
    fit2_boot = try(coxph(Surv(R, V, D2, type="counting") ~ X,weights = weights),T)
    fit3_boot = try(coxph(Surv(pmax(R,V), W, D3, type = "counting") ~ X + V, subset = D1,weights = weights),T)
    fitC_boot = try(coxph(Surv(R,V,1-D1-D2,type = "counting") ~ X,weights = weights),T)
    })
    if(any(c(class(fit1_boot),class(fit2_boot),class(fit3_boot),class(fitC_boot)) == "try-error"))
    {
      BootBaseSeed = BootBaseSeed + 1
      next
    }
    suppressWarnings({  
    H012_boot = stepfun(basehaz(fit1_boot,centered = F)[,2],c(0,basehaz(fit1_boot,centered = F)[,1]))
    H013_boot = stepfun(basehaz(fit2_boot,centered = F)[,2],c(0,basehaz(fit2_boot,centered = F)[,1]))
    H023_boot = stepfun(basehaz(fit3_boot,centered = F)[,2],c(0,basehaz(fit3_boot,centered = F)[,1]))
    H0C_boot = stepfun(basehaz(fitC_boot,centered = F)[,2],c(0,basehaz(fitC_boot,centered = F)[,1]))
    })
    
    elp13_boot = exp(X %*% coef(fit2_boot))
    elp23_boot = exp(cbind(X,V) %*% coef(fit3_boot))
    elpC_boot = exp(X %*% coef(fitC_boot))
    
    ZetaTerms_boot = getZetaTermsC(R = R, V = V, D1 = D1, D2 = D2, X23 =  as.matrix(X),
                                   elp13 = elp13_boot, elp23 = elp23_boot, elpC = elpC_boot,
                                   H013V = H013_boot(V), H023V = H023_boot(V),
                                   H023R = H023_boot(R), H0CV = H0C_boot(V), 
                                   H0CR = H0C_boot(R), b23 = coef(fit3_boot), numPairsPerElement = numPairs)
    
    OptFitBoot = try(optim(coef(fit1_boot),getPairwiseLogLike,X12 = X,D1 = D1,ZetaTerms = ZetaTerms_boot,
                           H012V = H012_boot(V),method = "BFGS", gr = getPairwiseDeriv, numPairsPerElement = numPairs,
                           weights=weights),T)
    if(class(OptFitBoot) == "try-error")
    {
      BootBaseSeed = BootBaseSeed + 1
      next
    }
    
    boot_betas[b,] = OptFitBoot$par
    b = b + 1
  }
  
  VarBoot1 = var(boot_betas)
  SDBoot1[i,] = sqrt(diag(VarBoot1))
 
  # getting bootstrap coverage rates:
  
  U1 = resPairwise[i,] + qnorm(0.975) * SDBoot1[i,]
  L1 = resPairwise[i,] - qnorm(0.975) * SDBoot1[i,]
  CoverBoot1[i,] = b_dis < U1 & b_dis > L1
  U2 = resPairwise[i,] + qnorm(0.975) * SDBoot2[i,]
  L2 = resPairwise[i,] - qnorm(0.975) * SDBoot2[i,]
  CoverBoot2[i,] = b_dis < U2 & b_dis > L2
  U3 = resPairwise[i,] + qnorm(0.975) * SDBoot3[i,]
  L3 = resPairwise[i,] - qnorm(0.975) * SDBoot3[i,]
  CoverBoot3[i,] = b_dis < U3 & b_dis > L3
  U3MAD = resPairwise[i,] + qnorm(0.975) * SDBoot3MAD[i,]
  L3MAD = resPairwise[i,] - qnorm(0.975) * SDBoot3MAD[i,]
  CoverBoot3MAD[i,] = b_dis < U3 & b_dis > L3

  i = i + 1
}

### summarizing the results

bmat = matrix(b_dis,nrow = N,ncol = length(b_dis),byrow = T)

RE = apply((resPL-bmat)^2,2,mean) / apply((resPairwise-bmat)^2,2,mean)
Dat = rbind(b_dis,apply(resPL,2,mean),apply(resPairwise,2,mean),
            apply(resPL,2,sd),apply(resPairwise,2,sd),
            RE,
            apply(SDBoot1,2,mean),
            apply(SDBoot2,2,mean),
            apply(SDBoot3,2,mean),
            apply(SDBoot3MAD,2,mean),
            apply(CoverBoot1,2,mean),
            apply(CoverBoot2,2,mean),
            apply(CoverBoot3,2,mean),
            apply(CoverBoot3MAD,2,mean))

xtable(Dat,digits = 2)


