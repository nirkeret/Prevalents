library(Rcpp)
library(survival)
library(MASS)
sourceCpp("Prevalents.cpp")

## Assume in the following that:
# The data are organized in a random order. If not, should be shuffled.
# V = first observed age (minimum of disease diagnosis age, death age or censoring)
# W = second observed age, after illness diagnosis (minimum of death age or censoring)
# D1 = Event indicator, takes 1 if V is an illness age.
# D2 = Event indicator, takes 1 if V is a death age.
# D3 = Event indicator, takes 1 if W is a death age.
# X = matrix of covariates.
# R = recruitment age.

n = nrow(X) #sample size
numPairs = 50 # user-defined number of sampled pairs per observation. 

#fitting marginal Cox models using partial likelihood 
  
fit1 = coxph(Surv(R, V, D1, type = "counting") ~ X)
fit2 = coxph(Surv(R, V, D2, type="counting") ~ X)
# V or functional forms of it can be incorporated as covariates in the illness -> death transition:
fit3 = coxph(Surv(pmax(R,V), W, D3, type = "counting") ~ X + V, subset = D1)
# If censoring is assumed random and should be estimated too:
fitC = coxph(Surv(R,V,1-D1-D2,type = "counting") ~ X)

#Getting the estimated baseline hazard functions:
H012_estimated = stepfun(basehaz(fit1,centered = F)[,2],c(0,basehaz(fit1,centered = F)[,1]))
H013_estimated = stepfun(basehaz(fit2,centered = F)[,2],c(0,basehaz(fit2,centered = F)[,1]))
H023_estimated = stepfun(basehaz(fit3,centered = F)[,2],c(0,basehaz(fit3,centered = F)[,1]))
H0C_estimated = stepfun(basehaz(fitC,centered = F)[,2],c(0,basehaz(fitC,centered = F)[,1]))
  
elp13 = exp(X %*% coef(fit2))
elp23 = exp(cbind(X,V) %*% coef(fit3))
elpC = exp(X %*% coef(fitC))
  
ZetaTerms = getZetaTerms(R = R, V = V, D1 = D1, D2 = D2, X = as.matrix(X), 
                         elp13 = elp13, elp23 = elp23, elpC = elpC,
                         H013V = H013_estimated(V), H023V = H023_estimated(V),
                         H023R = H023_estimated(R), H0CV = H0C_estimated(V),
                         H0CR = H0C_estimated(R), b23 = coef(fit3),numPairsPerElement = numPairs)

weights = rep(1,n)
resPairOpt = try(optim(coef(fit1),getPairwiseLogLike,X=X,D1=D1,
                       ZetaTerms=ZetaTerms,H012V = H012_estimated(V),
                       method = "BFGS", gr = getPairwiseDeriv, 
                       numPairsPerElement = numPairs, weights = weights),T)
resPair = resPairOpt$par

########################### variance estimation
## "Bootstrap 3" - a fast bootstrap procedure. Recommended when the sample size is big. 
## "Bootstrap 2" - "Piggyback" bootstrap procedure. Recommended when the sample size is small. 
## "Bootstrap 1" - Full weighted bootstrap. Usually not recommended because it is slow.

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

  ZetaTerms_boot = getZetaTerms(R = R, V = V, D1 = D1, D2 = D2, X =  as.matrix(X),
                                elp13 = elp13_boot, elp23 = elp23_boot, elpC = elpC_boot,
                                H013V = H013_boot(V), H023V = H023_boot(V),
                                H023R = H023_boot(R), H0CV = H0C_boot(V), 
                                H0CR = H0C_boot(R), b23 = b_dth2_samp, numPairsPerElement = numPairs)
  
  # Bootstrap 2: #
  OptFitBoot = try(optim(resPair,getPairwiseLogLike,X = X,D1 = D1,
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
  score = getPairwiseDeriv(resPair, X = X, D1 = D1, ZetaTerms = ZetaTerms_boot,
                           H012V = H012_boot(V),numPairsPerElement = numPairs,
                           weights=weights)
  hess = getPairwiseHessian(resPair, X = X, D1 = D1, ZetaTerms = ZetaTerms_boot,
                               H012V = H012_boot(V),numPairsPerElement = numPairs,weights = weights)
  boot_hess_scores[b,] = solve(hess) %*% score
  ####
  
  b = b + 1
}

VarBoot2 = cov(boot_betas)  # estimated variance matrix based on Bootstrap 2

## for Bootstrap 3: ##
#### Deriving the empirical cross-covariance matrix psi_ij and psi_il.  
#numPairsPerElement = number of sampled pairs per observations used to derive the point estimates
#numPairsPerElementVar = number of sampled pairs per observations for estimating the variance, can be the same, or smaller than numPairsPerElement.
CrossCov = getCrossCov(resPair, 1*D1, X, ZetaTerms = ZetaTerms, 
                       H012V = H012_estimated(V), numPairsPerElement = numPairs,
                       numPairsPerElementVar =  numPairs)
CrossCov = CrossCov + t(CrossCov) - diag(diag(CrossCov))

hess = getPairwiseHessian(b12 = resPair, D1 = D1, X = X, ZetaTerms = ZetaTerms,
                          H012V = H012_estimated(V),numPairsPerElement = numPairs, weights = weights)
VarBoot3 = solve(hess) %*% CrossCov %*% solve(hess) + cov(boot_hess_scores) # estimated variance matrix based on Bootstrap 3

############# Bootstrap 1 - full weighted bootstrap: ###########

BootBaseSeed = 1
N_Boot = 100 # number of bootstrap repetitions
boot_betas = matrix(NA,N_Boot,ncol(X)) 
b = 1

while(b <= N_Boot)
{
  set.seed(BootBaseSeed + b)
  weights = rexp(n)
  fit1_boot = try(coxph(Surv(R, V, D1, type = "counting") ~ X,weights = weights),T)
  fit2_boot = try(coxph(Surv(R, V, D2, type="counting") ~ X,weights = weights),T)
  fit3_boot = try(coxph(Surv(pmax(R,V), W, D3, type = "counting") ~ X + V, subset = D1,weights = weights),T)
  fitC_boot = try(coxph(Surv(R,V,1-D1-D2,type = "counting") ~ X,weights = weights),T)

  if(any(c(class(fit1_boot),class(fit2_boot),class(fit3_boot),class(fitC_boot)) == "try-error"))
  {
    BootBaseSeed = BootBaseSeed + 1
    next
  }

  H012_boot = stepfun(basehaz(fit1_boot,centered = F)[,2],c(0,basehaz(fit1_boot,centered = F)[,1]))
  H013_boot = stepfun(basehaz(fit2_boot,centered = F)[,2],c(0,basehaz(fit2_boot,centered = F)[,1]))
  H023_boot = stepfun(basehaz(fit3_boot,centered = F)[,2],c(0,basehaz(fit3_boot,centered = F)[,1]))
  H0C_boot = stepfun(basehaz(fitC_boot,centered = F)[,2],c(0,basehaz(fitC_boot,centered = F)[,1]))

  elp13_boot = exp(X %*% coef(fit2_boot))
  elp23_boot = exp(cbind(X,V) %*% coef(fit3_boot))
  elpC_boot = exp(X %*% coef(fitC_boot))

  ZetaTerms_boot = getZetaTerms(R = R, V = V, D1 = D1, D2 = D2, X =  as.matrix(X),
                              elp13 = elp13_boot, elp23 = elp23_boot, elpC = elpC_boot,
                              H013V = H013_boot(V), H023V = H023_boot(V),
                              H023R = H023_boot(R), H0CV = H0C_boot(V), 
                              H0CR = H0C_boot(R), b23 = coef(fit3_boot), numPairsPerElement = numPairs)

  OptFitBoot = try(optim(coef(fit1_boot),getPairwiseLogLike,X = X,D1 = D1,ZetaTerms = ZetaTerms_boot,
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

VarBoot1 = cov(boot_betas)  # estimated variance matrix based on Bootstrap 1
