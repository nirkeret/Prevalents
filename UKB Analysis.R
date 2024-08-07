library(Rcpp)
library(survival)
library(data.table)
sourceCpp("Prevalents.cpp")
UKB = fread("SynthUKB.csv",data.table = F)

ncol = ncol(UKB)
n = nrow(UKB)
set.seed(1)

#shuffling the order
UKB = UKB[sample(1:n),]

#scaling SNPs and PCs
UKB[,-(1:7)] = scale(UKB[,-(1:7)])

numPairs = 25
N_Boot = 500

cox_estimates = pairwise_estimates = pairwise_SE_boot2 = pairwise_SE_boot3 = matrix(NA,nrow = 24,ncol = 8)

prev = which(UKB$V < UKB$R)
inc = which(UKB$V > UKB$R & UKB$delta1 == 1)

for(i in 1:24)
{
  SNPind = i-1 + 9
  
  fit1 = coxph(Surv(R,V,delta1,type = 'counting') ~ sex + PC1 +
                 PC2 + PC3 + PC4 + PC5 + PC6 + UKB[,SNPind], data = UKB)
  fit2 = coxph(Surv(R,V,delta2,type = 'counting') ~ sex + PC1 +
                 PC2 + PC3 + PC4 + PC5 + PC6 + UKB[,SNPind], data = UKB)
  fit3 = coxph(Surv(pmax(R,V),W,delta3,type = 'counting') ~ sex + PC1 +
                 PC2 + PC3 + PC4 + PC5 + PC6 + UKB[,SNPind] + V, subset = which(delta1==1),data = UKB)
  
  H012_estimated = stepfun(basehaz(fit1,centered = F)[,2],c(0,basehaz(fit1,centered = F)[,1]))
  H013_estimated = stepfun(basehaz(fit2,centered = F)[,2],c(0,basehaz(fit2,centered = F)[,1]))
  H023_estimated = stepfun(basehaz(fit3,centered = F)[,2],c(0,basehaz(fit3,centered = F)[,1]))
  
  cox_estimates[i,] = coef(fit1)
  
  mtrx12 = model.matrix(~ sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 
                        + UKB[,SNPind],data = UKB)[,-1]
  
  mtrx13 = model.matrix(~ sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 
                        + UKB[,SNPind],data = UKB)[,-1]
  elp13 = as.vector(exp(mtrx13 %*% coef(fit2)))
  
  mtrx23 = model.matrix(~ sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 +
                          UKB[,SNPind] + V,data = UKB)[,-1]
  
  elp23 = as.vector(exp(mtrx23 %*% coef(fit3)))
  
  
  weights = rep(1,nrow(UKB))
  
  ZetaTerms = with(UKB,getZetaTerms(R, V, D1=delta1, D2=delta2, X23=mtrx23, 
                                    elp13=elp13, elp23=elp23, H013V = H013_estimated(V), H023V=H023_estimated(V), 
                                           H023R=H023_estimated(R), b23 = coef(fit3),numPairsPerElement=numPairs))
  
  resPair_opt = with(UKB,optim(coef(fit1),getPairwiseLogLike,
                                     X12=mtrx12,D1=delta1,
                                     ZetaTerms=ZetaTerms,H012V = H012_estimated(V),
                                     method = "BFGS", gr = getPairwiseDeriv,
                                     numPairsPerElement = numPairs,weights = weights))
  
  pairwise_estimates[i,] = resPair_opt$par
  
  numPairsVar = 25
  ZetaTerms_var_ind = rep(1:numPairsVar,nrow(UKB)) + rep(numPairs*(0:(nrow(UKB) - 1)),each = numPairsVar)
  
  CrossCovMat = with(UKB,getCrossCov(b12=pairwise_estimates[i,], D1=delta1, X12=mtrx12, ZetaTerms = ZetaTerms[ZetaTerms_var_ind], H012V = H012_estimated(V), numPairsPerElement=numPairs,numPairsPerElementVar=numPairsVar))
  CrossCovMat = CrossCovMat + t(CrossCovMat) - diag(diag(CrossCovMat))
  
  hess = with(UKB,getPairwiseHessian(b12=pairwise_estimates[i,], D1=delta1, X12=mtrx12, ZetaTerms=ZetaTerms, H012V = H012_estimated(V),numPairsPerElement=numPairs,weights=weights))
  estVarMat = solve(hess) %*% CrossCovMat %*% solve(hess)
  
  prev = with(UKB,which(V < R))
  ninc = n - length(prev)
  ord1 = with(UKB,order(c(V[-prev],R[-prev])))
  VRord1 = with(UKB,c(V[-prev],R[-prev])[ord1])
  ord3 = with(UKB,order(c(W[delta1],pmax(R,V)[delta1])))
  WRord = with(UKB,c(W[delta1],pmax(R,V)[delta1])[ord3])
  delta1ord = with(UKB,c(delta1[-prev],rep(F,ninc))[ord1])
  delta2ord = with(UKB,c(delta2[-prev],rep(F,ninc))[ord1])
  delta3ord = with(UKB,c(delta3[delta1],rep(F,sum(delta1)))[ord3])
  VorR = c(rep(T,ninc),rep(F,ninc))[ord1]
  WorR = with(UKB,c(rep(T,sum(delta1)),rep(F,sum(delta1)))[ord3])
  
  boot_betas_samp = matrix(NA,nrow = N_Boot,ncol = 8)
  boot_hess_scores = matrix(NA,N_Boot,length(coef(fit1)))
  seedbase = 0
  b = 1
  for(b in b:N_Boot)
  {
    while(T)
    {
      set.seed(b + seedbase)
      weights = rexp(n)
     
      b_dis_samp = mvrnorm(1,mu = coef(fit1),Sigma = vcov(fit1))
      b_dth1_samp = mvrnorm(1,mu = coef(fit2),Sigma = vcov(fit2))
      b_dth2_samp = mvrnorm(1,mu = coef(fit3),Sigma = vcov(fit3))
      
      elp12_boot = exp(mtrx12 %*% b_dis_samp)
      elp13_boot = exp(mtrx13 %*% b_dth1_samp)
      elp23_boot = exp(mtrx23 %*% b_dth2_samp)
      
      elp12_boot_ord = c(elp12_boot[-prev],elp12_boot[-prev])[ord1]
      elp13_boot_ord = c(elp13_boot[-prev],elp13_boot[-prev])[ord1]
      elp23_boot_ord = with(UKB,c(elp23_boot[delta1],elp23_boot[delta1])[ord3])
      
      weightsOrd1 = c(weights[-prev],weights[-prev])[ord1]
      weightsOrd3 = with(UKB,c(weights[delta1],weights[delta1])[ord3])
      
      jumps12_boot = getBreslow(elp12_boot_ord,VorR,delta1ord,weightsOrd1)
      jumps13_boot = getBreslow(elp13_boot_ord,VorR,delta2ord,weightsOrd1)
      jumps23_boot = getBreslow(elp23_boot_ord,WorR,delta3ord,weightsOrd3)
      
      H012_boot = with(UKB,stepfun(sort(V[-prev][delta1[-prev]]),c(0,cumsum(jumps12_boot))))
      H013_boot = with(UKB,stepfun(sort(V[delta2]),c(0,cumsum(jumps13_boot))))
      H023_boot = with(UKB,stepfun(sort(W[delta3]),c(0,cumsum(jumps23_boot))))
      
      
      ZetaTerms_sampboot = with(UKB,getZetaTerms(R=R, V=V, D1=delta1, D2=delta2, X23=mtrx23, 
                                                 elp13=elp13_boot, elp23=elp23_boot, H013V=H013_boot(V), H023V=H023_boot(V), 
                                                 H023R=H023_boot(R), b23=b_dth2_samp,numPairsPerElement=numPairs))
      #this part is required just for Boot2:
      
      bootOpt = try(with(UKB,optim(b_dis_samp,getPairwiseLogLike
                                   ,X12=mtrx12,D1=delta1,
                                   ZetaTerms=ZetaTerms_sampboot,
                                   H012V = H012_boot(V),
                                   method = "BFGS", 
                                   gr = getPairwiseDeriv,
                                   numPairsPerElement = numPairs,
                                   weights = weights)),T)
      if(class(bootOpt) == "try-error")
      {
        seedbase = seedbase + 1
        errorCount[i] = errorCount[i] + 1
        next
      }
      
      boot_betas_samp[b,] = bootOpt$par
      
      break
    }
    
    #This part is required for Boot3:
    
    score = with(UKB,getPairwiseDeriv(coef(fit1),X12=mtrx12,D1=delta1,ZetaTerms = ZetaTerms_sampboot,H012V = H012_boot(V),numPairsPerElement = numPairs,weights=rep(1,n)))
    hessian = with(UKB,getPairwiseHessian(coef(fit1),X12=mtrx12,D1=delta1,ZetaTerms=ZetaTerms_sampboot,H012V = H012_boot(V),numPairsPerElement = numPairs,weights=rep(1,n)))
    boot_hess_scores[b,] = solve(hessian) %*% score

  }
  
  Boot3Var = cov(boot_hess_scores) + estVarMat
  pairwise_SE_boot3[i,] = sqrt(diag(Boot3Var))
  
  Boot2Var = cov(boot_betas_samp) 
  pairwise_SE_boot2[i,] = sqrt(diag(Boot2Var))
  
  print(i)
}

rownames(pairwise_estimates) = rownames(pairwise_SE_boot2) = rownames(pairwise_SE_boot3) = colnames(UKB)[8:31]


