#include <Rcpp.h>
//#include<iostream>
//using namespace std;
using namespace Rcpp;

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

// In the following:
// V = vector of observed event times.
// R = vector of recruitment times.
// D1 = vector of event indicators. D1=1 if V is an observed event of type "disease".
// D2 = vector of event indicators. D2=1 if V is an observed event of type "death".
// X12,X23 = transition-specific covariate matrices (transitions 1->2 (healthy->disease) and 2->3 (disease->death)).
// elp12,elp13,elp23 = vector of exponents of transition-specific linear predictors (exp(x^b)).
// H012V,H013V,H023V,H0CV, H0CR, H023R = vectors of the transition-specific Breslow estimators at times V or R.
// b23 = partial-likelihood estimator of transition 2->3.
// numPairsPerElement = number of subsampled pairs per observation.
// numPairsPerElementVar = number of subsampled pairs per observation used for getting the cross-covariance matrix of the pairwise pseudolikelihood score terms, needed for one of the bootstrap procedures.
// ZetaTerms = vector of numPairsPerElement*n Zeta terms, as defined in the paper. Produced by the getZetaTerms/getZetaTermsC functions.
// weights = a vector of observation weights. Used for weighted bootstrap.

// [[Rcpp::export]]
// A function for getting the Zeta terms for all subsampled pairs, assuming censoring is not random.
NumericVector getZetaTerms(NumericVector R, NumericVector V, LogicalVector D1, LogicalVector D2, NumericMatrix X23 , NumericVector elp13, NumericVector elp23, NumericVector H013V, NumericVector H023V, NumericVector H023R, NumericVector b23, int numPairsPerElement) {
	int n = R.size();
	int nCov = X23.ncol();
	NumericVector resVec(n * numPairsPerElement);
	double loc = 0;
	for (int i = 0; i < n; ++i) {
		for (int p = 1; p <= numPairsPerElement; ++p) {
			int j = (i + p) % n;
			bool cond1 = (V[j] < R[i]) & !D1[j];
			bool cond2 = (V[i] < R[j]) & !D1[i];
			if (cond1 | cond2)  {
				resVec[loc] = 0;
			} else {
				double total = 0;
				double elp23_XiVj = 0;
				double elp23_XjVi = 0;
				for (int k = 0; k < nCov; ++k) {
					elp23_XiVj += b23[k] * X23(i,k);
					elp23_XjVi += b23[k] * X23(j,k);
				}
				elp23_XiVj += b23[nCov] * V[j];
				elp23_XjVi += b23[nCov] * V[i];
				elp23_XiVj = exp(elp23_XiVj);
				elp23_XjVi = exp(elp23_XjVi);
				total += (log(elp13[i]) - log(elp13[j])) * (D2[j] - D2[i]);
				total += (H013V[i] - H013V[j]) * (elp13[i] - elp13[j]);
				if (D1[j] & (R[i] > V[j])) {
					total += (H023V[j] - H023R[i]) * elp23_XiVj;
				}
				if (D1[i] & (R[i] > V[i])) {
					total += (H023R[i] - H023V[i]) * elp23[i];
				}
				if (D1[i] & (R[j] > V[i])) {
					total += (H023V[i] - H023R[j]) * elp23_XjVi;
				}
				if (D1[j] & (R[j] > V[j])) {
					total += (H023R[j] - H023V[j]) * elp23[j];
				}
				total = exp(total); 
				resVec[loc] = total; 
			}
			loc += 1;
		}
	}
	return resVec;
}

// [[Rcpp::export]]
// A function for getting the Zeta terms for all subsampled pairs, assuming censoring is random.
NumericVector getZetaTermsC(NumericVector R, NumericVector V, LogicalVector D1, LogicalVector D2, NumericMatrix X23 , NumericVector elp13, NumericVector elp23, NumericVector elpC, NumericVector H013V, NumericVector H023V, NumericVector H023R, NumericVector H0CV, NumericVector H0CR, NumericVector b23, int numPairsPerElement) {
	int n = R.size();
	int nCov = X23.ncol();
	NumericVector resVec(n * numPairsPerElement);
	double loc = 0;
	for (int i = 0; i < n; ++i) {
		for (int p = 1; p <= numPairsPerElement; ++p) {
			int j = (i + p) % n;
			bool cond1 = (V[j] < R[i]) & !D1[j];
			bool cond2 = (V[i] < R[j]) & !D1[i];
			if (cond1 | cond2)  {
				resVec[loc] = 0;
			} else {
				double total = 0;
				double elp23_XiVj = 0;
				double elp23_XjVi = 0;
				for (int k = 0; k < nCov; ++k) {
					elp23_XiVj += b23[k] * X23(i,k);
					elp23_XjVi += b23[k] * X23(j,k);
				}
				elp23_XiVj += b23[nCov] * V[j];
				elp23_XjVi += b23[nCov] * V[i];
				elp23_XiVj = exp(elp23_XiVj);
				elp23_XjVi = exp(elp23_XjVi);
				total += (log(elp13[i]) - log(elp13[j])) * (D2[j] - D2[i]);
				total += (H013V[i] - H013V[j]) * (elp13[i] - elp13[j]);
				if (D1[j] & (R[i] > V[j])) {
					total += (H023V[j] - H023R[i]) * elp23_XiVj;
				}
				if (D1[i] & (R[i] > V[i])) {
					total += (H023R[i] - H023V[i]) * elp23[i];
				}
				if (D1[i] & (R[j] > V[i])) {
					total += (H023V[i] - H023R[j]) * elp23_XjVi;
				}
				if (D1[j] & (R[j] > V[j])) {
					total += (H023R[j] - H023V[j]) * elp23[j];
				}
				if (D1[i] | D2[i]) {
					if (V[i] > R[j]) {
						total += (H0CR[j]-H0CV[i]) * elpC[j];
					}
					if (V[i] > R[i]) {
						total += (H0CV[i] - H0CR[i]) * elpC[i]; 
					}
				} else {
					total += log(elpC[j]) - log(elpC[i]) + H0CV[i] * (elpC[i] - elpC[j]) -H0CR[i] * elpC[i] + H0CR[j] * elpC[j];
				}
				if (D1[j] | D2[j]) {
					if (V[j] > R[i]) {
						total += (H0CR[i] - H0CV[j]) * elpC[i];
					}
					if (V[j] > R[j]) {
						total += (H0CV[j] - H0CR[j]) * elpC[j]; 
					}
				} else {
					total += log(elpC[i]) - log(elpC[j]) + H0CV[j] * (elpC[j] - elpC[i]) - H0CR[j] * elpC[j] + H0CR[i] * elpC[i];
				}
				total = exp(total); 
				resVec[loc] = total; 
			}
			loc += 1;
		}
	}
	return resVec;
}


// [[Rcpp::export]]
// Getting the pairwise pseudolikelihood based on subsampled pairs
double getPairwiseLogLike(NumericVector b12, IntegerVector D1, NumericMatrix X12, NumericVector ZetaTerms, NumericVector H012V, int numPairsPerElement, NumericVector weights) {
	int n = D1.size();
	int nCov = X12.ncol();
	double res = 0;
	double loc = 0;
	NumericVector lp12(n);
	NumericVector elp12(n);
	for (int i = 0; i < n; ++i) {  //Creating elp12 with given b12
		double lp12Cur = 0;
		for (int j = 0; j < nCov; ++j) {
			lp12Cur += b12[j] * X12(i,j); 
		}
		lp12[i] = lp12Cur;
		elp12[i] = exp(lp12Cur);
	}
	for (int i = 0; i < n; ++i) {
		for (int p = 1; p <= numPairsPerElement; ++p) {
			int j = (i + p) % n;
			double A = exp((lp12[i] - lp12[j]) * (D1[j] - D1[i])) * exp((H012V[i] - H012V[j]) * (elp12[i] - elp12[j])) * ZetaTerms[loc];
			res += weights(i) * weights(j) * log(1 + A);
			loc += 1;
		}
	}
	return res;
}

//[[Rcpp::export]]
// Getting the derivative of the pairwise pseudolikelihood based on subsampled pairs
NumericVector getPairwiseDeriv(NumericVector b12, IntegerVector D1, NumericMatrix X12, NumericVector ZetaTerms, NumericVector H012V,  int numPairsPerElement, NumericVector weights) {
	int n = D1.size();
	int nCov = X12.ncol();
	NumericVector resVec(nCov);
	double loc = 0;
	NumericVector lp12(n);
	NumericVector elp12(n);
	for (int i = 0; i < n; ++i) {  //Creating elp12 with given b12
		double lp12Cur = 0;
		for (int j = 0; j < nCov; ++j) {
			lp12Cur += b12[j] * X12(i,j); 
		}
		lp12[i] = lp12Cur;
		elp12[i] = exp(lp12Cur);
	}
	for (int i = 0; i < n; ++i) {
		for (int p = 1; p <= numPairsPerElement; ++p) {
			int j = (i + p) % n;
				double H012VDiff = H012V[i] - H012V[j];
				double A = exp((lp12[i] - lp12[j]) * (D1[j] - D1[i]) + H012VDiff * (elp12[i] - elp12[j])) * ZetaTerms[loc];
				double A_ratio = A / (1 + A);
				for (int k = 0; k < nCov; ++k) {
					resVec[k] += weights(i) * weights(j) * A_ratio * ((X12(i,k) - X12(j,k)) * (D1[j] - D1[i]) + H012VDiff * (elp12[i] * X12(i,k) - elp12[j] * X12(j,k)));
				}
			loc += 1;
		}
	}
	return resVec;
}


NumericVector getPairwiseDerivOnePair(NumericVector b12, IntegerVector D1, NumericVector X1, NumericVector X2, double ZetaTerms, NumericVector H012V) {
	int nCov = X1.size();
	NumericVector resVec(nCov);
	NumericVector lp12(2);
	NumericVector elp12(2);
	double lp12Cur = 0;
	for (int j = 0; j < nCov; ++j) {
		lp12Cur += b12[j] * X1(j); 
	}
	lp12[0] = lp12Cur;
	elp12[0] = exp(lp12Cur);
	lp12Cur = 0;
	for (int j = 0; j < nCov; ++j) {
		lp12Cur += b12[j] * X2(j); 
	}
	lp12[1] = lp12Cur;
	elp12[1] = exp(lp12Cur);
	double H012VDiff = H012V[0] - H012V[1];
	double A = exp((lp12[0] - lp12[1]) * (D1[1] - D1[0]) + H012VDiff * (elp12[0] - elp12[1])) * ZetaTerms;
	double A_ratio = A / (1 + A);
	for (int k = 0; k < nCov; ++k) {
			resVec[k] += A_ratio * ((X1(k) - X2(k)) * (D1[1] - D1[0]) + H012VDiff * (elp12[0] * X1(k) - elp12[1] * X2(k)));
		}
	return resVec;
}


// [[Rcpp::export]]
// Getting the cross-covariance matrix of the pairwise pseudolikelihood score terms.
NumericMatrix getCrossCov(NumericVector b12, IntegerVector D1, NumericMatrix X12, NumericVector ZetaTerms, NumericVector H012V, int numPairsPerElement, int numPairsPerElementVar) {
	int n = D1.size();
	int nCov = X12.ncol();
	NumericMatrix resMat(nCov);
	NumericMatrix resMat1(nCov);
	NumericMatrix resMat2(nCov);
	NumericMatrix PsiMat(numPairsPerElementVar,nCov);
	NumericVector Psi(nCov);
	double loc = 0;
	for (int i = 0; i < n; i++) {
		for (int p = 1; p <= numPairsPerElementVar; ++p) {
			int j = (i + p) % n;
			Psi = getPairwiseDerivOnePair(b12, {D1(i),D1(j)}, X12(i , _ ), X12(j , _ ), ZetaTerms(loc), {H012V(i), H012V(j)}); 
			PsiMat((p - 1) , _ ) = Psi;
			loc += 1;
		}
		for (int k1 = 0; k1 < nCov; k1++) {
			for (int k2 = k1; k2 < nCov; k2++) {
				for (int j = 0; j < numPairsPerElementVar; j++) {
					for (int m = 0; m < numPairsPerElementVar; m++) {
						if (m == j) {
							resMat1(k1,k2) += (PsiMat(j,k1)) * (PsiMat(m,k2)) ;
							continue;
						}
						resMat2(k1,k2) += (PsiMat(j,k1)) * (PsiMat(m,k2)) ;

					}	
				}
			}	
		}	
	}
	resMat += resMat1 * numPairsPerElement / numPairsPerElementVar;
	resMat += resMat2 * 2 * numPairsPerElement * (2*numPairsPerElement - 1) / (numPairsPerElementVar * (numPairsPerElementVar - 1));
	return resMat;
}


// [[Rcpp::export]]
// Getting the Hessian matrix of the pairwise log pseudolikelihood.
NumericMatrix getPairwiseHessian(NumericVector b12, LogicalVector D1, NumericMatrix X12, NumericVector ZetaTerms, NumericVector H012V,  int numPairsPerElement, NumericVector weights) {
	int n = D1.size();
	int nCov = X12.ncol();
	NumericMatrix resMat(nCov);
	double loc = 0;
	NumericVector lp12(n);
	NumericVector elp12(n);
	for (int i = 0; i < n; ++i) {  //Creating elp12 with given b12
		double lp12Cur = 0;
		for (int j = 0; j < nCov; ++j) {
			lp12Cur += b12[j] * X12(i,j); 
		}
		lp12[i] = lp12Cur;
		elp12[i] = exp(lp12Cur);
	}
	for (int i = 0; i < n; ++i) {
		for (int p = 1; p <= numPairsPerElement; ++p) {
			int j = (i + p) % n;
				double H012VDiff = H012V[i] - H012V[j];
				double part1 = exp((lp12(i) - lp12(j)) * (D1(j) - D1(i)) + H012VDiff * (elp12[i] - elp12[j]));
				double A = part1 * ZetaTerms[loc];
				NumericVector Aderiv1(nCov);
				for (int k = 0; k < nCov; ++k) {
					Aderiv1[k] = A * ((X12(i,k) - X12(j,k)) * (D1(j) - D1(i)) + H012VDiff * (X12(i,k) * elp12[i] - X12(j,k) * elp12[j]));
				}
				for (int k = 0; k < nCov; ++k) {
					for(int m = k; m < nCov; ++m) {
						if (A == 0) {
							continue;
						}
						double Aderiv2 = Aderiv1[k] * Aderiv1[m] / A + H012VDiff * A * (X12(i,k) * X12(i,m) * elp12[i] - X12(j,k) * X12(j,m) * elp12[j]);
						resMat(k,m) += weights(i) * weights(j) * (Aderiv2 / (1 + A) - Aderiv1[k] * Aderiv1[m] / pow(1 + A, 2));
					}
				}
			//}
			loc += 1;
		}
	}
	for (int k = 1; k < nCov; ++k) {
		for(int m = 0; m < k; ++m) {
			resMat(k,m) = resMat(m,k);
		}
	}
	return resMat;

}

// [[Rcpp::export]]
// A function that return the Breslow estimator "jumps", allowing usage of weights
// In the following, assume that we join the vectors V and R into a single vector, and then sort it.
// elp = a vector of size 2n, of the exponented linear predictors corresponding to the sorted V and R vector. 
// VorR = indicator vector of size 2n, taking the value 1 if elp[i] corresponds to a V time and 0 if it corresponds to an R time.
// D = a vector of event indicators of size 2n. Should be 0 for R times (VorR=0), and 1 for V times (VorR=1) that correspond to an observed event time. 
NumericVector getBreslow(NumericVector elp, LogicalVector VorR, LogicalVector D, NumericVector weights) {
	int nJumps = 0;
	int eventCounter = 0;
	double elpSum = 0;
	int n = D.size();
	for(int i = 0; i < n; ++i) {
		nJumps += D(i);
	}
	NumericVector res(nJumps);
	for(int i = 0; i < n; ++i)
	{
		if(VorR(i) == 0) {
			elpSum += weights(i) * elp(i);
		} else {
			if(D(i) == 1) {
			res(eventCounter) = weights(i)/elpSum;
			eventCounter += 1;
			}
		elpSum -= weights(i) * elp(i);
		}
	}
	return res;
}

// [[Rcpp::export]]
// Getting Likelihood2
// jumps12 = the jumps of the transition 1->2 Breslow estimator.
// T12sorted = vector of sorted unique event times of transition 1->2.
// H012T,H013T,H023T = the transition-specific Breslow estimators, evaluated at the T12sorted. 
// b23V = the estimated coefficient corresponding to time V when added as a covariate in the PL estimator for transition 2->3.
double getLogLike2(NumericVector b12, IntegerVector D1, NumericMatrix X12, NumericVector R, NumericVector H012V, NumericVector H012R, NumericVector H013R, NumericVector elp13, NumericVector elp23, NumericVector jumps12, NumericVector H012T, NumericVector H013T, NumericVector T12sorted, NumericVector H023R, NumericVector H023T, double b23V) {
	int n = D1.size();
	int n12 = T12sorted.size();
	int nCov = X12.ncol();
	double res = 0;
	NumericVector lp12(n);
	NumericVector elp12(n);
	for (int i = 0; i < n; ++i) {  //Creating elp12 with given b12
		double lp12Cur = 0;
		for (int j = 0; j < nCov; ++j) {
			lp12Cur += b12[j] * X12(i,j); 
		}
		lp12[i] = lp12Cur;
		elp12[i] = exp(lp12Cur);
	}
	for (int i = 0; i < n; ++i) {
		double part1 = D1(i) * lp12(i) - H012V(i) * elp12(i);
		double intsum = 0;
		for(int j = 0; j < n12; ++j) {
			if(T12sorted(j) > R(i)) {
				break;
			}
			intsum += jumps12(j) * elp12(i) * exp(-H012T(j)*elp12(i) - H013T(j)*elp13(i) + (H023T(j) - H023R(i)) * elp23(i)*exp(b23V*T12sorted(j)));
		} 
		double part2 = exp(-H012R(i)*elp12(i) - H013R(i)*elp13(i)) + intsum;
		res += part1 - log(part2);
	}
	return -res;
}

// [[Rcpp::export]]
// Getting the gradient of Likelihood2
NumericVector getLogLike2Deriv(NumericVector b12, IntegerVector D1, NumericMatrix X12,  NumericVector R, NumericVector H012V, NumericVector H012R, NumericVector H013R, NumericVector elp13, NumericVector elp23, NumericVector jumps12, NumericVector H012T, NumericVector H013T, NumericVector T12sorted, NumericVector H023R, NumericVector H023T, double b23V) {
	int n = D1.size();
	int n12 = T12sorted.size();
	int nCov = X12.ncol();
	NumericVector res(nCov);
	NumericVector lp12(n);
	NumericVector elp12(n);
	for (int i = 0; i < n; ++i) {  //Creating elp12 with given b12
		double lp12Cur = 0;
		for (int j = 0; j < nCov; ++j) {
			lp12Cur += b12[j] * X12(i,j); 
		}
		lp12[i] = lp12Cur;
		elp12[i] = exp(lp12Cur);
	}
	for (int i = 0; i < n; ++i) {
		double part1 = D1(i) - H012V(i) * elp12(i);
		double part2 = exp(-H012R(i)*elp12(i) - H013R(i)*elp13(i));
		double tmp = 0;
		double intsum1 = 0;
		double intsum2 = 0;
		for(int j = 0; j < n12; ++j) {
			if(T12sorted(j) > R(i)) {
				break;
			}
			tmp = jumps12(j) * elp12(i) * exp(-H012T(j)*elp12(i) - H013T(j)*elp13(i) + (H023T(j) - H023R(i)) * elp23(i) * exp(b23V*T12sorted(j)));
			intsum1 += tmp;
			intsum2 += tmp * (1 - H012T(j) * elp12(i));
		}
		tmp = part1 - (part2 * (-H012R(i) * elp12(i)) + intsum2) / (part2 + intsum1);
		for(int k = 0; k < nCov; ++k) {
			res(k) += tmp * X12(i,k); 
		}
	}
	return -res;
}