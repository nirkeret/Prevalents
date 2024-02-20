#include <Rcpp.h>
//#include<iostream>
//using namespace std;
using namespace Rcpp;

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))


// [[Rcpp::export]]
// A function for getting the Zeta terms for all subsampled pairs, assuming censoring is fixed
NumericVector getZetaTermsNoC(NumericVector R, NumericVector V, LogicalVector D1, LogicalVector D2, NumericMatrix X , NumericVector elp13, NumericVector elp23, NumericVector H013V, NumericVector H023V, NumericVector H023R, NumericVector b23, int numPairsPerElement) {
	int n = R.size();
	int nCov = X.ncol();
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
					elp23_XiVj += b23[k] * X(i,k);
					elp23_XjVi += b23[k] * X(j,k);
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
// A function for getting the Zeta terms for all subsampled pairs, assuming censoring is random
NumericVector getZetaTerms(NumericVector R, NumericVector V, LogicalVector D1, LogicalVector D2, NumericMatrix X , NumericVector elp13, NumericVector elp23, NumericVector elpC, NumericVector H013V, NumericVector H023V, NumericVector H023R, NumericVector H0CV, NumericVector H0CR, NumericVector b23, int numPairsPerElement) {
	int n = R.size();
	int nCov = X.ncol();
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
					elp23_XiVj += b23[k] * X(i,k);
					elp23_XjVi += b23[k] * X(j,k);
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
double getPairwiseLogLike(NumericVector b12, IntegerVector D1, NumericMatrix X, NumericVector ZetaTerms, NumericVector H012V, int numPairsPerElement, NumericVector weights) {
	int n = D1.size();
	int nCov = X.ncol();
	double res = 0;
	double loc = 0;
	NumericVector lp12(n);
	NumericVector elp12(n);
	for (int i = 0; i < n; ++i) {  //Creating elp12 with given b12
		double lp12Cur = 0;
		for (int j = 0; j < nCov; ++j) {
			lp12Cur += b12[j] * X(i,j); 
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
NumericVector getPairwiseDeriv(NumericVector b12, IntegerVector D1, NumericMatrix X, NumericVector ZetaTerms, NumericVector H012V,  int numPairsPerElement, NumericVector weights) {
	int n = D1.size();
	int nCov = X.ncol();
	NumericVector resVec(nCov);
	double loc = 0;
	NumericVector lp12(n);
	NumericVector elp12(n);
	for (int i = 0; i < n; ++i) {  //Creating elp12 with given b12
		double lp12Cur = 0;
		for (int j = 0; j < nCov; ++j) {
			lp12Cur += b12[j] * X(i,j); 
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
					resVec[k] += weights(i) * weights(j) * A_ratio * ((X(i,k) - X(j,k)) * (D1[j] - D1[i]) + H012VDiff * (elp12[i] * X(i,k) - elp12[j] * X(j,k)));
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
NumericMatrix getCrossCov(NumericVector b12, IntegerVector D1, NumericMatrix X, NumericVector ZetaTerms, NumericVector H012V, int numPairsPerElement, int numPairsPerElementVar) {
	int n = D1.size();
	int nCov = X.ncol();
	NumericMatrix resMat(nCov);
	NumericMatrix resMat1(nCov);
	NumericMatrix resMat2(nCov);
	NumericMatrix PsiMat(numPairsPerElementVar,nCov);
	NumericVector Psi(nCov);
	double loc = 0;
	for (int i = 0; i < n; i++) {
		for (int p = 1; p <= numPairsPerElementVar; ++p) {
			int j = (i + p) % n;
			Psi = getPairwiseDerivOnePair(b12, {D1(i),D1(j)}, X(i , _ ), X(j , _ ), ZetaTerms(loc), {H012V(i), H012V(j)}); 
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
NumericMatrix getPairwiseHessian(NumericVector b12, LogicalVector D1, NumericMatrix X, NumericVector ZetaTerms, NumericVector H012V,  int numPairsPerElement, NumericVector weights) {
	int n = D1.size();
	int nCov = X.ncol();
	NumericMatrix resMat(nCov);
	double loc = 0;
	NumericVector lp12(n);
	NumericVector elp12(n);
	for (int i = 0; i < n; ++i) {  //Creating elp12 with given b12
		double lp12Cur = 0;
		for (int j = 0; j < nCov; ++j) {
			lp12Cur += b12[j] * X(i,j); 
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
					Aderiv1[k] = A * ((X(i,k) - X(j,k)) * (D1(j) - D1(i)) + H012VDiff * (X(i,k) * elp12[i] - X(j,k) * elp12[j]));
				}
				for (int k = 0; k < nCov; ++k) {
					for(int m = k; m < nCov; ++m) {
						if (A == 0) {
							continue;
						}
						double Aderiv2 = Aderiv1[k] * Aderiv1[m] / A + H012VDiff * A * (X(i,k) * X(i,m) * elp12[i] - X(j,k) * X(j,m) * elp12[j]);
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