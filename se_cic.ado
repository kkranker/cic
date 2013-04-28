mata:
// Univariate density function
real scalar fden(real scalar y, real colvector Y, | real colvector wgt) {
	// 	this function estimates a univariate density function using kernel methods
	// 	the kernel function is the Epanechnikov kernel
	// 	the bandwidth is the optimal bandwith based on Silverman's rule of thumb
	//
	// 	INPUT
	// 	the input is an N vector of observations Y
	// 	and a scalar y where the density is to be estimated
	//  optionally, provide a frequency weight vector
	//     (for iweghts, wgt must be normalized to the number of observations, so that
	//      sum of weights = number of obserations)
	//
	// 	OUTPUT
	// the output is a scalar with the value of the estimated density

	real scalar h
	real colvector d, kd

	// Silverman optimal bandwidth
	if (args()==2) h=1.06*sqrt(variance(Y))*(length(Y)^(-.2))     // unweighted
	else           h=1.06*sqrt(variance(Y,wgt))*(sum(wgt)^(-.2))  // weighted

	// epanechnikov kernel
	d= abs((Y:-y):/h)
	kd=(d:<sqrt(5)):*(.75:-.15*d:^2)/sqrt(5)

	// return density
	if (args()==2) return(mean(kd/h))     // unweighted
	else           return(mean(kd/h,wgt)) // weighted
}


// test fden w/ this vector
Y= ( 1.11 \ 1.1 \ 1.2 \ 1.3 \ 5 \ 5.1 \ 5.2 \ 10 \ 10 \ 10 )

// the following three should all equal
// 0.072278

fden(1.1,Y)                           // no weights
fden(1.1,Y, J(10,1,1) )               // weight = 1
fden(1.1,Y[1..8], (J(7,1,1) \ 3) )    // weight = 1, except handle three "10" at bottom

// the following two should equal
// 0.076557

fden(1.1,(Y\Y))                       // stack Y on top of itself
fden(1.1,Y, J(10,1,2))                // weight = 2

// CUMULATIVE DISTRIBUTION FUNCTION
real scalar cdf_bar(real scalar y, real vector P, real vector YS)
{
	// given a cumulative distrubtion function (P) over the support points (YS),
	// returns the probability that a random variable
	// is less than a scalar value y
	if      (y<YS[1]+epsilon(y))          return(0)
	else if (y>YS[length(YS)]+epsilon(y)) return(1)
	else                                  return(P[colsum((YS:<(y-epsilon(y))))])
}

// Standard Error for CIC ROUTINE
struct cic_result scalar se_cic(real colvector Y00, real colvector Y01, real colvector Y10, real colvector Y11, real rowvector at, | real colvector W00, real colvector W01, real colvector W10, real colvector W11 )
{
	// Inputs:
	//   (1)-(4) Four column vectors with dependent variable
	//            - Y00 is data for the control group in pre-period
	//            - Y01 is data for the control group in post period
	//            - Y10 is data for the treatment group in post period
	//            - Y11 is data for the treatment group in post period
	//   (5)     Vector with k>=1 quantiles of interest, ranging from 0 to 1
	//   (6)-(9) (Optional) Column with fweights or iweights for Y00, Y01, Y10, and Y11 (respectively)
	//
	// Output: One structure (cic_result) with four row vectors.
	//   Each vector has (1+k) elements. The first element is the mean, followed by k results (one for each quantile in -at-).
	//   (1) result.se_con       = CIC ESTIMATOR WITH CONTINUOUS OUTCOMES, EQUATION 9
	//   (2) result.se_dci       = CIC MODEL WITH DISCRETE OUTCOMES (UNDER THE CONDITIONAL INDEPENDENCE ASSUMPTION), EQUATION 29
	//   (3) result.se_dcilowbnd = LOWER BOUND ESTIMATE OF DISCRETE CIC MODEL (WITHOUT CONDITIONAL INDEPENDENCE), EQUATION 25
	//   (4) result.se_dciuppbnd = UPPER BOUND ESTIMATE OF DISCRETE CIC MODEL (WITHOUT CONDITIONAL INDEPENDENCE), EQUATION 25

	// The code in cic() is somewhat convoluted because I am
	// calculating all four vectors simultaneously.
	// See the NOTE (below, at the bottom of the file)
	// for alternative routines that are more readily
	// accessible. The calculations here lead to
	// slower run-times.

/* input variable 'At' is unused */

	// Need all or none of args (6)-(9)
	if (args()>5 & args()!=9) _error(( "Expected 5 or 9 arguements, but received " + strofreal(args())))

	// Vector with support points for all four groups combined (YS) and the comparison-post group (YS01)
	real colvector YS, YS00, YS01, YS10, YS11
	YS00 = uniqrows(Y00)
	YS01 = uniqrows(Y01)
	YS10 = uniqrows(Y10)
	YS11 = uniqrows(Y11)
	YS   = uniqrows(Y00\Y01\Y10\Y11)
	if (length(YS)<2) _error("The dependent variable is a constant")

	// Vector with CDF functions of the four treat*post groups (F00,F01,F10,F11)
	// CDFs (w/ and w/out weights declared)
	real colvector f00, f01, f10, f11, F00, F01, F10, F11
	if (args()==5) {
		// CDFs without weights
		f00=prob(Y00,YS)
		f01=prob(Y01,YS)
		f10=prob(Y10,YS)
		f11=prob(Y11,YS)
	}
	else {
		// CDFs with weights
		f00=prob(Y00,YS,W00)
		f01=prob(Y01,YS,W01)
		f10=prob(Y10,YS,W10)
		f11=prob(Y11,YS,W11)
	}
	// because of rounding, sum of probabilities might be slightly different than one
	F00=runningsum(f00)
	F01=runningsum(f01)
	F10=runningsum(f10)
	F11=runningsum(f11)
	F00[length(F00)]=1
	F01[length(F01)]=1
	F10[length(F10)]=1
	F11[length(F11)]=1

	// Results will be returned into a structure w/ 4 vectors for con, dci, lower, upper
	struct cic_result scalar result
result=cic(Y00,Y01,Y10,Y11,at)
result.se_con=result.se_dci=result.se_dcilowbnd=result.se_dciuppbnd=J(1,1+length(at),.)

	// A. continuous estimator
	// A.0. preliminaries
	real colvector F00_10, F01invF00_10, f01F01invF00_10, P, PY00, PY01
	real scalar V00, V01, V10, V11, i
	F00_10=F01invF00_10=f01F01invF00_10=J(length(YS10),1,0)
	for(i=1; i<=length(YS10); ++i) {
		F00_10[i]=cdf(YS10[i],F00,YS)
		F01invF00_10[i]=cdfinv(F00_10[i],F01,YS)
		f01F01invF00_10[i]=fden(F01invF00_10[i],Y01)
	}
	// A.1. contribution of Y00
	P=J(length(YS00),1,0)
	for(i=1; i<=length(YS00); ++i) {
		PY00=((YS00[i]:<=YS10)-F00_10):/f01F01invF00_10
		P[i]=quadcross(PY00,select(f10,f10:>epsilon(1)))
	}
	V00=sum(P:^2:*select(f00,f00:>epsilon(1))):/length(Y00)
	// A.2. contribution of Y01
	P=J(length(YS01),1,0)
	for(i=1; i<=length(YS01); ++i) {
		PY01=-((cdf(YS01[i],F01,YS):<=F00_10):-F00_10):/f01F01invF00_10
		P[i]=quadcross(PY01,select(f10,f10:>epsilon(1)))
	}
	V01=sum((P:^2):*select(f01,f01:>epsilon(1))):/length(Y01)
	// A.3. contribution of Y10
	P=F01invF00_10:-quadcross(F01invF00_10,select(f10,f10:>epsilon(1)))
	V10=sum(P:^2:*select(f10,f10:>epsilon(1))):/length(Y10)
	// A.4. contribution of Y11
	P=YS11:-quadcross(YS,f11)
	V11=sum((P:^2):*select(f11,f11:>epsilon(1))):/length(Y11)
	// A.5 final result
	result.se_con[1]=sqrt(V00+V01+V10+V11)
"result.se_con = "; result.se_con

	// B. dci standard error
	// numerical approximation to delta method
	// four parts to variance
	// B.0. setup
	real colvector der00,der01,der10,der11,I00,I01,I10,I11,tI00,tI01,tI10,tI11,f00c,f01c,f10c,f11c,k_bar,t,dy11
	real scalar    delta,max00,max01,max10,max11,der00_c,der01_c,der10_c,der11_c
	der00=der01=der10=der11=J(length(YS),1,0)
	t=(1::length(YS))
	delta=0.0000001

	// B.1. contribution of Y00
	I00=f00:>epsilon(1)
	tI00=select(t,I00)
	max00=max(tI00)
	V00=(diag(f00)-quadcross(f00',f00'))/length(Y00)
	V00[max00,.]=J(1,length(YS),0)
	V00[.,max00]=J(length(YS),1,0)
	for(i=1; i<=(sum(I00)-1); ++i) {
		f00c=f00
		f00c[tI00[i]]=f00c[tI00[i]]+delta
		f00c[max00]  =f00c[max00]-delta
		der00_c=cic_dci(f00c,f01,f10,f11,YS,YS01,.)
		der00[tI00[i]]=(der00_c-result.dci[1])/delta
	}
	V00=quadcross(der00,V00)*der00
	// B.2. Contribution of Y01
	I01=f01:>epsilon(1)
	tI01=select(t,I01)
	max01=max(tI01)
	V01=(diag(f01)-quadcross(f01',f01'))/length(Y01)
	V01[max01,.]=J(1,length(YS),0)
	V01[.,max01]=J(length(YS),1,0)
	for(i=1; i<=sum(I01)-1; ++i) {
		f01c=f01
		f01c[tI01[i]]=f01c[tI01[i]]+delta
		f01c[max01]=f01c[max01]-delta
		der01_c=cic_dci(f00,f01c,f10,f11,YS,YS01,.)
		der01[tI01[i]]=(der01_c-result.dci[1])/delta
	}
	V01=quadcross(der01,V01)*der01
	// B.3. Contribution of Y10
	I10=f10:>epsilon(1)
	tI10=select(t,I10)
	max10=max(tI10)
	V10=(diag(f10)-quadcross(f10',f10'))/length(Y10)
	V10[max10,.]=J(1,length(YS),0)
	V10[.,max10]=J(length(YS),1,0)
	for(i=1; i<=sum(I10)-1; ++i) {
		f10c=f10
		f10c[tI10[i]]=f10c[tI10[i]]+delta
		f10c[max10]=f10c[max10]-delta
		der10_c=cic_dci(f00,f01,f10c,f11,YS,YS01,.)
		der10[tI10[i]]=(der10_c-result.dci[1])/delta
	}
	V10=quadcross(der10,V10)*der10
	// B.4. Contribution of Y11
	I11=f11:>epsilon(1)
	tI11=select(t,I11)
	max11=max(tI11)
	V11=(diag(f11)-quadcross(f11',f11'))/length(Y11)
	V11[max11,.]=J(1,length(YS),0)
	V11[.,max11]=J(length(YS),1,0)
	for(i=1; i<=sum(I11)-1; ++i) {
		f11c=f11
		f11c[tI11[i]]=f11c[tI11[i]]+delta
		f11c[max11]=f11c[max11]-delta
		der11_c=cic_dci(f00,f01,f10,f11c,YS,YS01,.)
		der11[tI11[i]]=(der11_c-result.dci[1])/delta
	}
	V11=quadcross(der11,V11)*der11
	// B.5 components dci variance
	result.se_dci[1]=sqrt(V00+V01+V10+V11)
"result.se_dci = "; result.se_dci

	// C. lower bound standard error
	k_bar=J(length(YS10),1,0)
	for(i=1; i<=length(YS10); ++i) {
		k_bar[i]=cdfinv(cdf(YS10[i],F00,YS),F01,YS)
	}
	k_bar=k_bar:-quadcross(k_bar,select(f10,f10:>epsilon(1)))
	dy11 =YS11 :-quadcross(YS11, select(f11,f11:>epsilon(1)))
	V10=sum((k_bar:^2):*select(f10,f10:>epsilon(1)))/length(Y10)
	V11=sum((dy11 :*dy11) :*select(f11,f11:>epsilon(1)))/length(Y11)
	result.se_dcilowbnd[1] = sqrt(V10+V11)
"result.se_dcilowbnd  = "; result.se_dcilowbnd

	// D. upper bound standard error
	k_bar=J(length(YS10),1,0)
	for(i=1; i<=length(YS10); ++i) {
		k_bar[i]=cdfinv(cdf_bar(YS10[i],F00,YS),F01,YS)
	}
	k_bar=k_bar:-quadcross(k_bar,select(f10,f10:>epsilon(1)))
	dy11 =YS11 :-quadcross(YS11, select(f11,f11:>epsilon(1)))
	V10=sum((k_bar:^2):*select(f10,f10:>epsilon(1)))/length(Y10)
	V11=sum((dy11:^2) :*select(f11,f11:>epsilon(1)))/length(Y11)
	result.se_dciuppbnd[1] = sqrt(V10+V11)
"result.se_dciuppbnd  = "; result.se_dciuppbnd

	// DONE.  RETURN STRUCTURE W/ ROW VECTORS CONTAINING POINT ESTIMATES.
	return(result)
} // end of cic

end // exit Mata
