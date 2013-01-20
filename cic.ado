clear all

*! $Id$
*! Changes-in-changes
*!
*! An implimentation of:
*! Athey, S. & G. W. Imbens. "Identification and Inference in Nonlinear Difference-in-Differences Models."
*!    	Econometrica, 74 (2), March 2006, pp. 431-497.
*!
*! Stata code by Keith Kranker
*! Based on Matlab code by S. Athey & G. W. Imbens, published on S. Athey's website
*! $Date$
*
*  Syntax:
*
*  cic y_var treat_var post_var [control_varlist] [if] [in] [using] [weight] [, options]
*
*  where
*    cmd                      is a regression command (e.g,. regress, logit)
*    y_var                    is the dependent variable
*    treat_var                is a dummy that equals 1 for the treatment group
*    post_var                 is a dummy that equals 1 for the period in which the treatment group is treated
*    control_varlist          is a list of control variables (optional)
*
*  and options are
*    at(numlist)              a list of percentiles for CIC results. default is at(10(10)90)
*    vce(none|                don't calculate standard errors, the default
*        delta|               use numerical
*        bootstrap[, bsopts]| use bootstrap (by default, 1000 reps stratified by treat/post) other options allowed
*        bspctile[, bsopts])  use Athey and Imben's of obtaining bootstrap standard errors ( se = (p(97.5) - p(2.5)) / (2*1.96) where p(N) is Nth percentile of bootstrap iterations
*                                 instead of using Stata's default method. The Stata manual ([R] bootstrap) addresses this issue.
*                                    "A better estimate is needed if you want to use the 2.5th and 97.5th percentiles of the distribution
*                                    to produce a 95% confidence interval. To extract many features simultaneously about the distribution,
*                                    an even better estimate is needed. Generally, replications on the order of 1,000 produce very good
*                                    estimates, but only 50�200 replications are needed for estimates of standard errors...."
*    untreated                counterfactual effect of the policy for the untreated group (Setion 3.2 of paper)
*    *                        any other options are passed through to the regression command (if no control varlist)

* weights not allowed with vce(boostrap)

program define cic, eclass byable(recall)
	version 11.2

	// parse arguments
	syntax varlist(default=none min=3 numeric fv ts) [if] [in] [fw iw]  ///
		[, at(numlist min=1 >=0 <=100 sort) ///
		Vce(passthru) ///
		UNTreated ///
		* ]
	marksample touse

	// first three variables need to be y, treat, and post
	gettoken y     varlist : varlist
	gettoken treat varlist : varlist
	gettoken post  varlist : varlist

	// prep to handle weights
	if !missing("`weight'") {
		tempvar wgtvar
		gen `wgtvar' `exp' if `touse'
		local wtexp_caller = `", "`wgtvar'" "'
		summ `wgtvar' if `touse' , meanonly
		local n=r(sum)
		qui count if `touse'
		local n_obs=r(N)
	}
	else {
		qui count if `touse'
		local n=r(N)
	}

	// parse percentiles
	if mi("`at'") local at "10(10)90" // default set (if undeclared)
	numlist "`at'"
	local at = r(numlist)

	// parse vce()
	// if vce(bootstrap), run bootstrap iterations with recursive cic call
	// this section is similar in function to the "_vce_parserun" command, except that I set default values for reps() and strata()
	if !mi("`vce'")  {
		cic_vce_parse `treat' `post' [`weight'`exp'], `vce'
		local vce = r(vce)
		local usebspctile = r(usebspctile)
		if inlist( "`vce'", "bootstrap", "bspctile" ) {
			// bootstrap the cic command (no SEs)
return list
			bootstrap _b, reps( `=r(reps)' ) `r(bsopts)' : cic `y' `treat' `post' `varlist' if `touse' [`weight'`exp'], at(`at')
return list
			if `usebspctile' {
			  if `=r(reps)' < 1000 nois di as error "Warning: More bootstrap repetitions might be needed with vce(bspctile)."
			  tempname ci_se
			  mata : bs_se( "e(ci_percentile)", "`ci_se'" )
			  ereturn repost V = `ci_se'
			  ereturn display
			}
			exit
		}
	}
	else local vce = "none"

	* adjust y for covariates using OLS regression
	if `: list sizeof varlist'!=0 {
		di as txt "Regression to adjust `y' for control variables:"
		regress `y' ib0.`treat'##ib0.`post' `varlist' if `touse' [`weight'`exp']
		tempvar yresid yadj
		qui predict `yresid' if `touse', residuals
		qui gen `yadj' = `yresid' + _b[_cons] + _b[1.`treat']*`treat' + _b[1.`post']*`post' + _b[1.`treat'#1.`post']*`treat'*`post' if `touse'
		label var `yresid' "`y' residuals after adjusting for covariates"
		label var `yadj' "`y' adjusted for covariates"
desc `y' `yadj'
summ `y' `yadj' if `touse'
	}
	else local yadj `y'

	// implement mata CIC routine
	if "`untreated'"=="" {
		// default - effect of treatment on the treated
		mata: cic_caller( "`y'", "`treat'", "`post'", "`touse'", "at", "`vce'" `wtexp_caller')
	}
	else {
		// option - effect of treatment on the untreated
		// just switch `treat' variable and use -e(b)
		tempvar untreated
		gen `untreated' = !`treat'
		mata: cic_caller( "`y'", "`untreated'", "`post'", "`touse'", "at", "`vce'" `wtexp_caller')
		matrix `cic_estimates' = -`cic_estimates'
	}

	// post results to ereturn
	ereturn clear
	ereturn post `cic_estimates' [`weight'`exp'], depname("`y'") obs(`n') esample(`touse')
	ereturn display
	* ereturn local cmd = "cic"
	ereturn local cmdline civ `0'
	if !missing("`weight'") ereturn scalar n_obs=`n_obs'
end // end of cic program definition


// subroutine to parse the vce() option
// this section is similar in function to the "_vce_parse" command, except that I set default values for reps() and strata()
program define cic_vce_parse, rclass
	syntax [varlist] [iw fw], vce(string asis)
	_parse comma vce 0 : vce
	if  inlist( "`vce'","bootstra","bootstr","bootst","boots","boot","boo","bo")     local vce "bootstrap"
	if  inlist( "`vce'","bspctil","bspcti","bspct","bspc","bsp","bsp","bs")          local vce "bspctile"
	if  inlist( "`vce'","delt","del","de","d")                                       local vce "delta"
	if  inlist( "`vce'","non","no","n")                                              local vce "none"
	if !inlist( "`vce'","delta","bootstrap","bspctile","none") {
		di as error "Only vce(delta), vce(bootstrap [, subopts]), vce(bspctile [, subopts]), and vce(none) allowed."
		error 198
	}
	return local vce `vce'
	if inlist("`vce'","bootstrap","bspctile") {
		local stratvars : copy local varlist
		syntax [, Reps(integer 1000) strata(string asis) notable *]
		if mi("`strata'") local strata : copy local stratvars
		return local  usebspctile=( "bspctile"=="`vce'" )
		return local  bsopts  strata(`strata') `=cond("bspctile"=="`vce'","notable","`table'")' `options'
		return scalar reps  = `reps'
		if !missing("`weight'") {
			di as error "Weights not allowed with vce(bootstrap)"
			error 198
		}
	}
	else if !mi("`0'") {
		di as error "suboptions are not allowed with vce(`vce')"
		error 198
	}
end


/* * * * *  BEGIN MATA BLOCK * * * * */
version 11.2
mata:
mata clear
mata set matastrict on
mata set matafavor speed


mata set matalnum on /* drop this later */

// STRUCTURE FOR RETURNING RESULTS
struct cic_result {
	real rowvector con, dci, dcilowbnd, dciuppbnd
}

// CIC CALLER -- THIS FUNCTION READS STATA DATA INTO MATA AND CALLS THE MAIN CIC ROUTINE
void cic_caller(string scalar y_var, string scalar treat_var, string scalar post_var, string scalar touse_var, string scalar at_local, string scalar vce, |string scalar wgt_var)
{
	// Inputs: Names of variables `y', `treat', `post', and `touse'
	//         Name of local macro containing quantiles of interest, ranging from 0 to 100
	//         Dummy indicating whether to calculation of SE is needed
	//         Name of variable with fweight or iweight (optional)
	// Output: Results are returned to a Stata matrix named `cic_estimates'

	// get data into mata
	// rows with missing data are already dropped by -marksample- in .ado file
	real matrix dta
	st_view(dta=., ., (y_var, treat_var, post_var), touse_var)
	if (args()==7) {
		// with weight variable
		real colvector wgt
		st_view(wgt=.,.,wgt_var,touse_var)
	}

	// quantiles
	real rowvector at
	at = strtoreal(tokens(st_local(at_local))) / 100

	// results will be returned into a structure w/ 4 vectors con, dci, lower, upper
	struct cic_result scalar result

	// call the main CIC routine
	if (args()==7) result = cic(dta,at,vce,wgt)  // with weights
	else           result = cic(dta,at,vce)       // without weights

	// save results into a Stata matrix (tempname `cic_estimates')
	stata("tempname cic_estimates")
	st_matrix(st_local("cic_estimates"), (result.con,result.dci,result.dcilowbnd,result.dciuppbnd))
	st_matrixcolstripe(st_local("cic_estimates"), ((J(1+length(at),1,"continuous") \ J(1+length(at),1,"discrete_ci") \ J(1+length(at),1,"dci_lower_bnd") \ J(1+length(at),1,"dci_upper_bnd")),J(4,1,strtoname(("mean" , ("p":+strofreal(at*100))))')))
}


// CIC ROUTINE
struct cic_result scalar cic(real matrix dta, real vector at, string scalar vce, |string scalar wgt)
{
	// Inputs: Matrix with data.
	//           - First column is dependent variable (Y),
	//           - Second column is a control (=0) vs. treatment (=1) dummy
	//           - Third column is a pre-period (=0) vs. post-period (=1) dummy
	//           - (Optional) Fourth column with fweights iweights
	//             (Rows are dropped if the two dummies do not equal 0 or 1)
	//         Vector with k>=1 quantiles of interest (at), ranging from 0 to 1
	// Output: Structure with four vectors.
    //         Each vector has (1+k) elements. The first element is the mean, followed by k results (one for each quantile in -at-).

	// Declare variables
	real colvector y, treat, post
	real colvector Y00, Y01, Y10, Y11
	real colvector YS, YS01
	real colvector F00, F01, F10, F11
	struct cic_result scalar result

	// Subviews with three columns, for convenience
	st_subview(y=.    ,dta,.,1)
	st_subview(treat=.,dta,.,2)
	st_subview(post=. ,dta,.,3)

	// Select the rows belonging to the treat*post groups
	st_select(Y00=.,y,(treat:==0 :& post:==0))
	st_select(Y01=.,y,(treat:==0 :& post:==1))
	st_select(Y10=.,y,(treat:==1 :& post:==0))
	st_select(Y11=.,y,(treat:==1 :& post:==1))

	// Check inputs
	if (!(cols(dta)==3|cols(dta)==4))                              _error( "Expecting input dataset with 3 or 4 variables (y, treat, post [, weight]")
	if (min((length(Y00),length(Y01),length(Y10),length(Y11)))==0) _error( "One or more of the four treat*post groups is empty." )
	if (min(at)<0 | max(at)>1)                                     _error( "at() must be between zero and one (inclusive)." )

	// Vector with support points for all four groups combined (YS) and the comparison-post group (YS01)
	YS   = uniqrows( y )
	YS01 = uniqrows( Y01 )

	// Vector with CDF functions of the four treat*post groups (F00,F01,F10,F11)
	// CDFs (w/ and w/out weights declared)
	if (cols(dta)==3) {
		// without weights
		F00=runningsum(prob(Y00,YS))
		F01=runningsum(prob(Y01,YS))
		F10=runningsum(prob(Y10,YS))
		F11=runningsum(prob(Y11,YS))
	}
	else if (cols(dta)==4) {
		// with weights declared (4th column of dta)
		F00=runningsum(prob(Y00,YS,select(dta[.,4],(treat:==0 :& post:==0))))
		F01=runningsum(prob(Y01,YS,select(dta[.,4],(treat:==0 :& post:==1))))
		F10=runningsum(prob(Y10,YS,select(dta[.,4],(treat:==1 :& post:==0))))
		F11=runningsum(prob(Y11,YS,select(dta[.,4],(treat:==1 :& post:==1))))
		F00[length(F00)]=1 // because of rounding, sum of weights might be slightly different than one
		F01[length(F01)]=1
		F10[length(F10)]=1
		F11[length(F11)]=1
	}

	// CIC ESTIMATOR WITH CONTINUOUS OUTCOMES, EQUATION 9
	result.con = cic_con(F00,F01,F10,F11,YS,YS01,at)

	// CIC MODEL WITH DISCRETE OUTCOMES (UNDER THE CONDITIONAL INDEPENDENCE ASSUMPTION), EQUATION 29
	result.dci = cic_dci(F00,F01,F10,F11,YS,YS01,at)

	// LOWER BOUND ESTIMATE OF DISCRETE CIC MODEL (WITHOUT CONDITIONAL INDEPENDENCE), EQUATION 25
	result.dcilowbnd = cic_lower(F00,F01,F10,F11,YS,YS01,at)

	// UPPER BOUND ESTIMATE OF DISCRETE CIC MODEL (WITHOUT CONDITIONAL INDEPENDENCE), EQUATION 25
	result.dciuppbnd = cic_upper(F00,F01,F10,F11,YS,YS01,at)
	
	return(result)
} // end of cic


// CIC ESTIMATOR WITH CONTINUOUS OUTCOMES, EQUATION 9
real vector cic_con(real vector F00, real vector F01, real vector F10, real vector F11, real vector YS, real vector YS01, real vector at)
{
	// this function calculates the continuous outcomes CIC estimator
	// first estimate the cdf of Y^N_11 using equation (9) in the paper and
	// then use that to calculate the average effect of the treatment
	real vector FCO, est_con
	real scalar i, F01y, F00invF01y, F10F00invF01y

	// for each y in the support of Y01, fill in FCO(y)=F_10(F^-1_00(F_01(y)))
	FCO=J(length(YS01),1,0)
	for(i=1; i<=length(YS01); ++i) {
		F01y=cdf(YS01[i],F01,YS)
		F00invF01y=cdfinv(F01y,F00,YS)
		F10F00invF01y=cdf(F00invF01y,F10,YS)
		FCO[i] = F10F00invF01y
	}
	FCO[length(FCO)]=1 // =1 at end

	// mean CIC estimate
	est_con=( (F11-(0 \ F11[1..(length(YS)-1)]))'*YS - (FCO-(0 \ FCO[1..(length(YS01)-1)]))'*YS01 )

	// quantile CIC estimate
	for(i=1; i<=length(at); ++i) {
		est_con = (est_con , ( cdfinv(at[i], F11, YS) - cdfinv(at[i], FCO, YS01) ) )
	}

	// matrix has mean estimate in first column, plus one column for each element of "at"
	return(est_con)
}


// CIC MODEL WITH DISCRETE OUTCOMES (UNDER THE CONDITIONAL INDEPENDENCE ASSUMPTION), EQUATION 29
real vector cic_dci(real vector F00, real vector F01, real vector F10, real vector F11, real vector YS, real vector YS01, real vector at)
{
	// this function calculates the discreate outcomes CIC estimator
	// first estimate the cdf of Y^N_11 using equation (29) in the paper and
	// then use that to calculate the average effect of the treatment
	real vector FDCI, FUB, FLB, est_dci
	real scalar i,F01y,F00invF01y,F10F00invF01y,F00invbF01y,F10F00invbF01y,F00F00invF01y,F00F00invbF01y,FDCI_weight

	// for each y in the support of Y01, fill in FCO(y)=F_10(F^-1_00(F_01(y)))
	FDCI=FUB=FLB=J(length(YS01),1,0)
	for(i=1; i<=length(YS01); ++i) {
		F01y=cdf(YS01[i],F01,YS)
		F00invF01y=cdfinv(F01y,F00,YS)
		F10F00invF01y=cdf(F00invF01y,F10,YS)
		F00invbF01y=cdfinv_brckt(F01y,F00,YS)
		F10F00invbF01y=cdf(F00invbF01y,F10,YS)
		F00F00invF01y =cdf(F00invF01y ,F00,YS)
		F00F00invbF01y=cdf(F00invbF01y,F00,YS)
		FLB[i]=F10F00invbF01y
		FUB[i]=F10F00invF01y
		if ((F00F00invF01y-F00F00invbF01y)>epsilon(1)) FDCI_weight=(F01y-F00F00invbF01y)/(F00F00invF01y-F00F00invbF01y)
		else                                           FDCI_weight=0
		FDCI[i]=FLB[i]+(FUB[i]-FLB[i])*FDCI_weight
	}
	FDCI[length(FDCI)]=1 // =1 at end

	// conditional independence estimate
	est_dci=( (F11-(0 \ F11[1..(length(YS)-1)]))'*YS - (FDCI-(0 \ FDCI[1..(length(YS01)-1)]))'*YS01 )


	// quantile CIC estimate
	for(i=1; i<=length(at); ++i) {
		est_dci = (est_dci , ( cdfinv(at[i], F11, YS) - cdfinv(at[i], FDCI, YS01) ) )
	}

	// matrix has mean estimate in first column, plus one column for each element of "at"
	return(est_dci)
}


// LOWER BOUND ESTIMATE OF DISCRETE CIC MODEL (WITHOUT CONDITIONAL INDEPENDENCE), EQUATION 25
real vector cic_lower(real vector F00, real vector F01, real vector F10, real vector F11, real vector YS, real vector YS01, real vector at)
{
	// this function calculates the discreate outcomes CIC estimator
	// first estimate the cdf of Y^N_11 using equation (29) in the paper and
	// then use that to calculate the average effect of the treatment
	real vector FLB, est_lower
	real scalar i,F01y, F00invbF01y,F10F00invbF01y

	// for each y in the support of Y01, fill in FCO(y)=F_10(F^-1_00(F_01(y)))
	FLB=J(length(YS01),1,0)
	for(i=1; i<=length(YS01); ++i) {
		F01y=cdf(YS01[i],F01,YS)
		F00invbF01y=cdfinv_brckt(F01y,F00,YS)
		F10F00invbF01y=cdf(F00invbF01y,F10,YS)
		FLB[i]=F10F00invbF01y;
	}
	FLB[length(FLB)]=1 // =1 at end

	// conditional independence estimate
	est_lower=( (F11-(0 \ F11[1..(length(YS)-1)]))'*YS - (FLB-(0 \ FLB[1..(length(YS01)-1)]))'*YS01 )

	// quantile CIC estimate
	for(i=1; i<=length(at); ++i) {
		est_lower = (est_lower , ( cdfinv(at[i], F11, YS) - cdfinv(at[i], FLB, YS01) ) )
	}

	// matrix has mean estimate in first column, plus one column for each element of "at"
	return(est_lower)
}


// UPPER BOUND ESTIMATE OF DISCRETE CIC MODEL (WITHOUT CONDITIONAL INDEPENDENCE), EQUATION 25
real vector cic_upper(real vector F00, real vector F01, real vector F10, real vector F11, real vector YS, real vector YS01, real vector at)
{
	// this function calculates the discreate outcomes CIC estimator
	// first estimate the cdf of Y^N_11 using equation (29) in the paper and
	// then use that to calculate the average effect of the treatment
	real vector FUB, est_upper
	real scalar i,F01y,F00invF01y,F10F00invF01y

	// for each y in the support of Y01, fill in FCO(y)=F_10(F^-1_00(F_01(y)))
	FUB=J(length(YS01),1,0)
	for(i=1; i<=length(YS01); ++i) {
		F01y=cdf(YS01[i],F01,YS)
		F00invF01y=cdfinv(F01y,F00,YS)
		F10F00invF01y=cdf(F00invF01y,F10,YS)
		FUB[i]=F10F00invF01y;
	}
	FUB[length(FUB)]=1 // =1 at end

	// conditional independence estimate
	est_upper=( (F11-(0 \ F11[1..(length(YS)-1)]))'*YS - (FUB-(0 \ FUB[1..(length(YS01)-1)]))'*YS01 )

	// quantile CIC estimate
	for(i=1; i<=length(at); ++i) {
		est_upper = (est_upper , ( cdfinv(at[i], F11, YS) - cdfinv(at[i], FUB, YS01) ) )
	}

	// matrix has mean estimate in first column, plus one column for each element of "at"
	return(est_upper)
}


// SAMPLE PROPORTIONS
real vector prob(real vector Y, real vector YS, |real vector wgt)
{
	// given a vector Y and a vector of support points YS
	// this function calculates the sample proportions at each of the support points
	// wgt is an (optional) set of weights for the vector Y
	real scalar n
	n = length(YS)
	if (args()==3) {
		// with weight variable
		return(rowsum((abs((YS:-J(n,1,Y'))):<=epsilon(1)):*J(n,1,wgt')):/J(n,1,colsum(wgt)))
	}
	else {
		// with equal weights
		return(rowsum(abs((YS:-J(n,1,Y'))):<=epsilon(1)):/length(Y))
	}
}


// CUMULATIVE DISTRIBUTION FUNCTION
real scalar cdf(real scalar y, real vector P, real vector YS)
{
	// given a cumulative distrubtion function (P) over the support points (YS),
	// returns the empirical cumulative distribution function at a scalar (y)
	if      (y< min(YS)) return(0)
	else if (y>=max(YS)) return(1)
	else                 return(P[colsum((YS:<=(y+epsilon(1))))])
}


// INVERSE OF CUMULATIVE DISTRIBUTION FUNCTION, EQUATION 8
real scalar cdfinv(real scalar p, real vector P, real vector YS)
{
	// given a cumulative distrubtion functin (P) over the support points (YS),
	// returns the inverse of the empirical cumulative distribution function at probability p (0<p<1)
	return(YS[min((length(YS)\select((1::length(YS)),(P:>=(p-epsilon(1))))))])
}


// INVERSE OF CUMULATIVE DISTRIBUTION FUNCTION, ALTERNATIVE FOR DISCREATE OUTCOMES, EQUATION 24
real scalar cdfinv_brckt(real scalar p, real vector P, real vector YS)
{
	// given a cumulative distrubtion functin (P) over the support points (YS),
	// returns the inverse of the empirical cumulative distribution function at probability p (0<p<1)
	// but if equals -oo, it returns min(YS)-100*(1+max(YS)-YS(min)) = 101*min(YS)-100*max(YS)-100
	if (p>=(P[1]-epsilon(1))) {
		return(YS[max(select((1::length(YS)),(P:<=(p+epsilon(1)))))])
	}
	else {
		return(101*YS[1]-100*YS[length(YS)]-100)
	}
}


// USE 95 PERCENTILES OF BOOTSTRAP ITERATIONS TO BACKOUT STANDARD ERRORS
void bs_se( string scalar in_ci, string scalar out_V )
{
	real vector bs_se
	bs_se = (1 /(2 * 1.96)) * (st_matrix(in_ci)[2,.] - st_matrix(in_ci)[1,.])
	st_matrix(out_V,diag(bs_se:*bs_se))
	st_matrixcolstripe(out_V,st_matrixcolstripe(in_ci))
	st_matrixrowstripe(out_V,st_matrixcolstripe(in_ci))
out_V
}

"check prob"
prob((1\1\2\3),(1\2\3\4\5))
prob((1\1\2\3),(1\2\3\4\5),(1\1\2\1))
prob((1\  2\3),(1\2\3\4\5),(2  \2\1))
prob((1\  2\3),(1\2\3\4\5),(1.5\1.5\.75))

"check cdfinv & cdfinv_brckt"
P  = (.1\.3\.6\.7\1.0)
YS = (1\2\3\4\5)
/* 1    = */ cdfinv(.05  , P, YS)
/* 1    = */ cdfinv(.1   , P, YS)
/* 2    = */ cdfinv(.2   , P, YS)
/* 2    = */ cdfinv(.2999, P, YS)
/* 2    = */ cdfinv(.3   , P, YS)
/* 3    = */ cdfinv(.3001, P, YS)
/* 5    = */ cdfinv(.9999, P, YS)
/* 5    = */ cdfinv(1    , P, YS)

/* -499 = */ cdfinv_brckt(.05  , P, YS)
/* 1    = */ cdfinv_brckt(.1   , P, YS)
/* 1    = */ cdfinv_brckt(.2   , P, YS)
/* 1    = */ cdfinv_brckt(.2999, P, YS)
/* 2    = */ cdfinv_brckt(.3   , P, YS)
/* 2    = */ cdfinv_brckt(.3001, P, YS)
/* 4    = */ cdfinv_brckt(.9999, P, YS)
/* 5    = */ cdfinv_brckt(1    , P, YS)


end
/* * * * *  END OF MATA BLOCK * * * * */




* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* test cic_vce_parse
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
cic_vce_parse, vce(boot, reps(1000) saving(myfile.dta, replace))
return list
cic_vce_parse, vce(boot, reps(25) strata(x y) cluster(treat))
return list
cic_vce_parse, vce(boot, reps(25) strata() cluster())
return list


* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
set tracedepth 2
if 0  set trace on
else  set trace off
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
local Nreps = 500
if 00     	local vce vce(bootstrap, reps(`Nreps'))
else if 0	local vce vce(bspctile , reps(`Nreps'))
else       	local vce ""
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* with athey and imbens data  suppliment
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

// load data
qui {
	infix y 1-4 after 7-8 high 9-10 male 11-12 marital 13-14 manu 17-18 ///
		constr 19-20 head 21-22  neck 23-24 upper_extr 25-26   trunk 27-28  ///
		low_back 29-30 lower_extr 31-32 occ_dis 33-34 state 39-41 v18 37-38 ///
		age 42-44 prev_earn 45-60 ///
		using "C:\Projects\Americhoice_Data\Matlab_CIC_Routine\mvd.dat"

	gen ind=(age<99) & (marital<8) & (male<9) & (manu<8) & (constr<8) & (v18<9)
	replace head=head==1
	replace neck=neck==1
	replace upper_extr=upper_extr==1
	replace trunk=trunk==1
	replace low_back=low_back==1
	replace lower_extr=lower_extr==1
	replace occ_dis=occ_dis==1
	gen ky=state==16
	keep if ky

	replace y=.25 if (y==0)  // add 0.25 to zero durations before taking logarithms
	gen ly=log(y)

	gen     high_after = .
	replace high_after = 1 if high==0 & after==0
	replace high_after = 2 if high==0 & after==1
	replace high_after = 3 if high==1 & after==0
	replace high_after = 4 if high==1 & after==1
	label define high_after 1 "Control Group, First Period"   ///
	                        2 "Control Group, Second Period"  ///
	                        3 "Treatment Group, First Period" ///
	                        4 "Treatment Group, Second Period"
	label val high_after high_after
}
set seed 1

// check tabulations
count
tabstat y ly , by(high_after) s(count mean sd min p25 p50 p75 p90 max) columns(s)  labelwidth(30) nototal format(%9.2f)

// DID estimate
reg y high##after
reg ly high##after

// CIC estimates from A&I Appendix
timer on 1
cic  y high after ,  at(25 50 75 90) `vce'
cic ly high after ,  at(50)          `vce'
cic  y high after ,  at(25 50 75 90) `vce' untreated
cic ly high after ,  at(50)          `vce' untreated
timer off 1
timer list
exit

// compare vce() option above to the bootstrap prefix
bootstrap , reps(`Nreps') strata(high after) : cic y  high after ,  at(25 50 75 90)

// check weights are working
gen testweight =(uniform()<.95) + (uniform()<.20)
tab testweight
cic y  high after [fw=testweight],  at(25 50 75 90)
drop if testweight==0
expand testweight
cic y  high after                ,  at(25 50 75 90)
cic y  high after                ,  at(25 50 75 90)

/*
// direct comparision of vce() options
set seed 1
cic  y high after ,  at(25 50 75 90) vce(bootstrap, reps(`Nreps'))
set seed 1
cic  y high after ,  at(25 50 75 90) vce(bspctile , reps(`Nreps'))
*/

// we should get an error if try to use weights
cap nois {
  cic y  high after [fw=testweight],  at(25 50 75 90) vce(bootstrap, reps(25) nodots)
}

// we should get an error with vce(delta) since I haven't coded it up yet
cap nois cic  y high after ,  at(25 50 75 90) vce(delta)


ereturn list
exit


* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* The following code can be used to test the program using "fake" data
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
sysuse nlsw88, clear
set seed 1
gen TREAT1 = uniform() < .5
replace wage = wage + TREAT1
gen POST1 = uniform() < .5
replace wage = wage - POST1

// bootstrap the sample conditional on Ngt for g; t = 0; 1
cic wage TREAT1 POST1 i.occupation,  at(10(10)90 99.5)



* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* The following code can be used to test the program using "fake" data
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
clear
set obs 4
gen post  = inlist(_n,1,3)
gen treat = inlist(_n,1,2)
local n_g=1000
expand `n_g'
bys p t: gen y = _n / `n_g'
gen     d = 1.75 - 1.5 * y if t==0 & p==0
replace d = 0.75 - 0.5 * y if t==0 & p==1
replace d = 0.80 - 0.4 * y if t==1 & p==0
replace d = 0.50 - 1.0 * y if t==1 & p==1

cic y treat post,  at(10(10)90)

replace y = ceil(y*10)/10
tab y
cic y treat post,  at(10(10)90)

replace y = ceil(y/50)*50
tab y
cic y treat post,  at(10(10)90)
