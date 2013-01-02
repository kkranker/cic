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
*    cmd                is a regression command (e.g,. regress, logit)
*    y_var              is the dependent variable
*    treat_var          is a dummy that equals 1 for the treatment group
*    post_var           is a dummy that equals 1 for the period in which the treatment group is treated
*    control_varlist    is a list of control variables (optional)
*
*  and options are
*    at(numlist)        a list of percentiles for CIC results. default is at(10(10)90)
*    vce(none|          don't calculate standard errors, the default
*        delta|         use numerical
*        bootstrap[, bsopts]) use bootstrap (by default, 1000 reps stratified by treat/post) other options allowed
*    untreated          counterfactual effect of the policy for the untreated group (Setion 3.2 of paper)
*    bsse               use Athey and Imben's bootstrap standard errors ( se = (p(97.5) - p(2.5)) / (2*1.96) where p(N) is Nth percentile of bootstrap iterations
*                           default is
*    *                  any other options are passed through to the regression command (if no control varlist)

* weights not allowed with vce(boostrap)

program define cic, eclass byable(recall)
	version 11.2

	// parse arguments
	syntax varlist(default=none min=3 numeric fv ts) [if] [in] [fw iw]  ///
		[, at(numlist min=1 >=0 <=100 sort) ///
		Vce(passthru) ///
		UNTreated ///
		bsse ///
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
		local wtexp_mata = `", "`wgtvar'" "'
		summ `wgtvar' if `touse' , meanonly
		local n=r(sum)
		qui count if `touse'
		local n_obs=r(N)
	}
	else {
		qui count if `touse'
		local n=r(N)
	}

	// check dummy variables
	cap assert (inlist(`treat',0,1) & inlist(`post',0,1)) if `touse'
	if _rc {
		di as error "-`treat'- and -`post'- must be dummy variables."
		error 2000
	}
	local reportdummyerror=0
	qui count if `treat'==0 & `post'==0 & `touse'
	if r(N)==0 local reportdummyerror++
	qui count if `treat'==0 & `post'==1 & `touse'
	if r(N)==0 local reportdummyerror++
	qui count if `treat'==1 & `post'==0 & `touse'
	if r(N)==0 local reportdummyerror++
	qui count if `treat'==1 & `post'==1 & `touse'
	if r(N)==0 local reportdummyerror++
	if `reportdummyerror'>1 {
		di as error "One or more of the treat*post groups is empty."
		error 2000
	}

	// parse percentiles
	if mi("`at'") local at "10(10)90" // default set (if undeclared)
	numlist "`at'"
	local at = r(numlist)

	// parse vce()
	// if vce(bootstrap), run bootstrap iterations with recursive cic call
	// this section is similar in function to the "_vce_parserun" command, except that I set default values for reps() and strata()
	if !mi("`vce'")  {
		cic_vce_parse `treat' `post' [`weight'`exp'], `vce' `bsse'
		local vce = r(vce)
		if "`vce'"=="bootstrap" {
			// bootstrap the cic command (no SEs)
return list
			bootstrap _b, reps( `=r(reps)' ) `r(bsopts)' : cic `y' `treat' `post' `varlist' if `touse' [`weight'`exp'], at(`at')
			if "bsse"=="`bsse'" {
			  if `=r(reps)' < 100 nois di as error "More than 100 bootstrap repetitions recommended for bsse option."
			  tempname ci_se
			  mata : bs_se( "e(ci_percentile)", "`ci_se'" )
			  ereturn repost V = `ci_se'
			  ereturn display
			  exit
			}
		}
	}

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
		mata: cic_caller( "`y'", "`treat'", "`post'", "`touse'", "at", `=("`vce'"=="delta")' `wtexp_mata')
	}
	else {
		// option - effect of treatment on the untreated
		// just switch `treat' variable and use -e(b)
		tempvar untreated
		gen `untreated' = !`treat'
		mata: cic_caller( "`y'", "`untreated'", "`post'", "`touse'", "at", `=("`vce'"=="delta")' `wtexp_mata')
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
	syntax [varlist] [iw fw], vce(string asis) [bsse]
	_parse comma vce 0 : vce
	if  inlist( "`vce'","bootstra","bootstr","bootst","boots","boot","boo","bo","b") local vce "bootstrap"
	if  inlist( "`vce'","delt","del","de","d")                                       local vce "delta"
	if  inlist( "`vce'","non","no","n")                                              local vce "none"
	if !inlist( "`vce'","delta","bootstrap","none") {
		di as error "Only vce(delta), vce(bootstrap [, subopts]), and vce(none) allowed."
		error 198
	}
	return local vce `vce'
	if "`vce'"=="bootstrap" {
		local stratvars : copy local varlist
		syntax [, Reps(integer 1000) strata(string asis) notable *]
		if mi("`strata'") local strata : copy local stratvars
		return local  bsopts  strata(`strata') `=cond("bsse"=="`bsse'","notable","`table'")' `options'
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


/* drop this later */ mata set matalnum on


// CIC CALLER -- THIS FUNCTION READS STATA DATA INTO MATA AND CALLS THE MAIN CIC ROUTINE
void cic_caller(string scalar y_var, string scalar treat_var, string scalar post_var, string scalar touse_var, string scalar at_local, real scalar calc_se , |string scalar wgt_var)
{
	// Inputs: Names of variables `y', `treat', `post', and `touse'
	//         Name of local macro containing quantiles of interest, ranging from 0 to 100
	//         Dummy indicating whether to calculation of SE is needed
	//         Name of variable `fweight' (optional)
	// Output: Results are returned to a Stata matrix named `cic_estimates'

	// declare variables
	real vector treat, post, y, at, result
	real vector YS, YS01
	real vector Y00, Y01, Y10, Y11
	real vector F00, F01, F10, F11

	// get data into mata
	y=treat=post=.
	st_view(y    ,. , y_var      , touse_var)
	st_view(treat,. , treat_var  , touse_var)
	st_view(post ,. , post_var   , touse_var)

	// quantiles
	at = strtoreal(tokens(st_local(at_local))) / 100
	assert((min(at)>=0) & (max(at)<=1))

	// select the rows belonging to the treat*post groups
	Y00=Y01=Y10=Y11=.
	st_select(Y00,y,(treat:==0 :& post:==0))
	st_select(Y01,y,(treat:==0 :& post:==1))
	st_select(Y10,y,(treat:==1 :& post:==0))
	st_select(Y11,y,(treat:==1 :& post:==1))

	// support points
	YS   = uniqrows( y )
	YS01 = uniqrows( Y01 )

	// CDFs (w/ and w/out weights declared)
	if (args()==7) {
		// with weight variable
		real vector wgt
		wgt=.
		st_view(wgt,.,wgt_var,touse_var)
		F00=runningsum(prob(Y00,YS,select(wgt,(treat:==0:&post:==0))))
		F01=runningsum(prob(Y01,YS,select(wgt,(treat:==0:&post:==1))))
		F10=runningsum(prob(Y10,YS,select(wgt,(treat:==1:&post:==0))))
		F11=runningsum(prob(Y11,YS,select(wgt,(treat:==1:&post:==1))))		
		F00[length(F00)]=1 // because of rounding, sum of weights might be slightly different than one
		F01[length(F01)]=1 
		F10[length(F10)]=1 
		F11[length(F11)]=1 
	}
	else {
		// with equal weights
		F00=runningsum(prob(Y00,YS))
		F01=runningsum(prob(Y01,YS))
		F10=runningsum(prob(Y10,YS))
		F11=runningsum(prob(Y11,YS))
	}

	// call the main CIC routine
	// and return result to Stata matrix named `cic_estimates'
	result = cic(F00, F01, F10, F11, YS, YS01, at)
	stata("tempname cic_estimates")
	st_matrix(st_local("cic_estimates"), result)

	// add labels
	string matrix at_lab
	at_lab = strtoname(("mean" , ("p":+strofreal(at * 100))))
	st_matrixcolstripe(st_local("cic_estimates"), ((J(1+length(at),1,"continuous") \ J(1+length(at),1,"discrete_ci") \ J(1+length(at),1,"dci_lower") \ J(1+length(at),1,"dci_upper")),J(4,1,at_lab')))
}


// CIC ROUTINE
real vector cic(real vector F00, real vector F01, real vector F10, real vector F11, real vector YS, real vector YS01, real vector at)
{
	// Inputs: Vector with CDF functions of the four treat*post groups (F00,F01,F10,F11)
	//         Vector with support points for all four groups combined (YS) and the comparison-post group (YS01)
	//         Vector with k>=1 quantiles of interest (at), ranging from 0 to 1
	// Output: Vector with 4*(1+k) results.  There are four sets; the first result in each set is the mean, followed by k results (one for each quantile).

	real vector est_con, est_dci, est_lower, est_upper

	// CIC ESTIMATOR WITH CONTINUOUS OUTCOMES, EQUATION 9
	est_con = cic_con(F00,F01,F10,F11,YS,YS01,at)

	// CIC MODEL WITH DISCRETE OUTCOMES (UNDER THE CONDITIONAL INDEPENDENCE ASSUMPTION), EQUATION 29
	est_dci = cic_dci(F00,F01,F10,F11,YS,YS01,at)

	// LOWER BOUND ESTIMATE OF DISCRETE CIC MODEL (WITHOUT CONDITIONAL INDEPENDENCE), EQUATION 25
	est_lower = cic_lower(F00,F01,F10,F11,YS,YS01,at)

	// UPPER BOUND ESTIMATE OF DISCRETE CIC MODEL (WITHOUT CONDITIONAL INDEPENDENCE), EQUATION 25
	est_upper = cic_upper(F00,F01,F10,F11,YS,YS01,at)

	return((est_con,est_dci,est_lower,est_upper))
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
	// this function calculates the sample proportions at
	// each of the support points
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
}


"check prob"
prob((1\1\2\3),(1\2\3\4\5))
prob((1\1\2\3),(1\2\3\4\5),(1\1\2\1))
prob((1\  2\3),(1\2\3\4\5),(2  \2\1))
prob((1\  2\3),(1\2\3\4\5),(1.5\1.5\.75))

"check cdfinv & cdfinv_brckt"
P = (.1\.3\.6\.7\1.0)
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
set tracedepth 3
if 0  set trace on
else  set trace off
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
local Nreps = 100
if 0  local vce vce(bootstrap, reps(`Nreps'))
else  local vce ""
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
cic  y high after ,  at(25 50 75 90) `vce'
cic ly high after ,  at(50)          `vce'
cic  y high after ,  at(25 50 75 90) `vce' untreated
cic ly high after ,  at(50)          `vce' untreated

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

// we should get an error if try to use weights
cap nois {
  cic y  high after [fw=testweight],  at(25 50 75 90) vce(bootstrap, reps(25) nodots)
}

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
