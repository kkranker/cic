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
		gen `wgtvar' `exp'
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
// add error/warning - no bootstrapping if weights.

	// check dummies
	cap assert (inlist(`treat',0,1) & inlist(`post',0,1)) if `touse'
	if _rc {
		di as error "-`treat'- and -`post'- must be dummy variables."
		error 2000
	}

	// parse percentiles
	if mi("`at'") local at "10(10)90"
	numlist "`at'"
	local at = r(numlist)

	// parse vce() 
	// if vce(bootstrap0, run bootstrap iterations with recursive cic call 
	// this section is similar in function to the "_vce_parserun" command, except that I set default values for reps() and strata()
	if !mi("`vce'")  {
		cic_vce_parse `treat' `post' [`weight'`exp'], `vce' `bsse'
		local vce = r(vce)
		if "`vce'"=="bootstrap" {			
			// bootstrap the cic command (no SEs)
			return list
			bootstrap _b, reps( `=r(reps)' ) `r(bsopts)' : cic `y' `treat' `post' `varlist' if `touse', at(`at')
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
	local run_se = `"`vce'"'=="delta"

	* adjust y for covariates using OLS regression
	if `: list sizeof varlist'!=0 {
		di as txt "Regression to adjust `y' for control variables:"
		regress `y' ib0.`treat'##ib0.`post' `varlist' if `touse' 
		tempvar yadj
		qui predict `yadj' if `touse', residuals
		qui replace `yadj' = `yadj' + _b[_cons] + _b[1.`treat']*`treat' + _b[1.`post']*`post' + _b[1.`treat'#1.`post']*`treat'*`post'
		label var `yadj' "`y' adjusted for covariates"
summ `y' `yadj' if `touse'
	}
	else local yadj `y'

	// implement mata CIC routine
	tempname cic_estimates
	if "`untreated'"=="" {
		// default - effect of treatment on the treated
		mata: cic( "`y'", "`treat'", "`post'", "`touse'", `= `wtexp_mata')
	}
	else {
		// option - effect of treatment on the untreated
		// just switch `treat' variable and use -e(b)
		tempvar untreated
		gen `untreated' = !`treat'
		mata: cic( "`y'", "`untreated'", "`post'", "`touse'", `=("`vce'"=="delta")' `wtexp_mata')
		matrix `cic_estimates' = -`cic_estimates'
	}
	
	
	// post results to ereturn
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
			"Weights not allowed with vce(bootstrap)"
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
mata set matastrict on


// CIC ROUTINE
void cic(string scalar y_var, string scalar treat_var, string scalar post_var, string scalar touse_var, real scalar calc_se , |string scalar wgt_var)
{
	// inputs: names of `y' `treat' `post' `fweight' and of local macro with percentiles

	// declare variables
	real vector treat, post, y, wgt, at
	real vector R00, R01, R10, R11
	real vector Y00, Y01, Y10, Y11
	real vector f00, f01, f10, f11
	real vector YS00, YS01, YS10, YS11, YS
	real vector F00, F01, F10, F11
	real vector est_con, est_dci, est_lower, est_upper
	real scalar i, tol

	// get data into mata
	st_view(y    ,. , y_var      , touse_var)
	st_view(treat,. , treat_var  , touse_var)
	st_view(post ,. , post_var   , touse_var)

	// quantiles
	at = strtoreal(tokens(st_local("at")))  // quantiles are stored in the stata macro variable `at' ; they are not passed as an arguement
	string matrix at_lab
	at_lab = strtoname(("mean" , ("p":+strofreal(at))))
	at = at / 100
// "at"
// at

	// identify which rows are in which group
	R00=R01=R10=R11=(.)
	for (i=1; i<=rows(treat); i++) {
		if      (treat[i]==0 & post[i]==0) R00 = ((R00==.) ? (i) : (R00 \ i))  // r is a list of rows (e.g., rows 2, 4, 5, 6, & 10)
		else if (treat[i]==0 & post[i]==1) R01 = ((R01==.) ? (i) : (R01 \ i))
		else if (treat[i]==1 & post[i]==0) R10 = ((R10==.) ? (i) : (R10 \ i))
		else                               R11 = ((R11==.) ? (i) : (R11 \ i))
	}
	// get four y vectors
	// Y00=y[R00]
	// Y01=y[R01]
	// Y10=y[R10]
	// Y11=y[R11]
	if (rowmin((rows(y[R00]), rows(y[R01]), rows(y[R10]), rows(y[R11])))==(0)) _error("One or more treat*post groups are empty")

	// support points
	YS   = uniqrows( y )
// YS00 = uniqrows( y[R00] )
	YS01 = uniqrows( y[R01] )
// YS10 = uniqrows( y[R10] )
// YS11 = uniqrows( y[R11] )

	// distributinso (w/ and w/out weights declared)
	if (args()==6) {
		st_view(wgt,. , wgt_var, touse_var)
		f00 = prob( y[R00], YS, wgt[R00] )
		f01 = prob( y[R01], YS, wgt[R01] )
		f10 = prob( y[R10], YS, wgt[R10] )
		f11 = prob( y[R11], YS, wgt[R11] )
	}
	else {
		f00 = prob( y[R00], YS )
		f01 = prob( y[R01], YS )
		f10 = prob( y[R10], YS )
		f11 = prob( y[R11], YS )
	}

// f's only called once each?  maybe just move into definition of F's

	// CDFs
	F00=runningsum(f00)
	F01=runningsum(f01)
	F10=runningsum(f10)
	F11=runningsum(f11)
// "(YS, f01, F01)"
//  (round(YS,.0001), f01, F01)

	tol = min((max(YS[2..rows(YS)]-YS[1..(rows(YS)-1)]),max(abs(f00-f01)),.00000001))/100000


	// CIC ESTIMATOR WITH CONTINUOUS OUTCOMES, EQUATION 9
	est_con = cic_con(F00,F01,F10,F11,YS,YS01,at,tol)

	// CIC MODEL WITH DISCRETE OUTCOMES (UNDER THE CONDITIONAL INDEPENDENCE ASSUMPTION), EQUATION 29
	est_dci = cic_dci(F00,F01,F10,F11,YS,YS01,at,tol)

	// LOWER BOUND ESTIMATE OF DISCRETE CIC MODEL (WITHOUT CONDITIONAL INDEPENDENCE), EQUATION 25
	est_lower = cic_lower(F00,F01,F10,F11,YS,YS01,at,tol)

	// UPPER BOUND ESTIMATE OF DISCRETE CIC MODEL (WITHOUT CONDITIONAL INDEPENDENCE), EQUATION 25
	est_upper = cic_upper(F00,F01,F10,F11,YS,YS01,at,tol)

// ( "est_con" \ "est_dci" \ "est_lower" \ "est_upper" )
// at_lab
// round((est_con \ est_dci \  est_lower \  est_upper ),.01)

 st_eclear()
 st_matrix(st_local("cic_estimates"), (est_con , est_dci ,  est_lower ,  est_upper ))
 st_matrixcolstripe(st_local("cic_estimates"), ((J(cols(est_con),1,"continuous") \ J(cols(est_con),1,"discrete_ci") \ J(cols(est_con),1,"dci_lower") \ J(cols(est_con),1,"dci_upper")),J(4,1,at_lab')))

} // end of m_cic


// CIC ESTIMATOR WITH CONTINUOUS OUTCOMES, EQUATION 9
real vector cic_con(real vector F00, real vector F01, real vector F10, real vector F11, real vector YS, real vector YS01, real vector at, real vector tol)
{
	// this function calculates the continuous outcomes CIC estimator
	// first estimate the cdf of Y^N_11 using equation (9) in the paper and
	// then use that to calculate the average effect of the treatment
	real vector FCO, est_con
	real scalar i, F01y, F00invF01y, F10F00invF01y

	// for each y in the support of Y01, fill in FCO(y)=F_10(F^-1_00(F_01(y)))
	FCO=J(rows(YS01),1,0)
	for(i=1; i<=rows(YS01); ++i) {
		F01y=cdf(YS01[i],F01,YS,tol)
		F00invF01y=cdfinv(F01y,F00,YS,tol)
		F10F00invF01y=cdf(F00invF01y,F10,YS,tol)
		FCO[i] = F10F00invF01y
	}
	FCO[rows(FCO)]=1 // =1 at end

	// mean CIC estimate
	est_con=( (F11-(0 \ F11[1..(rows(YS)-1)]))'*YS - (FCO-(0 \ FCO[1..(rows(YS01)-1)]))'*YS01 )

	// quantile CIC estimate
	for(i=1; i<=cols(at); ++i) {
		est_con = (est_con , ( cdfinv(at[i], F11, YS, tol) - cdfinv(at[i], FCO, YS01, tol) ) )
	}

	// matrix has mean estimate in first column, plus one column for each element of "at"
	return(est_con)
}


// CIC MODEL WITH DISCRETE OUTCOMES (UNDER THE CONDITIONAL INDEPENDENCE ASSUMPTION), EQUATION 29
real vector cic_dci(real vector F00, real vector F01, real vector F10, real vector F11, real vector YS, real vector YS01, real vector at, real vector tol)
{
	// this function calculates the discreate outcomes CIC estimator
	// first estimate the cdf of Y^N_11 using equation (29) in the paper and
	// then use that to calculate the average effect of the treatment
	real vector FDCI, FUB, FLB, est_dci
	real scalar i,F01y,F00invF01y,F10F00invF01y,F00invbF01y,F10F00invbF01y,F00F00invF01y,F00F00invbF01y,FDCI_weight

	// for each y in the support of Y01, fill in FCO(y)=F_10(F^-1_00(F_01(y)))
	FDCI=FUB=FLB=J(rows(YS01),1,0)
	for(i=1; i<=rows(YS01); ++i) {
		F01y=cdf(YS01[i],F01,YS,tol)
		F00invF01y=cdfinv(F01y,F00,YS,tol)
		F10F00invF01y=cdf(F00invF01y,F10,YS,tol)
		F00invbF01y=cdfinv_brckt(F01y,F00,YS,tol)
		F10F00invbF01y=cdf(F00invbF01y,F10,YS,tol)
		F00F00invF01y =cdf(F00invF01y ,F00,YS,tol)
		F00F00invbF01y=cdf(F00invbF01y,F00,YS,tol)
		FLB[i]=F10F00invbF01y
		FUB[i]=F10F00invF01y
		if ((F00F00invF01y-F00F00invbF01y)>tol)    FDCI_weight=(F01y-F00F00invbF01y)/(F00F00invF01y-F00F00invbF01y)
		else                                       FDCI_weight=0
		FDCI[i]=FLB[i]+(FUB[i]-FLB[i])*FDCI_weight
	}
	FDCI[rows(FDCI)]=1 // =1 at end

	// conditional independence estimate
	est_dci=( (F11-(0 \ F11[1..(rows(YS)-1)]))'*YS - (FDCI-(0 \ FDCI[1..(rows(YS01)-1)]))'*YS01 )


	// quantile CIC estimate
	for(i=1; i<=cols(at); ++i) {
		est_dci = (est_dci , ( cdfinv(at[i], F11, YS, tol) - cdfinv(at[i], FDCI, YS01, tol) ) )
	}

	// matrix has mean estimate in first column, plus one column for each element of "at"
	return(est_dci)
}


// LOWER BOUND ESTIMATE OF DISCRETE CIC MODEL (WITHOUT CONDITIONAL INDEPENDENCE), EQUATION 25
real vector cic_lower(real vector F00, real vector F01, real vector F10, real vector F11, real vector YS, real vector YS01, real vector at, real vector tol)
{
	// this function calculates the discreate outcomes CIC estimator
	// first estimate the cdf of Y^N_11 using equation (29) in the paper and
	// then use that to calculate the average effect of the treatment
	real vector FLB, est_lower
	real scalar i,F01y, F00invbF01y,F10F00invbF01y

	// for each y in the support of Y01, fill in FCO(y)=F_10(F^-1_00(F_01(y)))
	FLB=J(rows(YS01),1,0)
	for(i=1; i<=rows(YS01); ++i) {
		F01y=cdf(YS01[i],F01,YS,tol)
		F00invbF01y=cdfinv_brckt(F01y,F00,YS,tol)
		F10F00invbF01y=cdf(F00invbF01y,F10,YS,tol)
		FLB[i]=F10F00invbF01y;
	}
	FLB[rows(FLB)]=1 // =1 at end

	// conditional independence estimate
	est_lower=( (F11-(0 \ F11[1..(rows(YS)-1)]))'*YS - (FLB-(0 \ FLB[1..(rows(YS01)-1)]))'*YS01 )

	// quantile CIC estimate
	for(i=1; i<=cols(at); ++i) {
		est_lower = (est_lower , ( cdfinv(at[i], F11, YS, tol) - cdfinv(at[i], FLB, YS01, tol) ) )
	}

	// matrix has mean estimate in first column, plus one column for each element of "at"
	return(est_lower)
}


// UPPER BOUND ESTIMATE OF DISCRETE CIC MODEL (WITHOUT CONDITIONAL INDEPENDENCE), EQUATION 25
real vector cic_upper(real vector F00, real vector F01, real vector F10, real vector F11, real vector YS, real vector YS01, real vector at, real vector tol)
{
	// this function calculates the discreate outcomes CIC estimator
	// first estimate the cdf of Y^N_11 using equation (29) in the paper and
	// then use that to calculate the average effect of the treatment
	real vector FUB, est_upper
	real scalar i,F01y,F00invF01y,F10F00invF01y

	// for each y in the support of Y01, fill in FCO(y)=F_10(F^-1_00(F_01(y)))
	FUB=J(rows(YS01),1,0)
	for(i=1; i<=rows(YS01); ++i) {
		F01y=cdf(YS01[i],F01,YS,tol)
		F00invF01y=cdfinv(F01y,F00,YS,tol)
		F10F00invF01y=cdf(F00invF01y,F10,YS,tol)
		FUB[i]=F10F00invF01y;
	}
	FUB[rows(FUB)]=1 // =1 at end

	// conditional independence estimate
	est_upper=( (F11-(0 \ F11[1..(rows(YS)-1)]))'*YS - (FUB-(0 \ FUB[1..(rows(YS01)-1)]))'*YS01 )

	// quantile CIC estimate
	for(i=1; i<=cols(at); ++i) {
		est_upper = (est_upper , ( cdfinv(at[i], F11, YS, tol) - cdfinv(at[i], FUB, YS01, tol) ) )
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
	// wgt is an (optional) et of weights for the vector Y
	real scalar n, tol
	n   = rows(YS)
	tol = colmin(abs(YS[2..n]-YS[1..(n-1)]))/2  // as tolerance, use half distance between nearest two support points
	if (args()==2) return( rowsum( abs( (YS :- J(n,1,Y')) ) :<= tol                 ) :/ rows(Y)           )
	else           return( rowsum((abs( (YS :- J(n,1,Y')) ) :<= tol ) :* J(n,1,wgt')) :/ J(n,1,colsum(wgt)))
} // end of prob


// CUMULATIVE DISTRIBUTION FUNCTION
real scalar cdf(real scalar y, real vector P, real vector YS, real scalar tol)
{
	// given a cumulative distrubtion function (P) over the support points (YS),
	// returns the empirical cumulative distribution function at a scalar (y)
	if      (y< min(YS)) return(0)
	else if (y>=max(YS)) return(1)
	else                 return(P[colsum((YS:<=(y+tol)))])
}


// INVERSE OF CUMULATIVE DISTRIBUTION FUNCTION, EQUATION 8
real scalar cdfinv(real scalar p, real vector P, real vector YS, real scalar tol)
{
	// given a cumulative distrubtion functin (P) over the support points (YS),
	// returns the inverse of the empirical cumulative distribution function at probability p (0<p<1)
	real scalar r
	r=1
	while (P[r]<(p-tol) & r<rows(YS)) {
		r++
	}
	return(YS[r])
}


// INVERSE OF CUMULATIVE DISTRIBUTION FUNCTION, ALTERNATIVE FOR DISCREATE OUTCOMES, EQUATION 24
real scalar cdfinv_brckt(real scalar p, real vector P, real vector YS, real scalar tol)
{
	// given a cumulative distrubtion functin (P) over the support points (YS),
	// returns the inverse of the empirical cumulative distribution function at probability p (0<p<1)
	// but if equals -oo, it returns min(YS)-100*(1+max(YS)-YS(min)) = 101*min(YS)-100*max(YS)-100
	real scalar r, n
	n=rows(YS)
	r=1
	if (p>=P[1]-tol) {
		while (P[r]<=(p+tol) & r<n) {
			r++
		}
		return(YS[r-1])
	}
	else {
		return(101*YS[1]-100*YS[n]-100)
	}
}


// USE 95 PERCENTILES OF BOOTSTRAP ITERATIONS TO BACKOUT STANDARD ERRORS
void bs_se( string scalar in_ci, string scalar out_V )
{
	real    vector bs_se
	bs_se = (1 /(2 * 1.96)) * (st_matrix(in_ci)[2,.] - st_matrix(in_ci)[1,.])
	st_matrix(out_V,diag(bs_se:*bs_se))
	st_matrixcolstripe(out_V,st_matrixcolstripe(in_ci))
	st_matrixrowstripe(out_V,st_matrixcolstripe(in_ci))
}


"check prob"
prob((1\1\2\3),(1\2\3\4\5))
prob((1\1\2\3),(1\2\3\4\5),(1\1\2\1))
prob((1\  2\3),(1\2\3\4\5),(2  \2\1))

"check cdfinv & cdfinv_brckt"
tol = .0000000000000000000001
cdfinv(.05    , (.1\.3\.6\.7\1.0) , (1\2\3\4\5), tol)
cdfinv(.2     , (.1\.3\.6\.7\1.0) , (1\2\3\4\5), tol)
cdfinv(.3     , (.1\.3\.6\.7\1.0) , (1\2\3\4\5), tol)
cdfinv(.300001, (.1\.3\.6\.7\1.0) , (1\2\3\4\5), tol)

cdfinv_brckt(.05    , (.1\.3\.6\.7\1.0) , (1\2\3\4\5), tol)
cdfinv_brckt(.2     , (.1\.3\.6\.7\1.0) , (1\2\3\4\5), tol)
cdfinv_brckt(.299999, (.1\.3\.6\.7\1.0) , (1\2\3\4\5), tol)
cdfinv_brckt(.3     , (.1\.3\.6\.7\1.0) , (1\2\3\4\5), tol)
cdfinv_brckt(.300001, (.1\.3\.6\.7\1.0) , (1\2\3\4\5), tol)


end
/* * * * *  END OF MATA BLOCK * * * * */

cic_vce_parse, vce(boot, reps(1000) saving(myfile.dta, replace))
return list
cic_vce_parse, vce(boot, reps(25) strata(x y) cluster(treat))
return list
cic_vce_parse, vce(boot, reps(25) strata() cluster())
return list


set seed 1


* The following code can be used to test the program using "fake" data
// exit

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* with athey and imbens data  suppliment
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

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
tabstat y ly , by(high_after) s(count mean sd min p25 p50 p75 p90 max) columns(s)  labelwidth(30) nototal format(%9.2f)

count

reg y high##after
reg ly high##after

set tracedepth 2
if 01 set trace on 
else  set trace off
if 01 local vce vce(bootstrap, reps(50))
else  local vce ""

cic  y high after ,  at(25 50 75 90) `vce'
cic ly high after ,  at(50)          `vce'
cic  y high after ,  at(25 50 75 90) `vce' untreated
cic ly high after ,  at(50)          `vce' untreated


exit


cic y  high after ,  at(25 50 75 90)



cic y  high after ,  at(25 50 75 90) vce(bootstrap, reps(25) nodots)


bootstrap , reps(25) strata(high after) ///
	: cic y  high after ,  at(25 50 75 90)



// bootstrap the sample conditional on Ngt for g; t = 0; 1

gen testweight =1 + (uniform()<.20)
cic y  high after [fw=testweight],  at(25 50 75 90)


ereturn list


* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
sysuse nlsw88, clear
set seed 1
gen TREAT1 = uniform() < .5
replace wage = wage + TREAT1
gen POST1 = uniform() < .5
replace wage = wage - POST1

keep if uniform()<.012

expand = 2*(uniform()<.1)

// bootstrap the sample conditional on Ngt for g; t = 0; 1
cic wage TREAT1 POST1 i. occupation ,  at(10(10)90 99.5)

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
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
