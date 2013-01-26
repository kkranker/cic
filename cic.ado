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
*    y_var                    is the dependent variable
*    treat_var                is a dummy that equals 1 for the treatment group
*    post_var                 is a dummy that equals 1 for the period in which the treatment group is treated
*    control_varlist          is a list of control variables (optional)
*
*  and options are as follows:
*
*  MAIN
*    at(numlist)              a list of percentiles for CIC results. default is at(10(10)90)
*    vce(none|                don't calculate standard errors, the default
*        delta|               use numerical
*        bootstrap[, bsopts]| use bootstrap (by default, 1000 reps stratified by treat/post) other options allowed
*        bspctile[, bsopts])  use Athey and Imben's of obtaining bootstrap standard errors ( se = (p(97.5) - p(2.5)) / (2*1.96) where p(N) is Nth percentile of bootstrap iterations
*                                 instead of using Stata's default method. The Stata manual ([R] bootstrap) addresses this issue.
*                                    "A better estimate is needed if you want to use the 2.5th and 97.5th percentiles of the distribution
*                                    to produce a 95% confidence interval. To extract many features simultaneously about the distribution,
*                                    an even better estimate is needed. Generally, replications on the order of 1,000 produce very good
*                                    estimates, but only 50–200 replications are needed for estimates of standard errors...."
*    untreated                counterfactual effect of the policy for the untreated group (Setion 3.2 of paper)
*
*  REPORTING
*      level(passthru)              set confidence level; default is level(95)
*      notable                      suppress table of results
*      noheader                     suppress table header
*      nolegend                     suppress table legend
*      display_options              control spacing and display of omitted variables and base and empty cells
*								       display_options:  noomitted, vsquish, noemptycells, baselevels, allbaselevels;
*                                                        see help estimation options##display_options
*
*  BSOPTS SUBOPTIONS
*     reps(#)                      perform # bootstrap replications; default is reps(200)
*     saving(filename[,replace])   save bootstrap results to filename (optionally, replace specifies that filename be overwritten, if it exists.)
*     accel(vector)                acceleration values for each statistic
*     mse                          use MSE formula for variance estimation
*     nodots                       suppress the replication dots
*  See [R] bootstrap postestimation for features available after estimation.



* Weights may be iweights or fweights.  Weights are not allowed with vce(boostrap).

program define cic, eclass byable(recall)
	version 11.2
	ereturn clear

	// parse arguments
	syntax varlist(default=none min=3 numeric fv ts) [if] [in] [fw iw]  ///
		[, at(numlist min=1 >=0 <=100 sort) ///
		Vce(passthru) ///
		UNTreated ///
		level(passthru) notable NOHeader NOLegend * ] // Reporting options
	marksample touse
	 _get_diopts diopts, `options'
	local diopts `diopts' `table' `header' `legend'

	// first three variables need to be y, treat, and post
	gettoken y     varlist : varlist
	gettoken treat varlist : varlist
	gettoken post  varlist : varlist
// byable(recall) working right?

// _error( "also, switch bootstrap --> fullbootstrap and bspctile as a followup to bstat ")
// error -- add check non-integer weights allowed with cic bootstrap
// add documentation that *  fweights, but not iweights, work with vce( ??????????????? )
// in documentation, talk about why bootstrap: prefix might be needed (longitudinal data, sample by id instead of pre/post groups)
// vce(????) is equivalent to bootstrap:

// vce(bootstrap, [bsopts]) is equivalent to
//    . bootstrap _b, strata(treat post) [bsopts]: cic y treat post ... , vce(none)
// but slower because vce(bootstrap) is implimented in META and runs with less overhead.
// However, the bootstrap prefix is more flexible due the availability of size(), strata(), cluster(), idcluster() and other variables.

	// prep to handle weights
	if !missing("`weight'") {
		tempvar wgtvar
		gen `wgtvar'`exp' if `touse'
		local wtexp_caller = `", "`wgtvar'" "'
		summ `wgtvar' if `touse' , meanonly
		local n=round(r(sum))
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
	if mi(`"`vce'"') local vce vce(none)
	cic_vce_parse, `vce'
	local vce     = r(vce)
	local bsreps  = r(bsreps)
	local nodots  = r(nodots)
	local sedelta = r(sedelta)
	if !missing(r(bsopts)) local bsopts  = r(bsopts)
	if !missing(r(saving)) local bsopts  = r(saving)
	
	// adjust y for covariates (OLS regression)
	if `: list sizeof varlist'!=0 {
		di as txt "Regression to adjust `y' for control variables:"
		regress `y' ib0.`treat'##ib0.`post' `varlist' if `touse' [`weight'`exp'], `diopts' `level'
		tempvar yresid yadj
		qui predict `yresid' if `touse', residuals
		qui gen `yadj' = `yresid' + _b[_cons] + _b[1.`treat']*`treat' + _b[1.`post']*`post' + _b[1.`treat'#1.`post']*`treat'*`post' if `touse'
		label var `yresid' "`y' residuals after adjusting for covariates"
		label var `yadj' "`y' adjusted for covariates"
	}
	else local yadj `y'

	// implement mata CIC routine
	tempname eresults
	if "`untreated'"=="" {
		// default - effect of treatment on the treated
		mata: cicresult = cic_caller( "`yadj'", "`treat'",     "`post'", "`touse'", "at", `bsreps', `nodots', `sedelta' `wtexp_caller')
		mat `eresults' = e(results)
	}
	else {
		// option - effect of treatment on the untreated (just switch `treat' variable and use -e(results))
		tempvar untreated
		gen `untreated' = (`treat'==0)
		mata: cicresult = cic_caller( "`yadj'", "`untreated'", "`post'", "`touse'", "at", `bsreps', `nodots', `sedelta' `wtexp_caller')
		mat `eresults' = e(results)
		matrix `eresults' = -`eresults'
	}

	// post results to ereturn
	if `bsreps' {
		// option - save bs reps into a specified file (`saving' contains ", replace" if it was provided)
		if !mi(`"`saving'"') copy `bstempfile' `saving'
		
		// post boostrapped estimates with standard errors using bstat
		bstat  using `bstempfile', stat( `eresults' ) `bsopts' `diopts' `level' title(Changes in Changes (CIC) Estimation)

		if "bspctile"=="`vce'" {
			if `=r(reps)' < 1000 nois di as error "Warning: More bootstrap repetitions might be needed with vce(bspctile)."
			tempname ci_se
			mata : bs_se( "e(ci_percentile)", "`ci_se'" )
			ereturn repost V = `ci_se'
			ereturn display
		}

	}
	else {
		// otherwise just use ereturn to get pretty estimates table
	    di as txt _n "Changes in Changes (CIC) Estimation"
		ereturn post `eresults' [`weight'`exp'], depname(`y') obs(`n') esample(`touse') dof(`=`n'-colsof(`eresults')') `level'
		ereturn display
	}
 	ereturn local cmd = "cic"
 	ereturn local cmdline cic `0'
	if  !missing("`untreated'")  di as txt "(Effect of Treatment on the Treated Group)"
	else                         di as txt "(Effect of Treatment on the Untreated Group)"
	if !missing("`untreated'")   ereturn local effecton "Effect of Treatment on the Treated Group"
	else                         ereturn local effecton "Effect of Treatment on the Untreated Group"

// graph results
/*
getmata

matrix b = e(b)
matrix b = b'
keep in 1/3
svmat b, names("b")
*/
ereturn list
char list

end // end of cic program definition


// subroutine to parse the vce() option
// this section is similar in function to the "_vce_parse" command, except that I set default values for reps() and strata()
program define cic_vce_parse, rclass
	version 11.2
	syntax , vce(string asis)
	_parse comma vce 0 : vce
	if  inlist( "`vce'","bootstra","bootstr","bootst","boots","boot","boo","bo")     local vce "bootstrap"
	if  inlist( "`vce'","bspctil","bspcti","bspct","bspc","bsp","bsp","bs")          local vce "bspctile"
	if  inlist( "`vce'","delt","del","de","d")                                       local vce "delta"
	if  inlist( "`vce'","non","no","n")                                              local vce "none"

	if !inlist("`vce'","bootstrap","bspctile") & !mi("`0'") {
		di as error "suboptions are not allowed with vce(`vce')"
		error 198
	}

	return local vce `vce'
	if inlist("`vce'","bootstrap","bspctile") {
		syntax [, Reps(integer 200) saving(string asis) NODots *]
		return scalar bsreps  = `reps'
		return scalar nodots  = ( "nodots"=="`dots'" )
		return scalar sedelta = 0
		return local  saving  : copy local saving
		return local  bsopts  : copy local options
	}
	else if ("`vce'"=="delta") {
		return scalar bsreps  = 0
		return scalar nodots  = 0
		return scalar sedelta = 1

	}
	else if ("`vce'"=="none") {
		return scalar bsreps  = 0
		return scalar nodots  = 0
		return scalar sedelta = 0
	}
	else {
		di as error "Only vce(delta), vce(bootstrap [, subopts]), vce(bspctile [, subopts]), and vce(none) allowed."
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
	// main results
	real rowvector con, dci, dcilowbnd, dciuppbnd

	// store a few inputs
	real rowvector at
}


// CIC CALLER -- THIS FUNCTION READS STATA DATA INTO MATA AND CALLS THE MAIN CIC ROUTINE
struct cic_result scalar cic_caller(string scalar y_var, string scalar treat_var, string scalar post_var, string scalar touse_var, string scalar at_local, real scalar bsreps, real scalar nodots, real scalar sedelta, |string scalar wgt_var)
{
	// Inputs:
	//   1-4. Names of variables `y', `treat', `post', and `touse'
	//   5.   Name of local macro containing quantiles of interest, ranging from 0 to 100 (`at')
	//   6.   Number of bootstrap reps (=0 to skip)
	//   7.   Do not show bootstrapping dots (=1 to suppress dots)
	//   8.   Standard error for conditional independence based on numerical differentiation (=0 to skip)
	//   9.   (Optional) Name of variable with fweight or iweight
	// Output: Structure with four vectors. Each vector has (1+k) elements. The first element is the mean, followed by k results (one for each quantile in -at-).

	// Get data into mata
	real colvector y, treat, post
	st_view(y=.    ,.,y_var    ,touse_var)  // note that rows with missing data are already dropped by -marksample- in .ado file
	st_view(treat=.,.,treat_var,touse_var)
	st_view(post=. ,.,post_var ,touse_var)

	// Quantiles
	real rowvector at
	at = strtoreal(tokens(st_local(at_local))) / 100
	if (min(at)<0 | max(at)>1) _error( "at() must be between zero and one (inclusive)." )

	// Select the rows belonging to the treat*post groups
	real colvector Y00, Y01, Y10, Y11
	st_select(Y00=.,y,(treat:==0 :& post:==0))
	st_select(Y01=.,y,(treat:==0 :& post:==1))
	st_select(Y10=.,y,(treat:==1 :& post:==0))
	st_select(Y11=.,y,(treat:==1 :& post:==1))
	if (min((length(Y00),length(Y01),length(Y10),length(Y11)))==0) _error( "One or more of the four treat*post groups is empty." )

	// Select the rows of wgt belonging to the treat*post groups
	if (args()==9) {
		real colvector wgt, W00, W01, W10, W11
		st_view(wgt=.,.,wgt_var,touse_var)
		st_select(W00=.,wgt,(treat:==0 :& post:==0))
		st_select(W01=.,wgt,(treat:==0 :& post:==1))
		st_select(W10=.,wgt,(treat:==1 :& post:==0))
		st_select(W11=.,wgt,(treat:==1 :& post:==1))
	}

	// Results will be returned into a structure w/ 4 vectors for con, dci, lower, upper
	struct cic_result scalar result

	// Call the main CIC routine
	if (args()==8) result=cic(Y00,Y01,Y10,Y11,at)                 // without weights
	else           result=cic(Y00,Y01,Y10,Y11,at,W00,W01,W10,W11) // with weights

	// Return results into a Stata matrix too
	string matrix colfulllabels
	st_matrix("e(at)",(result.at))
	st_matrix("e(results)", (result.con,result.dci,result.dcilowbnd,result.dciuppbnd))

	// Label the columns of the matrix
	colfulllabels=((J(1+length(at),1,"continuous") \ J(1+length(at),1,"discrete_ci") \ J(1+length(at),1,"dci_lower_bnd") \ J(1+length(at),1,"dci_upper_bnd")),J(4,1,strtoname(("mean" , ("p":+strofreal(at*100))))'))
	st_matrixcolstripe("e(results)", colfulllabels)
	st_local("cic_coleq"   ,invtokens(colfulllabels[.,1]'))
	st_local("cic_colnames",invtokens(colfulllabels[.,2]'))
	st_numscalar( "e(N_strata)", 4)
	if (args()==8) st_numscalar( "e(N)"       , rows(y))
	else           st_numscalar( "e(N)"       , round(sum(wgt)))
	st_numscalar( "e(N00)"     , rows(Y00))
	st_numscalar( "e(N01)"     , rows(Y01))
	st_numscalar( "e(N10)"     , rows(Y10))
	st_numscalar( "e(N11)"     , rows(Y11))

	// Bootstrapping
	if (bsreps & sedelta) _error( "bsreps and sedelta not allowed at the same time.")
	else if (bsreps) {
		real scalar b
		struct cic_result scalar bs_loop
		real matrix bsdata
		bsdata=J(bsreps,4*(1+length(at)),.)

		// header for dots
		if (!nodots) {
			printf( "{txt}\nBootstrap replications ({res}%g{txt})\n", bsreps)
			display( "{txt}{hline 4}{c +}{hline 3} 1 " +
				"{hline 3}{c +}{hline 3} 2 " + "{hline 3}{c +}{hline 3} 3 " +
				"{hline 3}{c +}{hline 3} 4 " + "{hline 3}{c +}{hline 3} 5 ")
		}

		if ((args()==8) & (round(wgt)!=wgt)) _error( "CIC bootstrapping does not work with iweights." )

		for(b=1; b<=bsreps; ++b) {
			if (args()==8) bs_loop=cic(drawsmpl(Y00),drawsmpl(Y01),drawsmpl(Y10),drawsmpl(Y11),at)                 // without weights
			else           bs_loop=cic(drawsmpl(Y00,W00),drawsmpl(Y01,W01),drawsmpl(Y10,W10),drawsmpl(Y11,W11),at) // with frequency weights

			// save into return structure
			bsdata[b,.]  =(bs_loop.con,bs_loop.dci,bs_loop.dcilowbnd,bs_loop.dciuppbnd)

			// show dots
			if (!nodots) {
				if (missing((bs_loop.con,bs_loop.dci,bs_loop.dcilowbnd,bs_loop.dciuppbnd))) printf( "{err}x{txt}")
				else printf( ".")
				if (!mod(b,50)) printf( " %5.0f\n",b)
				displayflush()
			}
		} // end loop through bs iterations
		if (!nodots & mod(b-1,50)) display("")

		// save bootstrap iterations in a temporary .dta file (named `bstempfile')
		stata( "preserve" )
		  string rowvector bstempfile, bstempvars
		  // clear data (after preserve) and fill with bsdata matrix
		  st_dropvar(.)
		  st_addobs(rows(bsdata))
		  bstempvars=strtoname("est":+strofreal(1::cols(bsdata)))'
		  st_store(.,st_addvar("double",bstempvars), bsdata)
		  // setup file for bstat command
		  st_global( "_dta[bs_version]" , "3")
		  if (args()==8) st_global( "_dta[N]", strofreal(rows(y)))
		  else           st_global( "_dta[N]", strofreal(round(sum(wgt))))
		  st_global( "_dta[N_strata]"   , "4")
		  st_global( "_dta[strata]"     , (treat_var + " " + post_var))
		  st_global( "_dta[command]"    , "cic")
		  for(b=1; b<=cols(bsdata); ++b) {
			 st_global( (bstempvars[1,b]+"[colname]"), colfulllabels[b,2])
			 st_global( (bstempvars[1,b]+"[coleq]")  , colfulllabels[b,1])
		  }
		  // save as `bstempfile'
		  bstempfile=st_tempfilename()
		  st_local( "bstempfile",bstempfile)
		  stata(( "qui save " + bstempfile ))
		stata( "restore" )
	} // done bootstrapping
	else if (sedelta) {
_error( "Code for sedelta=0 not written." )
	}
	

	// DONE.  Pass results back to caller
	return(result)
}


// >>>>>>>>>>  check that column names and such are the same as the ouput in the example in the NOTE (below)

// CIC ROUTINE
struct cic_result scalar cic(real colvector Y00, real colvector Y01, real colvector Y10, real colvector Y11, real rowvector at, | real colvector W00, real colvector W01, real colvector W10, real colvector W11 )
{
	// Inputs:
	//   (1)-(4) Four column vectors with data.
	//            - Y00 is control group in pre-period
	//            - Y01 is control group in post period
	//            - Y10 is treatment group in post period
	//            - Y11 is treatment group in post period
	//   (5)     Vector with k>=1 quantiles of interest, ranging from 0 to 1
	//   (6)-(9) (Optional) Column with fweights or iweights for Y00, Y01, Y10, and Y11 (respectively)
	//
	// Output: One structure (cic_result) with four row vectors.
	//   Each vector has (1+k) elements. The first element is the mean, followed by k results (one for each quantile in -at-).
	//   (1) result.con       = CIC ESTIMATOR WITH CONTINUOUS OUTCOMES, EQUATION 9
	//   (2) result.dci       = CIC MODEL WITH DISCRETE OUTCOMES (UNDER THE CONDITIONAL INDEPENDENCE ASSUMPTION), EQUATION 29
	//   (3) result.dcilowbnd = LOWER BOUND ESTIMATE OF DISCRETE CIC MODEL (WITHOUT CONDITIONAL INDEPENDENCE), EQUATION 25
	//   (4) result.dciuppbnd = UPPER BOUND ESTIMATE OF DISCRETE CIC MODEL (WITHOUT CONDITIONAL INDEPENDENCE), EQUATION 25

	// The code in cic() is somewhat convoluted because I am
	// calculating all four vectors simultaneously.
	// See the NOTE (below, at the bottom of the file)
	// for alternative routines that are much more clear
	// and accessible. The alternatives produce
	// identical results. However redundant calculations lead to
	// slower run-times.

	// Need all or none of args (6)-(9)
	if (args()>5 & args()!=9) _error(( "Expected 5 or 9 arguements, but received " + strofreal(args())))

	// Vector with support points for all four groups combined (YS) and the comparison-post group (YS01)
	real colvector YS, YS01
	YS01 = uniqrows(Y01)
	YS   = uniqrows(Y00\YS01\Y10\Y11)

	// Vector with CDF functions of the four treat*post groups (F00,F01,F10,F11)
	// CDFs (w/ and w/out weights declared)
	real colvector F00, F01, F10, F11
	if (args()==5) {
		// CDFs without weights
		F00=runningsum(prob(Y00,YS))
		F01=runningsum(prob(Y01,YS))
		F10=runningsum(prob(Y10,YS))
		F11=runningsum(prob(Y11,YS))
	}
	else {
		// CDFs with weights
		F00=runningsum(prob(Y00,YS,W00))
		F01=runningsum(prob(Y01,YS,W01))
		F10=runningsum(prob(Y10,YS,W10))
		F11=runningsum(prob(Y11,YS,W11))
		// because of rounding, sum of weights might be slightly different than one
		F00[length(F00)]=1
		F01[length(F01)]=1
		F10[length(F10)]=1
		F11[length(F11)]=1
	}

	// First estimate the cdf of Y^N_11 using equation (9) in the paper and
	// then use that to calculate the average effect of the treatment.
	// For each y in the support of Y01, fill in FCO(y)=F_10(F^-1_00(F_01(y))).
	real vector FCO,FLB,FUB,FDCI,FDCI_weight
	real scalar i,F01y,F00invF01y,F00invbF01y,F00F00invF01y,F00F00invbF01y
	FCO=FDCI=FLB=FUB=J(length(YS01),1,0)
	for(i=1; i<=length(YS01); ++i) {
		F01y=cdf(YS01[i],F01,YS)
		F00invF01y=cdfinv(F01y,F00,YS)
		F00invbF01y=cdfinv_brckt(F01y,F00,YS)
		F00F00invF01y =cdf(F00invF01y,F00,YS)
		F00F00invbF01y=cdf(F00invbF01y,F00,YS)
		FCO[i]=FUB[i]=cdf(F00invF01y,F10,YS)
		FLB[i]=cdf(F00invbF01y,F10,YS)
		if ((F00F00invF01y-F00F00invbF01y)>epsilon(1)) FDCI_weight=(F01y-F00F00invbF01y)/(F00F00invF01y-F00F00invbF01y)
		else                                           FDCI_weight=0
		FDCI[i]=FLB[i]+(FUB[i]-FLB[i])*FDCI_weight
	}
	FCO[length(FCO)]=FDCI[length(FDCI)]=FLB[length(FLB)]=FUB[length(FUB)]=1   // =1 in last row

	// Results will be returned into a structure w/ 4 vectors for con, dci, lower, upper
	struct cic_result scalar result

	// CIC ESTIMATOR WITH CONTINUOUS OUTCOMES, EQUATION 9
	// calculate the continuous outcomes CIC estimator
	// matrix has mean estimate in first column, plus one column for each element of "at"
	// mean CIC estimate
	result.con=( (F11-(0 \ F11[1..(length(YS)-1)]))'*YS - (FCO-(0 \ FCO[1..(length(YS01)-1)]))'*YS01 )
	// quantile CIC estimates
	for(i=1; i<=length(at); ++i) {
		result.con = (result.con , ( cdfinv(at[i], F11, YS) - cdfinv(at[i], FCO, YS01) ) )
	}


	// CIC MODEL WITH DISCRETE OUTCOMES (UNDER THE CONDITIONAL INDEPENDENCE ASSUMPTION), EQUATION 29
	// calculate the discreate outcomes CIC estimator
	// matrix has mean estimate in first column, plus one column for each element of "at"
	// conditional independence estimate
	result.dci=( (F11-(0 \ F11[1..(length(YS)-1)]))'*YS - (FDCI-(0 \ FDCI[1..(length(YS01)-1)]))'*YS01 )
	// quantile CIC estimates
	for(i=1; i<=length(at); ++i) {
		result.dci = (result.dci , ( cdfinv(at[i], F11, YS) - cdfinv(at[i], FDCI, YS01) ) )
	}


	// LOWER BOUND ESTIMATE OF DISCRETE CIC MODEL (WITHOUT CONDITIONAL INDEPENDENCE), EQUATION 25
	// calculate the discreate outcomes CIC estimator
	// conditional independence estimate
	result.dcilowbnd =( (F11-(0 \ F11[1..(length(YS)-1)]))'*YS - (FLB-(0 \ FLB[1..(length(YS01)-1)]))'*YS01 )
	// quantile CIC estimates
	for(i=1; i<=length(at); ++i) {
		result.dcilowbnd  = (result.dcilowbnd  , ( cdfinv(at[i], F11, YS) - cdfinv(at[i], FLB, YS01) ) )
	}


	// UPPER BOUND ESTIMATE OF DISCRETE CIC MODEL (WITHOUT CONDITIONAL INDEPENDENCE), EQUATION 25
	// calculate the discreate outcomes CIC estimator
	// matrix has mean estimate in first column, plus one column for each element of "at"
	// conditional independence estimate
	result.dciuppbnd=( (F11-(0 \ F11[1..(length(YS)-1)]))'*YS - (FUB-(0 \ FUB[1..(length(YS01)-1)]))'*YS01 )
	// quantile CIC estimates
	for(i=1; i<=length(at); ++i) {
		result.dciuppbnd = (result.dciuppbnd , ( cdfinv(at[i], F11, YS) - cdfinv(at[i], FUB, YS01) ) )
	}

	// DONE.  RETURN STRUCTURE W/ FOUR ROW VECTORS.
	// Each vector has mean estimate in first column, plus one column for each element of "at"
	return(result)
} // end of cic


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


// FOR BOOTSTRAPPING, DRAW RANDOM SAMPLE WITH REPLACEMENT
real colvector drawsmpl(real colvector x, |real colvector wgt)
{
	// Inputs: 1. Vector we're drawing from
	//         2. (Optional) Frequency weights
	// Output: Vector with a simple random sample
	real scalar N
	if ((args()==1)|(wgt==0)) {
		N = rows(x)
		return(x[ceil(runiform(N,1):*N),1])
	}
	else {
		real colvector exp
		real scalar i,j,w
		N = sum(wgt) // assume weights are integers (I check this once in the main CIC function--before calling this sub-routine hundreds of times)
		// I thought it might be faster to just "expand" the dataset, then draw from it
		exp=J(N,1,.)
		j=1
		for (i=1;i<=rows(x);i++) {
			for (w=1;w<=wgt[i];w++) {
				exp[j]=x[i]; j++
			}
		}
		return(exp[ceil(runiform(N,1):*N),1])
	}
}


// FOR BOOTSTRAPPING, USE 95 PERCENTILES OF BOOTSTRAP ITERATIONS TO BACKOUT STANDARD ERRORS
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

YS
drawsmpl(YS)
drawsmpl(YS,(10\1\1\1\0))

mata describe
mata memory


end
/* * * * *  END OF MATA BLOCK * * * * */


cd "C:\Users\keith\Desktop\CIC\"


* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* test cic_vce_parse
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
cic_vce_parse, vce(boot, reps(1000) saving(myfile.dta, replace))
return list
cap nois cic_vce_parse, vce(none, reps(25))
return list
cic_vce_parse, vce(bspctile, reps(25) mse accel(myvector) saving("c:\temp\test.dta", replace))
return list


* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
set tracedepth 1
if 0  set trace on
else  set trace off
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
local Nreps = 50
if 0   	local vce vce(bootstrap, reps(`Nreps'))
else if 0	local vce vce(bspctile , reps(`Nreps'))
else       	local vce " "
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
		using A_I_Matlab\mvd.dat

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

cap log close
log using cid_test_aid_data.log, replace

mac list _Nreps _vce

// Table 1
count
tabstat y ly , by(high_after) s(count mean sd min p25 p50 p75 p90 max) columns(s)  labelwidth(30) nototal format(%9.2f)

// DID estimate
reg y high##after
reg ly high##after

// CIC estimates from A&I Appendix
// Table 2
timer on 2
cic  y high after ,  at(25 50 75 90) `vce' 
cic ly high after ,  at(50)          `vce'
timer off 2
timer list 2

// Table 3
timer on 3
cic  y high after ,  at(25 50 75 90) `vce' untreated
cic ly high after ,  at(50)          `vce' untreated
timer off 3
timer list 3

log close

// compare vce() option above to the bootstrap prefix
timer on 10
cic y high after, vce(bootstrap, reps(`Nreps'))
timer off 10
timer list 10
est store X11
estat bootstrap , all // estat bootstrap does work

set trace on
estimates replay X11
exit
timer on 12
bootstrap , reps(`Nreps') strata(high after) : cic y high after
timer off 12
timer list 12


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

exit

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


exit



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// NOTE:
// The code in cic() is somewhat convoluted because I am
// calculating all four vectors simultaneously.
// The alternative routines, provided here, are much
// more clear and accessible. These alternatives produce
// identical results to cic(). However redundant calculations
// lead to slower runtimes.
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

version 11.2
mata:
mata clear
mata set matastrict on
mata set matafavor speed

struct cic_result scalar cic_seperate(real colvector Y00, real colvector Y01, real colvector Y10, real colvector Y11, real rowvector at, | real colvector W00, real colvector W01, real colvector W10, real colvector W11 )
{
	// Inputs:
	//   (1)-(4) Four column vectors with data.
	//            - Y00 is control group in pre-period
	//            - Y01 is control group in post period
	//            - Y10 is treatment group in post period
	//            - Y11 is treatment group in post period
	//   (5)     Vector with k>=1 quantiles of interest, ranging from 0 to 1
	//   (6)-(9) (Optional) Column with fweights or iweights for Y00, Y01, Y10, and Y11 (respectively)
	//
	// Output: One structure (cic_result) with four row vectors.
	//   Each vector has (1+k) elements. The first element is the mean, followed by k results (one for each quantile in -at-).
	//   (1) result.con       = CIC ESTIMATOR WITH CONTINUOUS OUTCOMES, EQUATION 9
	//   (2) result.dci       = CIC MODEL WITH DISCRETE OUTCOMES (UNDER THE CONDITIONAL INDEPENDENCE ASSUMPTION), EQUATION 29
	//   (3) result.dcilowbnd = LOWER BOUND ESTIMATE OF DISCRETE CIC MODEL (WITHOUT CONDITIONAL INDEPENDENCE), EQUATION 25
	//   (4) result.dciuppbnd = UPPER BOUND ESTIMATE OF DISCRETE CIC MODEL (WITHOUT CONDITIONAL INDEPENDENCE), EQUATION 25

	// Need all or none of args (6)-(9)
	if (args()>5 & args()!=9) _error(( "Expected 5 or 9 arguements, but received " + strofreal(args())))

	// Declare variables
	real colvector YS, YS01
	real colvector F00, F01, F10, F11
	struct cic_result scalar result

	// Vector with support points for all four groups combined (YS) and the comparison-post group (YS01)
	YS01 = uniqrows(Y01)
	YS   = uniqrows(Y00\YS01\Y10\Y11)

	// Vector with CDF functions of the four treat*post groups (F00,F01,F10,F11)
	// CDFs (w/ and w/out weights declared)
	if (args()==5) {
		// CDFs without weights
		F00=runningsum(prob(Y00,YS))
		F01=runningsum(prob(Y01,YS))
		F10=runningsum(prob(Y10,YS))
		F11=runningsum(prob(Y11,YS))
	}
	else {
		// CDFs with weights
		F00=runningsum(prob(Y00,YS,W00))
		F01=runningsum(prob(Y01,YS,W01))
		F10=runningsum(prob(Y10,YS,W10))
		F11=runningsum(prob(Y11,YS,W11))
		// because of rounding, sum of weights might be slightly different than one
		F00[length(F00)]=1
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

	// DONE.  RETURN STRUCTURE W/ FOUR ROW VECTORS.
	// Each vector has mean estimate in first column, plus one column for each element of "at"
	return(result)
}



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

// Confirm cic() == cic_seperate()
struct cic_result scalar result1, result2
if (args()==6) result1=cic(Y00,Y01,Y10,Y11,at)                 // without weights
else           result1=cic(Y00,Y01,Y10,Y11,at,W00,W01,W10,W11) // with weights

if (args()==6) result2=cic_seperate(Y00,Y01,Y10,Y11,at)                 // without weights
else           result2=cic_seperate(Y00,Y01,Y10,Y11,at,W00,W01,W10,W11) // with weights

if (result1==result2) "Elements are equal"
else {
	result1.con;result1.dci;result1.dcilowbnd;result1.dciuppbnd;result2.con;result2.dci;result2.dcilowbnd;result2.dciuppbnd
	_error( "Elements not equal")
}

end
