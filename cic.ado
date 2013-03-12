clear all

*! *! version $Id$
*! Changes-in-changes
*!
*! An implimentation of:
*! Athey, S. & G. W. Imbens. "Identification and Inference in Nonlinear Difference-in-Differences Models."
*!     Econometrica, 74 (2), March 2006, pp. 431-497.
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
*                                  This is implimented according to parametric approach outlined in section 5.1. of AI.
*                                  "apply the CIC estimator to the residuals from an ordinary least squares regression with the effects of the dummy variables added back in." (p. 466)
*
*  and options are as follows:
*
*  MAIN
*    at(numlist)              a list of percentiles for CIC results. default is at(10(10)90)
*    vce(none|                don't calculate standard errors, the default
*        delta|               use numerical
*        bootstrap[, bsopts]| use bootstrap (by default, 1000 reps stratified by treat/post) other options allowed
*    did                      calculated traditional DID and quantile DID (always on if there are any control variables)
*    untreated                counterfactual effect of the policy for the untreated group (Setion 3.2 of paper)
*    round(integer)           round dependent variable to nearest r (=0 for no rounding, the default)
*                             this rounding is performed after adjusting for covariates, if applicable
*
*  REPORTING
*      level(passthru)              set confidence level; default is level(95)
*      notable                      suppress table of results
*      noheader                     suppress table header
*      nolegend                     suppress table legend
*      display_options              control spacing and display of omitted variables and base and empty cells
*                                      display_options:  noomitted, vsquish, noemptycells, baselevels, allbaselevels;
*                                                        see help estimation options##display_options
*
*  BSOPTS SUBOPTIONS
*     reps(#)                      perform # bootstrap replications; default is reps(200)
*     saving(filename[,replace])   save bootstrap results to filename (optionally, replace specifies that filename be overwritten, if it exists.)
*     sepercentile                 obtain bootstrap standard errors from percentiles of bootstrap estimates instead of using Stata's default method.
*                                      standard error = (p(97.5) - p(2.5)) / (2*1.96), where p(N) is Nth percentile of bootstrap iterations.
*                                      this is the method used in Athey and Imbens' MATLAB code.
*     accel(vector)                acceleration values for each statistic
*     mse                          use MSE formula for variance estimation
*     nodots                       suppress the replication dots
*     pop(#)                       total sample size (used for scaling iweights for bootstrap replications)
*                                      the sample in each group is calculated as the sum of the iweights for observations in the group, divided by the sum of the iweights for all observations, and multiplied by the population size
*  See [R] bootstrap postestimation for features available after estimation.


* Weights may be iweights or fweights.

* vce(bootstrap, [bsopts]) is equivalent to
*    . bootstrap _b, strata(treat post) [bsopts]: cic y treat post ... , vce(none)
* but slower because vce(bootstrap) is implimented in META and runs with less overhead.
* However, the bootstrap prefix is more flexible due the availability of size(), strata(), cluster(), idcluster() and other options.
* in documentation, talk about why bootstrap: prefix might be needed (longitudinal data, sample by id instead of pre/post groups)

* When by is used, only the last group is saved in ereturn.



* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* TO DO
* 1. SE option
* 2. QDID - add labels and return to e(b) and 
* 3. Add arguments so you can run only some of the
*    options (e.g., continuous only or discrete or qreg)
* 4. byable(recall) working right?
* 5. add error/documentation that *  fweights, but not iweights, work with vce( ??????????????? )
* 6. look at qreg and make similar matrix labels
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *



program cic, properties(mi) eclass byable(onecall)
	version 11.2
	if replay() {
		if ("`e(cmd)'"!="cic") error 301
		if _by()               error 190
		else Replay `0'
		exit
	}
	if _by() local BY `"by `_byvars'`_byrc0':"'
	`BY' Estimate `0'
	ereturn local cmdline `"cic `0'"'
 	ereturn local cmd     "cic"
	ereturn local title   "Changes in Changes (CIC) Model"
end

//"

program define Estimate, eclass byable(recall)
	version 11.2

	// parse arguments
	syntax varlist(default=none min=3 numeric fv ts) [if] [in] [fweight iweight]  ///
		[, at(numlist min=1 >=0 <=100 sort) ///
		Vce(passthru) ///
		did ///
		Qdid ///
		UNTreated ///
		ROUnd(real 0) ///
		level(passthru) notable NOHeader NOLegend * ] // Reporting options
	marksample touse  // note that rows are dropped for (1) if/in (2) zero weight (3) missing data (and other reasons, see "help mark")
	_get_diopts diopts, `options'
	local diopts `diopts' `table' `header' `legend'
	_rmcoll `varlist' [`weight'`exp'] if `touse', expand
	local varlist `r(varlist)'

	// first three variables need to be y, treat, and post
	gettoken y     varlist : varlist
	gettoken treat varlist : varlist
	gettoken post  varlist : varlist

	// parse percentiles
	if mi("`at'") local at "10(10)90" // default set (if undeclared)
	numlist "`at'"
	local at = r(numlist)

	// parse vce()
	if mi(`"`vce'"') local vce vce(none)
	cic_vce_parse, `vce'
	local vce     = r(vce)
	local bsreps  = r(bsreps)
	local bsiwpop = r(bsiwpop)
	local dots    = r(dots)
	if !missing(r(bsopts))       local bsopts       = r(bsopts)
	if !missing(r(saving))       local saving       = r(saving)
	local sepercentile = (r(sepercentile)==1)

	// prep to handle weights
	if !missing("`weight'") {
		tempvar wgtvar
		gen `wgtvar'`exp' if `touse'
		summ `wgtvar' if `touse' , meanonly
		local n=round(r(sum))
		if "`wgtvar'"=="iweight" {
			if ("`vce'"=="bootstrap" & `bsiwpop'==0) {
				di as error "If bootstrapping standard errors with iweights, you must use the vce(bootstrap, pop(N)) suboption to declare the population size."
				error 198
			}
			local wtexp_caller = `", "`wgtvar'", `bsiwpop' "'
		}
		else local wtexp_caller = `", "`wgtvar'", `bsiwpop' "'
	}
	else {
		qui count if `touse'
		local n=r(N)
	}


	// adjust y for covariates (OLS regression)
	local runDID  = (("`did'"=="did") | (`: list sizeof varlist'!=0))

if ( `runDID' ) regress `y' ibn.`treat'#ibn.`post' `varlist' if `touse' [`weight'`exp'], `diopts' `level' nocons

	// implement mata CIC routine
	ereturn clear
	tempname mata_b
	if "`untreated'"=="" {
		// default - effect of treatment on the treated
		mata: cic_caller( "`y' `treat'   `post' `varlist'", "`touse'", "at", `runDID', 1, `bsreps', `dots', `round' `wtexp_caller')
	}
	else {
		// option - effect of treatment on the untreated (just switch `treat' variable and use -e(results))
		tempvar untreat
		gen `untreat' = (`treat'==0)
		mata: cic_caller( "`y' `untreat' `post' `varlist'", "`touse'", "at", `runDID', 0, `bsreps', `dots', `round' `wtexp_caller')
	}

	// post results to ereturn
	di as txt _n "Changes in Changes (CIC) Model"
	if (`bsreps'>0 & !mi(`bsreps')) {
		// option - save bs reps into a specified file (`saving' contains ", replace" if it was provided)
		if !mi(`"`saving'"') copy `bstempfile' `saving'

		// post boostrapped estimates with standard errors using bstat
		`=cond(`sepercentile',"quietly","")' /// quietly if need to call sepercentile
		bstat  using `bstempfile', stat( `mata_b' ) `bsopts' `diopts' `level'

		if ( `sepercentile' ) {
			if (e(level)!=95) {
				ereturn display
				di as error "-sepercentile- sub-option only works with level(95).  Standard errors displayed above are from Stata's default method of producing standard errors."
				error 198
			}
			if `bsreps' < 1000 nois di as error "Warning: More bootstrap repetitions might be needed with the -sepercentile- bootstrapping sub-option.."
			tempname ci_se
			mata : bs_se( "e(ci_percentile)", "`ci_se'" )
			ereturn repost V = `ci_se'
			ereturn display
		}
	}
	else {
		// otherwise just use ereturn to get pretty estimates table
		di as txt _col(49) "Number of obs      =" as res %8.0fc = e(N)
		ereturn post `mata_b' [`weight'`exp'], depname(`y') obs(`n') esample(`touse') dof(`=`n'-colsof(`mata_b')') `level'
		ereturn display
	}
 	ereturn local depvar  "`y'"
 	ereturn local vce     "`vce'"
	if (`: list sizeof varlist'!=0 | `runDID') {
		ereturn scalar k_eq =  5
		ereturn local  eqnames continuous discrete_ci dci_lower_bnd dci_upper_bnd did
	}
	else {
	 	ereturn scalar k_eq =  4
		ereturn local  eqnames continuous discrete_ci dci_lower_bnd dci_upper_bnd
	}
	if "`untreated'"=="" ereturn local footnote "Effect of Treatment on the Treated Group"
	else                 ereturn local footnote "Effect of Treatment on the Untreated Group"
	di as txt "(" e(footnote) ")"

mata: mata describe

end // end of cic program definition


// subroutine to replay estimates
// this section is similar in function to the "_vce_parse" command, except that I set default values for reps() and strata()
program Replay
	syntax [, notable noHeader  noRULES OR GROUPED *]
	_get_diopts diopts, `options'
	local diopts `diopts' `table' `header' `legend'
	_prefix_display, `diopts' `table' `header' `legend'
	di as txt "(" e(footnote) ")"
end

// subroutine to parse the vce() option
// this section is similar in function to the "_vce_parse" command, except that I set default values for reps() and strata()
program define cic_vce_parse, rclass
	version 11.2
	syntax , vce(string asis)
	_parse comma vce 0 : vce
	if  inlist( "`vce'","bootstra","bootstr","bootst","boots","boot","boo","bo","b") local vce "bootstrap"
	if  inlist( "`vce'","delt","del","de","d")                                       local vce "delta"
	if  inlist( "`vce'","non","no","n")                                              local vce "none"

	if ("`vce'"!="bootstrap") & !mi("`0'") {
		di as error "suboptions are not allowed with vce(`vce')"
		error 198
	}

	return local vce `vce'
	if ("`vce'"=="bootstrap") {
		syntax [, Reps(integer 200) pop(integer 0) SAving(string asis) NODots SEPercentile *]
		return scalar bsreps  = `reps'
		return scalar bsiwpop = `pop'
		return scalar dots  = ( "nodots"!="`dots'" )
		return scalar sepercentile = ( "sepercentile"=="`sepercentile'" )
		return local  saving  : copy local saving
		return local  bsopts  : copy local options
	}
	else if ("`vce'"=="delta") {
		return scalar bsreps  = -1
		return scalar dots    = 0
	}
	else if ("`vce'"=="none") {
		return scalar bsreps  = 0
		return scalar dots    = 0
	}
	else {
		di as error "Only vce(delta), vce(bootstrap [, subopts]), and vce(none) allowed."
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

// STRUCTURES FOR RETURNING RESULTS
struct cic_result {
	real rowvector con, dci, dcilowbnd, dciuppbnd
}
struct did_ols_result {
	pointer(real colvector) Y
	real colvector coef
	real scalar did
	string matrix labels, didlabel
}


// CIC CALLER -- THIS FUNCTION READS STATA DATA INTO MATA AND CALLS THE MAIN CIC ROUTINE
void cic_caller(string rowvector varlist, string scalar touse_var, string scalar at_local, real scalar did, real scalar tot, real scalar bsreps, real scalar dots, real scalar round, |string scalar wgt_var, real scalar popsize)
{
	// Inputs:

	//   1.   Name of variables, in the following order:
	//         - dependent variable, `y'
	//         - treatment dummy (0/1 for control/treatment groups)
	//         - time period dummy (0/1 for pre/post periods)
	//         - (Optional) control variables
	//   2.   Name of variable indicating which rows to include, `touse'
	//   3.   Name of local macro containing quantiles of interest, ranging from 0 to 100 percent (`at')
	//   4.   0/1 to estimate difference in difference (DID) and quantile DID model.  (always estimated if there are control variables)
	//   5.   0/1 if estimates effect of treatment on the treated (=1) or effect of treatment on the untreated (=0)
	//               if untreated, you must provide a variable that =0 if treated and =1 if untreated
	//               setting tot==0 simply flips the sign of the cic() output
	//   6.   How to calculate standard errors:
	//         - if >1, bootstrapped standard errors and bsreps = number of bootstrap repetitions)
	//         - if 0, no standard errors
	//         - if =(-1), standard error for conditional independence based on numerical differentiation
	//   7.   0/1 to hide bootstrapping dots (=1 to show dots)
	//   8.  Round y to the nearest ___.  (set =0 for no rounding).
	//   (Optional)
	//   9.  Name of variable with fweight or iweight
	//   10. Population size for scaling iweights.  Set to zero if fweights.  (This is only used if bootstrapping SEs.)
	// Output: Output is returned to stata through various st_*() functions.

	// Read y, treat and post into mata
	real colvector y, treat, post
	varlist = tokens(varlist)
	st_view(y    =.,.,varlist[1],touse_var)  // note that rows with missing data are already dropped by -marksample- in .ado file
	st_view(treat=.,.,varlist[2],touse_var)
	st_view(post =.,.,varlist[3],touse_var)
	if ((uniqrows(treat),uniqrows(post))!=(0,0\1,1)) {
		_error( "treat and post must be dummy variables equal to 0 and 1" )
	}
	real scalar N
	N   = rows(y)

	// read control variables into mata (if need to run DID model)
	if (length(varlist)>3 | did ) {
		did = 1  // always run DID regression if control variables present
		real matrix rhs
		st_view(rhs =.,.,invtokens(varlist[2..length(varlist)]),touse_var)  // note that rows with missing data are already dropped by -marksample- in .ado file
	}

	// read weights into mata (if there are any)
	if (args()!=8) {
		real colvector wgt
		st_view(wgt=.,.,wgt_var,touse_var)
		if (args()==9) popsize=0
	}
	else wgt=1

	// Quantiles
	real rowvector at
	at = strtoreal(tokens(st_local(at_local))) / 100
	if (min(at)<0 | max(at)>1) _error( "at() must be between 0% and 100% (inclusive)." )

	// Results will be returned into a structures (defined above)
	struct cic_result scalar result
	struct did_ols_result scalar didresult

	// DID regression
	// After this section, upper-case Y is now the dependent variable for
	// the cic() function. It is a pointer. It points to (lower case) y
	// or a (temporary) variable that is adjusted for covariates and/or rounded.
	pointer(real colvector) scalar Y
	Y = &y
	if (did) {
		didresult = did_ols(y, rhs, wgt, round, varlist)
		swap(Y,didresult.Y)
	}
	else if (round!=0) Y = &round(y,round)

"sizeof(y)"
sizeof(y)
"sizeof(Y)"
sizeof(Y)
"sizeof(*Y)"
sizeof(*Y)
if (did) "(*didresult.Y)[1..10]"
if (did)  (*didresult.Y)[1..10]
"(*Y)[1..10]"
 (*Y)[1..10]
"rows in y"
rows((y))
"rows in *Y"
rows((*Y))
"unique rows in y"
rows(uniqrows(y))
"unique rows in *Y"
rows(uniqrows(*Y))

	// Permutation vectors identifying the four treat*post groups
	real colvector p00, p01, p10, p11
	st_select(p00=.,(1::N),(treat:==0 :& post:==0))
	st_select(p01=.,(1::N),(treat:==0 :& post:==1))
	st_select(p10=.,(1::N),(treat:==1 :& post:==0))
	st_select(p11=.,(1::N),(treat:==1 :& post:==1))

	// Number of observations
	real scalar N00, N01, N10, N11
	N00 = rows(p00);
	N01 = rows(p01)
	N10 = rows(p10)
	N11 = rows(p11)
	if (min((N00,N01,N10,N11))<1) _error( "One or more of the four treat*post groups is empty.")
	if (min((N00,N01,N10,N11))<2 & bsreps>0) _error( "One or more group has size less than 2. There will be no variation in bootstrap draws.")

	// Call the quantile DID routine
	if (did) {
		real rowvector qdid_result
		if (args()==8) qdid_result=qdid((*Y)[p00],(*Y)[p01],(*Y)[p10],(*Y)[p11],at) // without weights
		else           qdid_result=qdid((*Y)[p00],(*Y)[p01],(*Y)[p10],(*Y)[p11],at,wgt[p00],wgt[p01],wgt[p10],wgt[p11]) // with weights
	}
"Output from qdid is saved in qdid_result"
"Need to add labels and return qdid_result to e(b)"
qdid_result


//	// Select the rows belonging to the treat*post groups
//	real colvector Y00, Y01, Y10, Y11
//	st_select(Y00=.,*Y,(treat:==0 :& post:==0))
//	st_select(Y01=.,*Y,(treat:==0 :& post:==1))
//	st_select(Y10=.,*Y,(treat:==1 :& post:==0))
//	st_select(Y11=.,*Y,(treat:==1 :& post:==1))
//
//	// Select the rows of wgt belonging to the treat*post groups
//	if (args()!=8) {
//		real colvector W00, W01, W10, W11
//		st_select(W00=.,wgt,(treat:==0 :& post:==0))
//		st_select(W01=.,wgt,(treat:==0 :& post:==1))
//		st_select(W10=.,wgt,(treat:==1 :& post:==0))
//		st_select(W11=.,wgt,(treat:==1 :& post:==1))
//	}
"Begin Main CIC call:"

	// Call the main CIC routine
	if (args()==8) result=cic((*Y)[p00],(*Y)[p01],(*Y)[p10],(*Y)[p11],at) // without weights
	else           result=cic((*Y)[p00],(*Y)[p01],(*Y)[p10],(*Y)[p11],at,wgt[p00],wgt[p01],wgt[p10],wgt[p11]) // with weights

	// return results into a Stata matrix named st_local("mata_b") with lables
	if (did) {
		if (tot) st_matrix(st_local("mata_b"),  (didresult.coef',didresult.did,result.con,result.dci,result.dcilowbnd,result.dciuppbnd))
		else     st_matrix(st_local("mata_b"), -(didresult.coef',didresult.did,result.con,result.dci,result.dcilowbnd,result.dciuppbnd))
	}
	else {
		if (tot) st_matrix(st_local("mata_b"),  (result.con,result.dci,result.dcilowbnd,result.dciuppbnd))
		else     st_matrix(st_local("mata_b"), -(result.con,result.dci,result.dcilowbnd,result.dciuppbnd))
	}

	// matrix labels for `mata_b'
	string matrix colfulllabels
	colfulllabels=((J(1+length(at),1,"continuous") \ J(1+length(at),1,"discrete_ci") \ J(1+length(at),1,"dci_lower_bnd") \ J(1+length(at),1,"dci_upper_bnd")),J(4,1,strtoname(("mean" , ("q":+strofreal(at*100))))'))
	if (did) colfulllabels = (didresult.labels \ didresult.didlabel \ colfulllabels )
	st_matrixcolstripe(st_local("mata_b"), colfulllabels)
	st_local("cic_coleq"   ,invtokens(colfulllabels[.,1]'))
	st_local("cic_colnames",invtokens(colfulllabels[.,2]'))

	// return
	st_matrix("e(at)",at)
	st_numscalar( "e(N_strata)", 4)
	if (args()==8) {
		st_numscalar( "e(N)"       , N)
		st_numscalar( "e(N00)"     , N00)
		st_numscalar( "e(N01)"     , N01)
		st_numscalar( "e(N10)"     , N10)
		st_numscalar( "e(N11)"     , N11)
	}
	else if (popsize) {
		st_numscalar( "e(N)"       , round(popsize))
		st_numscalar( "e(N_obs)"   , N)
	}
	else {
		st_numscalar( "e(N)"       , round(sum(wgt)))
		st_numscalar( "e(N_obs)"   , N)
	}
	st_numscalar( "e(N_support)",rows(uniqrows(*Y)))

"start BS now:"
	// Bootstrapping
	if (bsreps>0) {
		real scalar b
		real colvector bs_wgt
		struct did_ols_result scalar bs_didresult
		struct cic_result     scalar bs_cicresult

		// pointer to y
		// a new pointer is needed for dependent variable since it might be adjusted for covariates with boostrap sample
		pointer(real colvector) scalar bs_Y
		bs_Y = Y

		// empty matrix to store results
		real matrix bsdata
		if (did) bsdata=J(bsreps,4*(1+length(at))+length(didresult.coef)+1,.)
		else     bsdata=J(bsreps,4*(1+length(at)),.)

		// Before loop, extra setup needed for drawing a sample with unequal weights
		if (args()!=8) {
			// weight variables with cumulative sum of the weights from each group
			real colvector cumsum00, cumsum01, cumsum10, cumsum11
			cumsum00 = quadrunningsum(wgt[p00])
			cumsum01 = quadrunningsum(wgt[p01])
			cumsum10 = quadrunningsum(wgt[p10])
			cumsum11 = quadrunningsum(wgt[p11])

			// scalars with population size for each group
			real scalar popsize00, popsize01, popsize10, popsize11
			if (popsize) {
				// With importance weights, use the fraction of popsize (e.g., popsize00 = round(cumsum00[n00]/colsum(wgt)*popsize))
				real scalar sumwgt
				sumwgt = quadcolsum(wgt)
				popsize00 = round(cumsum00[N00]/sumwgt*popsize) // the number of obs. in each group is rounded to the nearest integer
				popsize01 = round(cumsum01[N01]/sumwgt*popsize)
				popsize10 = round(cumsum10[N10]/sumwgt*popsize)
				popsize11 = round(cumsum11[N11]/sumwgt*popsize)
			}
			else {
				// With frequency weights, this is the weighted number of individuals in the group (e.g., popsize00 = cumsum00[n00])
				popsize00 = cumsum00[N00]
				popsize01 = cumsum01[N01]
				popsize10 = cumsum10[N10]
				popsize11 = cumsum11[N11]
				if (round(popsize00)!=popsize00 | round(popsize01)!=popsize01 | round(popsize10)!=popsize10 | round(popsize11)!=popsize11) "When drawing bootstrap sample with frequency weights, found non-integer fweights in one or more groups."
			}
			cumsum00 = cumsum00/cumsum00[N00] // normalize to sum to one within groups
			cumsum01 = cumsum01/cumsum01[N01]
			cumsum10 = cumsum10/cumsum10[N10]
			cumsum11 = cumsum11/cumsum11[N11]
			if (min((popsize00,popsize01,popsize10,popsize11))<2) "One or more groups has size less than 2. There will be no variation in bootstrap draws."
		}

		// header for dots
		if (dots) {
			printf( "{txt}\nBootstrap replications ({res}%g{txt})\n", bsreps)
			display( "{txt}{hline 4}{c +}{hline 3} 1 " +
				"{hline 3}{c +}{hline 3} 2 " + "{hline 3}{c +}{hline 3} 3 " +
				"{hline 3}{c +}{hline 3} 4 " + "{hline 3}{c +}{hline 3} 5 ")
		}

		// Bootstrapping replications
		for(b=1; b<=bsreps; ++b) {

			if (args()!=8 | did) {
if (b==1) "weights"
				// if estimating DID model or have a weighted sample, the bootstrap sample
				// is drawn by calculating a new weigting vector.

				// Draw bootstrap sample
				// bs_draw_wgt() produces a vector with frequency weights in the unweighted case
				// or a replacement weight vector (iweights or fweights) in the weighted case

				if (args()==8) bs_wgt = bs_draw_wgt(p00, p01, p10, p11, N00, N01, N10, N11)
				else           bs_wgt = bs_draw_wgt(p00, p01, p10, p11, N00, N01, N10, N11, cumsum00, cumsum01, cumsum10, cumsum11, popsize00, popsize01, popsize10, popsize11)

				// calculate DID and adjust for covariates w/ bootstrap sample
				if (did) {
					bs_didresult = did_ols(y, rhs, bs_wgt, round, varlist)
					swap(bs_Y,bs_didresult.Y)
				}

				// call cic() with bootstrap sample
				bs_cicresult=cic((*bs_Y)[p00],(*bs_Y)[p01],(*bs_Y)[p10],(*bs_Y)[p11],at,bs_wgt[p00],bs_wgt[p01],bs_wgt[p10],bs_wgt[p11])
			}
			else {

if (b==1) "no weights"
				// in the simple case of no regression adjustment and no weights, simply
				// call cic() with a random draw of dependent variable
				bs_cicresult=cic(bs_draw_nowgt((*bs_Y)[p00]),bs_draw_nowgt((*bs_Y)[p01]),bs_draw_nowgt((*bs_Y)[p10]),bs_draw_nowgt((*bs_Y)[p11]),at)
			}

			// save estimates into a matrix with one row per bootstrap sample
			if (did) {
				if (tot==1) bsdata[b,.]  =  (bs_didresult.coef',bs_didresult.did,bs_cicresult.con,bs_cicresult.dci,bs_cicresult.dcilowbnd,bs_cicresult.dciuppbnd)
				else        bsdata[b,.]  = -(bs_didresult.coef',bs_didresult.did,bs_cicresult.con,bs_cicresult.dci,bs_cicresult.dcilowbnd,bs_cicresult.dciuppbnd)
			}
			else {
				if (tot==1) bsdata[b,.]  =  (bs_cicresult.con,bs_cicresult.dci,bs_cicresult.dcilowbnd,bs_cicresult.dciuppbnd)
				else        bsdata[b,.]  = -(bs_cicresult.con,bs_cicresult.dci,bs_cicresult.dcilowbnd,bs_cicresult.dciuppbnd)
			}

			// show dots
			if (dots) {
				if (missing(bsdata[b,.])) printf( "{err}x{txt}")
				else printf( ".")
				if (!mod(b,50)) printf( " %5.0f\n",b)
				else if (b==bsreps & mod(b-1,50)) display("") // end of dots
				displayflush()
			}
		} // end loop through bs iterations
"done with bootstrapping loosps"

		// save bootstrap iterations in a temporary .dta file (named `bstempfile')
		stata( "preserve" )
		  string rowvector bstempfile, bstempvars
		  // clear data (after preserve) and fill with bsdata matrix
		  st_dropvar(.)
		  st_addobs(rows(bsdata))
		  bstempvars=strtoname("_bs_":+strofreal(1::cols(bsdata)))'
		  st_store(.,st_addvar("double",bstempvars), bsdata)
		  // setup file for bstat command
		  st_global( "_dta[bs_version]" , "3")
		  if (args()==8)    st_global( "_dta[N]", strofreal(N))
		  else if (popsize) st_global( "_dta[N]", strofreal(popsize))
		  else              st_global( "_dta[N]", strofreal(round(sum(wgt))))
		  st_global( "_dta[N_strata]"   , "4")
		  st_global( "_dta[strata]"     , (varlist[2] + " " + varlist[1]))
		  st_global( "_dta[command]"    , "cic")
		  if (did) st_global( "_dta[k_eq]", "6")
		  else     st_global( "_dta[k_eq]", "4")
		  st_global( "_dta[k_extra]"    , "0")

		  for(b=1; b<=cols(bsdata); ++b) {
			if (did) {
				if (tot==1) st_global( (bstempvars[1,b]+"[observed]")  , strofreal( (bs_didresult.coef',bs_didresult.did,bs_cicresult.con,bs_cicresult.dci,bs_cicresult.dcilowbnd,bs_cicresult.dciuppbnd)[1,b]))
				else        st_global( (bstempvars[1,b]+"[observed]")  , strofreal(-(bs_didresult.coef',bs_didresult.did,bs_cicresult.con,bs_cicresult.dci,bs_cicresult.dcilowbnd,bs_cicresult.dciuppbnd)[1,b]))
			}
			else {
				if (tot==1) st_global( (bstempvars[1,b]+"[observed]")  , strofreal( (bs_cicresult.con,bs_cicresult.dci,bs_cicresult.dcilowbnd,bs_cicresult.dciuppbnd)[1,b]))
				else        st_global( (bstempvars[1,b]+"[observed]")  , strofreal(-(bs_cicresult.con,bs_cicresult.dci,bs_cicresult.dcilowbnd,bs_cicresult.dciuppbnd)[1,b]))
			}
			 st_global( (bstempvars[1,b]+"[expression]"), ( "["+colfulllabels[b,1]+"]_b["+colfulllabels[b,2]+"]"))
			 st_global( (bstempvars[1,b]+"[coleq]")     , colfulllabels[b,1])
			 st_global( (bstempvars[1,b]+"[colname]")   , colfulllabels[b,2])
			 st_global( (bstempvars[1,b]+"[is_eexp]")   , "1" )
		  }

		  // save as `bstempfile'
		  bstempfile=st_tempfilename()
		  st_local( "bstempfile",bstempfile)
		  stata(( "qui save " + bstempfile ))
		stata( "restore" )
	} // done bootstrapping
	else if (bsreps==-1) {
_error( "Code for sedelta not written." )
	}
	else if (bsreps==0) "Specify vce() option to calculate standard errors."
	else _error( "bsreps invalid.")
"end of cic_caller"


} // end of cic_caller; everthing is returned to Stata with st_*() commands.


// >>>>>>>>>>  check that column names and such are the same as the ouput in the example in the NOTE (below)

// CIC ROUTINE
struct cic_result scalar cic(real colvector Y00, real colvector Y01, real colvector Y10, real colvector Y11, real rowvector at, | real colvector W00, real colvector W01, real colvector W10, real colvector W11 )
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
	if (length(YS)<2) _error("The dependent variable is a constant")

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


// TRADITIONAL DIFFERENCES IN DIFFERENCES REGERSSION (OLS)
struct did_ols_result scalar did_ols(real colvector y, real matrix rhs, real colvector wgt, real scalar round, |string rowvector varlist)
{
	// Inputs:
	// (1) y, the dependent variable
	// (2) rhs, matrix of independent variables
	//      - first column is treat
	//      - second column is post
	//      - (optional) remaining columns are covariates
	// (3) wgt, a column vector of fweights or iweights (can set to scalar =1 for no weights)
	// (4) round, a scalar indicating the nearest unit for rounding Y (=0 for no rounding)
	// (5) (optional) vector with variable list (columns corresponding to names of (y,rhs)
	//
	// Output: One structure (did_ols_result) with:
	// 1. *Y, pointer to adjusted variable (pointing to either a temporary variable or input y itself)
	// 2. a vector with coefficients from the DID regression
	// 3. (if varlist provided) labels for coefficients compatible for st_matrixcolstripe()


	// Nx4 matrix with dummies indicating group membership to p00, p01, p10, p11 (respectively)
	real matrix D
	D = ( ((-rhs[.,1]:+1):*(-rhs[.,2]:+1)), ((-rhs[.,1]:+1):*(rhs[.,2])), ((rhs[.,1]):*(-rhs[.,2]:+1)), ((rhs[.,1]):*(rhs[.,2])))

	// OLS DID regression
	struct did_ols_result scalar didresult
	if (cols(rhs)==2) didresult.coef = invsym(quadcross(D,wgt,D))*quadcross(D,wgt,y)
	else              didresult.coef = invsym(quadcross((D,rhs[.,3..cols(rhs)]),wgt,(D,rhs[.,3..cols(rhs)])))*quadcross((D,rhs[.,3..cols(rhs)]),wgt,y)

	// DIFF-IN-DIFF estimate
	didresult.did = didresult.coef[4] - didresult.coef[2] - didresult.coef[3] + didresult.coef[1]

	// Dependent variable (potentially adjusted for covariates or rounded)
	if (cols(rhs)>2) {
		// adjust for covariates
		//     yadj = y - X * _b[X]
		//          = D * _b[D] + resid    (yadj is also rounded if round!=0)
		//
		// notice that control variables are columns 3 to cols(rhs) of the matrix rhs,
		// but are in rows 5 to rows(didresult.coef) of the regression's independent variables
		// because the treat/post dummies (the first two variables in rhs) were turned
		// into four group dummies in the regression
		didresult.Y = &round( y - rhs[.,3..cols(rhs)]*didresult.coef[5..rows(didresult.coef),1] , round)
	}
	else if (round!=0) {
		// no covariaters but need to round y
		didresult.Y = &round(y,round)
	}
	else {
		// no covariaters or rounding, just point to input vector
		didresult.Y = &y
	}

	// labels for didresult.coef
	if (args()==5) {
		if (cols(rhs)==2) didresult.labels = (J(rows(didresult.coef),1,"did_model"),( ( "0."+varlist[2]+"#0."+varlist[3]) \( "0."+varlist[2]+"#1."+varlist[3]) \( "1."+varlist[2]+"#0."+varlist[3]) \( "1."+varlist[2]+"#1."+varlist[3])))
		else              didresult.labels = (J(rows(didresult.coef),1,"did_model"),( ( "0."+varlist[2]+"#0."+varlist[3]) \( "0."+varlist[2]+"#1."+varlist[3]) \( "1."+varlist[2]+"#0."+varlist[3]) \( "1."+varlist[2]+"#1."+varlist[3]) \ varlist[4..length(varlist)]'))
	}
	// labels for didresult.coef
	didresult.didlabel = ( "did", "did" )

	return(didresult)
}


// QUANTILE DID MODEL, EQUATION 22
real rowvector qdid(real colvector Y00, real colvector Y01, real colvector Y10, real colvector Y11, real rowvector at, | real colvector W00, real colvector W01, real colvector W10, real colvector W11 )
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
	// Output: Vector with (k) elements. (one for each quantile in -at-).
	real rowvector qdid; qdid = J(1,length(at),.)
	real scalar i

	// Need all or none of args (6)-(9)
	if (args()>5 & args()!=9) _error(( "Expected 5 or 9 arguements, but received " + strofreal(args())))

	if (args()==5) {
		// No weights
		for(i=1; i<=length(at); ++i) {
			qdid[i] = cumdfinv(Y11,at[i])-cumdfinv(Y10:+mean(Y01):-mean(Y00),at[i])
		}
	}
	else {
		// With weights
		for(i=1; i<=length(at); ++i) {
			qdid[i] = cumdfinv(Y11,at[1,i],W11)-cumdfinv(Y10:+mean(Y01,W01):-mean(Y00,W00),at[1,i],W10)
		}
	}
	return(qdid)
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
		return(rowsum((abs((YS:-J(n,1,Y'))):<=epsilon(1)):*J(n,1,wgt')):/J(n,1,quadcolsum(wgt)))
	}
	else {
		// without weights
		return(rowsum(abs((YS:-J(n,1,Y'))):<=epsilon(1)):/length(Y))
	}
}


// CUMULATIVE DISTRIBUTION FUNCTION
real scalar cdf(real scalar y, real vector P, real vector YS)
{
	// given a cumulative distrubtion function (P) over the support points (YS),
	// returns the empirical cumulative distribution function at a scalar (y)
	if      (y< YS[1])          return(0)
	else if (y>=YS[length(YS)]) return(1)
	else                        return(P[colsum((YS:<=(y+epsilon(1))))])
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
	// given a cumulative distribution function (P) over the support points (YS),
	// returns the inverse of the empirical cumulative distribution function at probability p (0<p<1)
	// but if equals -oo, it returns min(YS)-100*(1+max(YS)-YS(min)) = 101*min(YS)-100*max(YS)-100
	if (p>=(P[1]-epsilon(1))) {
		return(YS[max(select((1::length(YS)),(P:<=(p+epsilon(1)))))])
	}
	else {
		return(101*YS[1]-100*YS[length(YS)]-100)
	}
}


// EMPIRICAL DISTRIBUTION FUNCTION
real scalar cumdfinv(real colvector X, real scalar p, |real colvector wgt)
{
	// given a vector of observations (X),
	// returns the empirical distribution of X evaluated at a point (p).
	// optionally, the vector X can have weights (wgt)
	if      (p<=epsilon(1))   return(min(X))
	else if (p>=1-epsilon(1)) return(max(X))
	else if (args()==2) {
		// without weights
		return(sort(X,1)[floor(length(X)*p+1-epsilon(1)),1])

		// Note that floor(length(X)*p+1-epsilon(1)) is smallest integer larger than length(X)*p
		// e.g., if length(X)=10, p=0.34 then floor(3.4+1-2.2e-16)=4
		//       if length(X)=10, p=0.30 then floor(3.0+1-2.2e-16)=3
	}
	else {
		// with weights
		real matrix xs, sum_wgt
		xs = sort((X,wgt),1)
		sum_wgt = runningsum(xs[.,2]) :/ colsum(xs[.,2])
		sum_wgt[rows(xs),1]=1 // total of weight column is sum_wgt[rows(xs),1], set to 1

		// return the observation from fist row
		// where sum_wgt larger than p
		return(xs[colmax(select((1::rows(xs)),(sum_wgt:<=p))),1])

		// The weighted e.c.d.f. (empirical cumulative distribution function) Fn is defined so that,
		// for any real number y, the value of Fn(y) is equal to the total weight of all entries of x
		// that are less than or equal to y. That is Fn(y) = sum(weights[x <= y]).
		// Thus Fn is a step function which jumps at the values of x. The height of the jump at a point
		// y is the total weight of all entries in x number of tied observations at that value.
		// cumdfinv() returns the inverse of Fn()
	}
}


// FOR BOOTSTRAPPING, DRAW RANDOM SAMPLE WITH REPLACEMENT
real colvector bs_draw_nowgt(real colvector x)
{
	// Input:  Vector we're drawing rows from
	// Output: Vector with a simple random sample
	// (This function is adequate when x is a vector,
	// and when it is a permutatin vector)
	return(x[ceil(runiform(rows(x),1):*rows(x)),1])
}


// FOR BOOTSTRAPPING, DRAW RANDOM SAMPLE WITH REPLACEMENT
real colvector bs_draw_wgt(real colvector      p00, real colvector      p01, real colvector      p10, real colvector      p11,
                           real scalar         N00, real scalar         N01, real scalar         N10, real scalar         N11,
                         | real colvector cumsum00, real colvector cumsum01, real colvector cumsum10, real colvector cumsum11,
                           real scalar   popsize00, real scalar   popsize01, real scalar   popsize10, real scalar   popsize11)
{
	// Case 1: Unweighted
	// Inputs: 1-4.   Four (4) permutation vectors identifying the rows of data belonging to each group
	//         5-8.   Four (4) scalars with the number of observations in each group (e.g., N00=rows(p00))
	// Output: A frequency weight vector
	//
	// Case 2: Freqency or Importance Weights (fweights or iweights)
	// Inputs: 1-8.   All the inputs provided in Case 1 (unweighted)
	//         9-12.  Four (4) weight variables with cumulative sum of the weights from each group, normalized to sum to one within groups (e.g., cumsum00 = quadrunningsum(wgt[p00]); cumsum00 = cumsum00/cumsum00[N00])
	//         13-16. Four (4) scalars with population size for each group.
	//                  - With frequency weights, this is the weighted number of individuals in the group (e.g., popsize00 = cumsum00[n00])
	//                  - With importance weights, use the fraction of popsize (e.g., popsize00 = round(cumsum00[N00]/colsum(wgt)*popsize))
	// Output: A weight vector that replaces the input vector in cic_caller (wgt)
	//
	// This specialzed program was written to minimize processing time for CIC bootstrapping.
	// The basic principle is that anything calculated more than once should be caclulated
	// just once in cic_caller(), leaving only the tasks needed for each draw.
	// Case 2 code was modeled on mm_upswr() in the moremata package (version 1.0.0 by Ben Jann).
	real colvector u, reweight
	real scalar i, j
	reweight=J(N00+N01+N10+N11,1,0)

	if (args()==8) {
		// no weights, just a simple random draw from each group

		// random draw of permutation vectors
		u  = ( p00[ceil(runiform(N00,1):*N00),1] \
		       p01[ceil(runiform(N01,1):*N01),1] \
			   p10[ceil(runiform(N10,1):*N10),1] \
			   p11[ceil(runiform(N11,1):*N11),1] )

		// count number of times each row was drawn
		for (i=1;i<=rows(u);i++) {
			reweight[u[i]] = reweight[u[i]]+1
		}
	} // end unweighted section

	else if (args()==16) {
		// fweights or iweights
		real colvector r

		// 1st group
		u = sort(runiform(popsize00,1),1)   // random draw
		r = J(N00,1,0)
		j=1
		for (i=1;i<=popsize00;i++) {
			while (u[i]>cumsum00[j]) j++    // use cumulative distribution to get counts
			r[j] = r[j]+1                   // r will contain the number of observations drawn for each row
		}
		reweight[p00] = r                   // return results for this group

		// 2nd group
		u = sort(runiform(popsize01,1),1)
		r = J(N01,1,0)
		j=1
		for (i=1;i<=popsize01;i++) {
			while (u[i]>cumsum01[j]) j++
			r[j] = r[j]+1
		}
		reweight[p01] = r

		// 3rd group
		u = sort(runiform(popsize10,1),1)
		r = J(N10,1,0)
		j=1
		for (i=1;i<=popsize10;i++) {
			while (u[i]>cumsum10[j]) j++
			r[j] = r[j]+1
		}
		reweight[p10] = r

		// 4th group
		u = sort(runiform(popsize11,1),1)
		r = J(N11,1,0)
		j=1
		for (i=1;i<=popsize11;i++) {
			while (u[i]>cumsum11[j]) j++
			r[j] = r[j]+1
		}
		reweight[p11] = r

	} // end weights section
	else _error( "Expecting 8 arguments (unweighted) or 16 arguments (weighted), but received " + strofreal(args()) )
	return(reweight)
}


// AFTER BSTAT STATA PROGRAM, USE 95 PERCENTILES OF BOOTSTRAP ITERATIONS TO BACK-OUT STANDARD ERRORS
void bs_se( string scalar in_ci, string scalar out_V )
{
	// Inputs: 1. Vector with
	//         2. Stata matrix name for output
	// Output: Stata matrix with square of standard errors on diagonal, zero' off diagonal
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


// YS
// bs_draw_nowgt(YS)
// bs_draw_nowgt(YS,(10\1\1\1\0))


// check draw sample
YS = (1\2\3\4\5\6\7\8)
N00=N01=N10=N11=2
for (i=1; i<=1000; i++) {
	bs_draw = bs_draw_wgt((1\2),(3\4),(5\6),(7\8),N00, N01, N10, N11)
	if (i==1) bs_draw
	if (i==1) avg = bs_draw'
	else      avg = (avg \ bs_draw')
}
meanvariance(avg)'
colmin(avg)'
colmax(avg)'

popsize=rows(YS)
tempwgt = (1\9\1\9\1\9\1\9)
cumsum00 = quadrunningsum(tempwgt[1\2])
cumsum01 = quadrunningsum(tempwgt[3\4])
cumsum10 = quadrunningsum(tempwgt[5\6])
cumsum11 = quadrunningsum(tempwgt[7\8])
popsize00 = round(cumsum00[N00]) // the number of obs. in each group is rounded to the nearest integer
popsize01 = round(cumsum01[N01])
popsize10 = round(cumsum10[N10])
popsize11 = round(cumsum11[N11])
cumsum00 = cumsum00/cumsum00[N00] // normalize to sum to one within groups
cumsum01 = cumsum01/cumsum01[N01]
cumsum10 = cumsum10/cumsum10[N10]
cumsum11 = cumsum11/cumsum11[N11]

cumsum00
popsize00

for (i=1; i<=1000; i++) {
	bs_draw = bs_draw_wgt((1\2),(3\4),(5\6),(7\8),N00, N01, N10, N11,  cumsum00, cumsum01, cumsum10, cumsum11, popsize00, popsize01, popsize10, popsize11)
	if (i==1) bs_draw
	if (i==1) avg = bs_draw'
	else      avg = (avg \ bs_draw')
}
meanvariance(avg)'
colmin(avg)'
colmax(avg)'

popsize=1000
cumsum00 = quadrunningsum(tempwgt[1\2])
cumsum01 = quadrunningsum(tempwgt[3\4])
cumsum10 = quadrunningsum(tempwgt[5\6])
cumsum11 = quadrunningsum(tempwgt[7\8])
popsize00 = round(cumsum00[N00]/colsum(tempwgt)*popsize) // the number of obs. in each group is rounded to the nearest integer
popsize01 = round(cumsum00[N00]/colsum(tempwgt)*popsize)
popsize10 = round(cumsum00[N00]/colsum(tempwgt)*popsize)
popsize11 = round(cumsum00[N00]/colsum(tempwgt)*popsize)
cumsum00 = cumsum00/cumsum00[N00] // normalize to sum to one within groups
cumsum01 = cumsum01/cumsum01[N01]
cumsum10 = cumsum10/cumsum10[N10]
cumsum11 = cumsum11/cumsum11[N11]


cumsum00
popsize00

for (i=1; i<=1000; i++) {
	bs_draw = bs_draw_wgt((1\2),(3\4),(5\6),(7\8),N00, N01, N10, N11,  cumsum00, cumsum01, cumsum10, cumsum11, popsize00, popsize01, popsize10, popsize11)
	if (i==1) bs_draw
	if (i==1) avg = bs_draw'
	else      avg = (avg \ bs_draw')
}
meanvariance(avg)'
colmin(avg)'
colmax(avg)'



mata describe
mata memory


end
/* * * * *  END OF MATA BLOCK * * * * */



program define cicgraph, sortpreserve

	// parse arguements
	syntax [, ///
		Ci(name) /// {normal|percentile|bc|bca|*} ///
			/// ci() expecting normal, percentile, bc, or bca. However it will still
			/// work if there is a conforming matrix named e(ci_`ci')
			/// by default, -normal- is used if matrix e(ci_normal) exists
		Equations(namelist) /// {continuous, discrete_ci, dci_lower_bnd, and/or dci_upper_bnd}
		Name(name) /// graph names; if >1  name in equations(), names get .  graphs in memory will be replaced.
		*]

	if ("`e(cmd)'"!="cic") | (lower(e(vcetype))!="bootstrap") error 301

	// pull coefficients
	tempname b int
	tempvar coef eqn p pctile mean meanll meanul meanaxis
	matrix `b'  = e(b)
	matrix `b'  = `b''
	svmat `b' , names("`coef'")

	// cic option
	cap qui di colsof(e(ci_normal))
	if !_rc & mi("`ci'")         local ci "normal"
	if regexm("`ci'","^ci_(.*)") local ci = regexs(1) // if someone types ci_normal instead of normal, drop the "ci_"
	if !mi("`ci'") {
		if !inlist("`ci'","normal","percentile","bc","bca") di as txt "ci() expecting normal, percentile, bc, or bca. However it will work if there is a conforming matrix named e(ci_`ci')"
		matrix `int' =  e(ci_`ci')
		matrix `int' = `int''
		local level =  e(level)
		svmat `int',
	}
	else {
		gen `int'1=.
		gen `int'2=.
	}

	local rows = rowsof(`b')
	if (c(N) < `rows') set obs `rows'

	local names : rownames `b'
	local eqns  : roweq `b'

	local ylab : var lab `e(depvar)'
	if mi(`"`ylab'"') local ylab = e(depvar)

	qui {
		gen `eqn' = ""
		gen `p' = ""
		gen `pctile' = .
		forvalues i = 1/`rows' {
			local thiseqn  : word `i' of `eqns'
			local thisname : word `i' of `names'
			replace `eqn' = "`thiseqn'"  in `i'
			replace `p'   = "`thisname'" in `i'
			if regexm("`thisname'","^q([0-9]*)_?([0-9]*)") replace `pctile' = real(regexs(1)+"."+regexs(2)) in `i'
		}

		bys `eqn' (`p'): gen `mean'     = `coef'[1] if inlist(_n,1,_N)
		by  `eqn'      : gen `meanll'   = `int'1[1] if inlist(_n,1,_N)
		by  `eqn'      : gen `meanul'   = `int'2[1] if inlist(_n,1,_N)
		gen     `meanaxis' = 0
		replace `meanaxis' = 100 if `p'!="mean"
	}

	// graph results
	if mi("`equations'") local equations "continuous"
	local c=1
	foreach eqnname of local equations {
		if      "`eqnname'"=="continuous"    local eqnlabel "CIC model with continuous outcomes"
		else if "`eqnname'"=="discrete_ci"   local eqnlabel "CIC model with discrete outcomes"
		else if "`eqnname'"=="dci_lower_bnd" local eqnlabel "Discrete CIC model lower bound"
		else if "`eqnname'"=="dci_upper_bnd" local eqnlabel "Discrete CIC model upper bound"
		else {
			di as error "eqnname() should be continuous, discrete_ci, dci_lower_bnd, dci_upper_bnd"
			error 198
		}

		if !mi("`ci'") local addlegendlabels label(4 "Quantiles `level'% CI") label(2 "Mean `level'% CI")
		if !mi("`ci'") local addlegendorder  "4 2"

		graph twoway ///
			(scatter `mean'        `meanaxis', sort pstyle(p2) connect(L) msymbol(none) lwidth(*1.25)) ///
			(scatter `meanll'      `meanaxis', sort pstyle(p2) connect(L) msymbol(none) lwidth(*.85)  lpattern(dash)) ///
			(scatter `meanul'      `meanaxis', sort pstyle(p2) connect(L) msymbol(none) lwidth(*.85)  lpattern(dash)) ///
			(rcap    `int'1 `int'2 `pctile',   sort pstyle(p1) lcolor(*.85) lwidth(*.85)) ///
			(scatter `coef'        `pctile',   sort pstyle(p1) msize(*1.15) connect(L)) ///
				if `eqn'=="`eqnname'", ///
				legend(order(5 1 `addlegendorder') cols(2) label(5 "Quantiles") label(1 "Mean") `addlegendlabels' ) ///
				xtitle( "Quantile" ) ytitle(`"`ylab'"') title( "`eqnlabel'") subtitle( "`=e(footnote)'" ) ///
				name(`name'`=cond(`c'>1,"`c'","")',replace) `options'
		local ++c
	}
end // end of cic program definition

cd "C:\Users\keith\Desktop\CIC\"


* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* test cic_vce_parse
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
cic_vce_parse, vce(boot, reps(1000) saving(myfile.dta, replace) sepercent)
return list
cap nois cic_vce_parse, vce(none, reps(25))
return list
cic_vce_parse, vce(boot, reps(25) mse accel(myvector) saving("c:\temp\test.dta", replace))
return list


* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
set tracedepth 3
if 0  set trace on
else  set trace off
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
local Nreps = 15
if 01     	local vce vce(bootstrap, reps(`Nreps'))
else       	macro drop _vce
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
	label define high_after 1 "Control Group, 1st Period"   ///
	                        2 "Control Group, 2nd Period"   ///
	                        3 "Treatment Group, 1st Period" ///
	                        4 "Treatment Group, 2nd Period"
	label val high_after high_after
}
set seed 1

cap log close
log using cid_test_aid_data.log, replace

mac list _Nreps _vce

* Temp stuff
gen tempweight = 1
replace tempweight = 2 in 1
// bys high after: replace tempweight = 50 if _n<=20

egen agegroup = cut(age), group(7)


* Basic
timer on 21
cic  y high after, at(5(10)95) `vce' did
timer off 21

* With control variables
timer on 22
cap nois cic  y high after i.agegroup, did `vce' round(.25)
timer off 22

* Test recall
cic
cicgraph , name(r0)

* With weights
timer on 23
cap nois cic ly high after [fw=tempweight], did  `vce'
timer off 23

* With control variables and weights
timer on 24
cap nois cic ly high after i.agegroup [fw=tempweight], did `vce' round(.25)
timer off 24
timer list

exit

// Table 1
count
tabstat y ly , by(high_after) s(count mean sd min p25 p50 p75 p90 max) columns(s)  labelwidth(30) nototal format(%9.2f)

// DID estimate
reg y high##after
reg ly high##after

cic  ly high after ,  at(5(10)95) `vce'
ereturn list
*set trace on

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


// graphs
cic  y high after ,  at(1 5(2.5)90) `vce'
cicgraph,  name(g) e(continuous discrete_ci dci_lower_bnd dci_upper_bnd)

// compare vce() option above to the bootstrap prefix

set seed 1
timer on 10
cic y high after, vce(bootstrap, reps(`Nreps'))
timer off 10
timer list 10
ereturn list
estat bootstrap , all // estat bootstrap does work

cic // test replay works
est store a

timer on 12
cap nois bootstrap, reps(`Nreps') strata(high after) : cic y high after
timer off 12
timer list 12

// test if selected vars works
cap nois bootstrap [continuous]_b[mean], reps(20) strata(high after) : cic y high after

// test jacknife
jacknife: cic y high after if uniform()<.1

// check weights are working
gen testweight =(uniform()<.95) + (uniform()<.20)
tab testweight
cic y  high after [fw=testweight],  at(25 50 75 90)
drop if testweight==0
expand testweight
cic y  high after                ,  at(25 50 75 90)
cic y  high after                ,  at(25 50 75 90)


// direct comparision of vce() options
est restore a

set seed 1
cic y high after, vce(bootstrap, reps(`Nreps') sepercentile)



// we should get an error if try to use weights
cap nois {
  cic y  high after [fw=testweight],  at(25 50 75 90) vce(bootstrap, reps(25) nodots)
}

// we should get an error with vce(delta) since I haven't coded it up yet
cap nois cic  y high after ,  at(25 50 75 90) vce(delta)


ereturn list



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
local n_g=500
expand `n_g'
bys p t: gen y = _n / `n_g'
gen     d = 1.75 - 1.5 * y if t==0 & p==0
replace d = 0.75 - 0.5 * y if t==0 & p==1
replace d = 0.80 - 0.4 * y if t==1 & p==0
replace d = 0.50 - 1.0 * y if t==1 & p==1
replace d = round(d,.01)
cic d treat post,  at(10(10)90) vce(b)

exit

*set trace on
cicgraph , name(r1)
cicgraph , ci(ci_percentile) name(r2)

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

mata:
mata set matastrict on

// PARELLEL FUNCTION TO CIC() WITH MORE LEGIBLE CODE
// CIC estimates only -- No bootstrapping, SEs, or DID
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

	// Results will be returned into a structure w/ 4 vectors for con, dci, lower, upper
	struct cic_result scalar result

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
