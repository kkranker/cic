*! *! version $Id$
*! Plot Changes-in-Changes Model Estimates
*!
*! Following cic, this program creates a twoway plot of the model estimates.
*!
*! $Date$

program define cicgraph, sortpreserve rclass

	// parse arguements
	syntax [, ///
		Ci(name) /// Matrix name with intervals to use {normal|percentile|bc|bca|*} ///
			/// ci() expecting normal, percentile, bc, or bca. However cicgraph will still
			/// work if there is a conforming matrix named e(ci_`ci').
			/// cicgraph might give strange error messages if e(ci_`ci')
			/// ci(none) suppresses confidence intervals
			/// by default, -normal- is used if matrix e(ci_normal) exists. 
			///      if it e(ci_normal) not exist, options e(ci_percentile) if available
			/// This will work whether or not you include "ci_"  (e.g., "ci_normal" and "normal" are both valid)
		Equations(namelist) /// {continuous, discrete_ci, dci_lower_bnd, dci_upper_bnd, and/or qdid}
		                    /// default is all the equations available
		Name(string) /// graph names; default is name(cicgraph)
			/// if >1  name in equations(), name is treated as a prefix.
			/// default is name(cicgraph, replace)
		*]

	cap assert ("`e(cmd)'"=="cic")
	if _rc {
		di as error "cicgraph was written as a post-estimation command for cic."
		error 301
	}

	// pull coefficients
	tempname b int
	tempvar coef eqn p pctile mean meanll meanul meanaxis
	matrix `b'  = e(b)
	matrix `b'  = `b''
	svmat `b' , names("`coef'")

	// cic option
	if mi("`ci'") {
		cap qui di colsof(e(ci_normal))
		if !_rc local ci "normal"
		else {
			cap qui di colsof(e(ci_percentile))
			if !_rc local ci "ci_percentile"
		}
	}
	if regexm("`ci'","^ci_(.*)") local ci = regexs(1) // if someone types ci_normal instead of normal, drop the "ci_"
	cap qui di colsof(e(ci_normal))
	if !mi("`ci'") {
		if !inlist("`ci'","normal","percentile","bc","bca") di as txt "ci() expecting normal, percentile, bc, or bca. However it will work if there is a conforming matrix named e(ci_`ci')"
		matrix `int' =  e(ci_`ci')
		matrix `int' = `int''
		local level =  e(level)
		svmat `int',
	}
	else {
		qui gen `int'1=.
		qui gen `int'2=.
	}

	// size of e(b)
	local rows = rowsof(`b')
	if (c(N) < `rows') set obs `rows'
	if !mi("`ci'") {
		if (rowsof(`int')!=rowsof(`b') | colsof(`int')!=2) {
			di as error "e(ci_`ci') has the incorrect number of rows or columns."
			error 198
		}
	}

	// labels
	local names : rownames `b'
	local eqns  : roweq `b'
	local ylab : var lab `e(depvar)'
	if mi(`"`ylab'"') local ylab = e(depvar)

	// add temporary variables containing e(b)' equation names and quantiles
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
	if mi("`equations'") {
		local equations : list uniq eqns
		local equations : subinstr local equations "did_model" ""
		local equations : subinstr local equations "did"       ""
	}
	if mi("`name'") local name "cicgraph, replace"
	_parse comma name name_rhs:  name
	local c=1
	foreach eqnname of local equations {
		if      "`eqnname'"=="continuous"    local eqnlabel "Continuous CIC model"
		else if "`eqnname'"=="discrete_ci"   local eqnlabel "Discrete CIC model (under the conditional independence assumption)"
		else if "`eqnname'"=="dci_lower_bnd" local eqnlabel "Lower bound for the discrete CIC model (without conditional independence)"
		else if "`eqnname'"=="dci_upper_bnd" local eqnlabel "Upper bound for the discrete CIC model (without conditional independence)"
		else if "`eqnname'"=="qdid"          local eqnlabel "Quantile DID model"
		else {
			di as error "eqnname() should be continuous, discrete_ci, dci_lower_bnd, dci_upper_bnd, and/or qdid"
			error 198
		}
		if wordcount(`"`equations'"') == 1 local graphname `name'
		else                               local graphname `name'`c'
		local graphnamelist `graphnamelist' `graphname'

		if !mi("`ci'") local addlegendlabels label(4 "`level'% CI") label(2 "`level'% CI")
		if !mi("`ci'") local addlegendorder  "4 2"

		graph twoway ///
			(scatter `mean'        `meanaxis', sort pstyle(p2) connect(L) msymbol(none) lwidth(*1.25)) ///
			(scatter `meanll'      `meanaxis', sort pstyle(p2) connect(L) msymbol(none) lwidth(*.85)  lpattern(dash)) ///
			(scatter `meanul'      `meanaxis', sort pstyle(p2) connect(L) msymbol(none) lwidth(*.85)  lpattern(dash)) ///
			(rcap    `int'1 `int'2 `pctile',   sort pstyle(p1) lcolor(*.85) lwidth(*.85)) ///
			(scatter `coef'        `pctile',   sort pstyle(p1) msize(*1.15) connect(L)) ///
				if `eqn'=="`eqnname'", ///
				legend(order(5 1 `addlegendorder') cols(2) label(5 "CIC at quantiles") label(1 "Mean CIC") `addlegendlabels' ) ///
				xtitle( "Quantile" ) ytitle(`"`ylab'"') note( "`eqnlabel'" "`=e(footnote)'" ) ///
				name(`graphname' `name_rhs') `options'
		local ++c
	}
	return local cmd       cicgraph
	return local name      `graphnamelist'
	return local equations `equations'
	if mi("`ci'") return local ci none
	else          return local ci ci_`ci'
	
end // end of cicgraph program definition
