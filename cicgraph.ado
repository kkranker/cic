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
			/// by default, -normal- is used if matrix e(ci_normal) exists
			/// This will work whether or not you include "ci_"  (e.g., "ci_normal" and "normal" are both valid)
		Equations(namelist) /// {continuous, discrete_ci, dci_lower_bnd, and/or dci_upper_bnd}
		Name(name) /// graph names; default is name(cicgraph)
			/// if >1  name in equations(), name is treated as a prefix.
			/// graphs in memory will be replaced.
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

	// size of e(b)
	local rows = rowsof(`b')
	if (c(N) < `rows') set obs `rows'
	if (!mi("`ci'") & (rowsof(`int')!=rowsof(`b') | colsof(`int')!=2)) {
		di as error "e(ci_`ci') has the incorrect number of rows or columns."
		error 198
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
	if mi("`equations'") local equations continuous discrete_ci dci_lower_bnd dci_upper_bnd
	if mi("`name'")      local name "cicgraph"
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
				name(`name'`=cond(`: word count `equations''>1,"`c'","")',replace) `options'
		local ++c
		local graphnamelist `graphnamelist' `name'`=cond(`: word count `equations''>1,"`=`c'","")'
	}
	return local cmd       cicgraph
	return local name      `graphnamelist'
	return local equations `equations'
	return local ci        ci_`ci'
end // end of cicgraph program definition
