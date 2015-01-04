*! $Id$
*! Examples of Stata implimetation of the changes-in-changes model (cic.ado)
*! Stata code by Keith Kranker
*! Last updated $Date$

clear all
cap    nois cd "C:\Users\keith\Desktop\cic"
if _rc nois cd "C:\Users\kkranker\Documents\Dissertation\Stata-Changes-in-Changes"
include cic.ado



* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
set tracedepth 3
if 0  set trace on
else  set trace off
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
local Nreps = 200
if 01     	local vce "vce(bootstrap, reps(\`Nreps'))"
else       	macro drop _vce
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *



* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* Benchmark againt results in Athey and Imbens supplimental appendix
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

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
log using cic_benchmark_testing.log, replace

mac list _Nreps _vce

// Table 1
count
tabstat y ly , by(high_after) s(count mean sd min p25 p50 p75 p90 max) columns(s)  labelwidth(30) nototal format(%9.2f)

// DID estimate
reg y high##after
reg ly high##after

// cic estimates from A&I Appendix
// Table 2
timer on 2
cic all  y high after ,  at(25 50 75 90) `vce' did
cic all ly high after ,  at(50)          `vce' did
timer off 2
timer list 2



// Table 3
timer on 3
cic all  y high after ,  at(25 50 75 90) `vce' untreated did
cic all ly high after ,  at(50)          `vce' untreated did
timer off 3
timer list 3

// graphs
cic all  y high after ,  at(1 5(2.5)90) `vce'
ereturn list
cicgraph,  name(g) e(continuous discrete_ci dci_lower_bnd dci_upper_bnd)


// VCE via delta metho
cic all  y high after, vce(delta)     at(25 50 75 90)
cic all  y high after, vce(delta) did at(25 50 75 90)


* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* Test misc features with this dataset
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

* Basic
cic all  y high after, at(5(10)95) `vce' did

* With control variables
egen agegroup = cut(age), group(7)
cic all  y high after i.agegroup, did `vce' round(.25)

* Test recall
cic
cicgraph , name(r0)

* With weights
gen tempweight = 1
replace tempweight = 2 in 1
cic all ly high after [fw=tempweight], did  `vce'

* With control variables and weights
timer on 24
cic all ly high after i.agegroup [fw=tempweight], did `vce' round(.25)
timer off 24
timer list

// compare vce() option above to the bootstrap prefix
set seed 1
timer on 10
cic all y high after, vce(bootstrap, reps(`Nreps'))
timer off 10
timer list 10
ereturn list
estat bootstrap , all // estat bootstrap does work

cic // test replay works
est store a

timer on 12
cap nois bootstrap, reps(`Nreps') strata(high after) : cic all y high after
timer off 12
timer list 12

// test if selected vars works
cap nois bootstrap [continuous]_b[mean], reps(20) strata(high after) : cic all y high after

// test jacknife
jacknife: cic all y high after if uniform()<.1

// check fweights are working
// I should get the same results if I use fweights or if I expand the datset
preserve
  gen testweight =(uniform()<.95) + (uniform()<.20)
  tab testweight
  cic all y  high after [fw=testweight],  at(25 50 75 90)  `vce'
  drop if testweight==0
  expand testweight
  cic all y  high after                ,  at(25 50 75 90)  `vce'
restore

// check other weights are working
preserve
  gen  testweight = max(0,rnormal(1.25,.05))
  summ testweight
  cic all y  high after [iw=testweight],  at(25 50 75 90)  `vce'
  cic all y  high after [aw=testweight],  at(25 50 75 90)  `vce'
  cap nois cic all y  high after [pw=testweight],  at(25 50 75 90)  `vce'
  cap nois cic all y  high after [pw=testweight],  at(25 50 75 90)  vce(none)
restore

// check it works with svy bootstrap.
// This example is just quick-and-dirty.
// It does not account for the fact that you might want to include strata(high after) in some applications.
preserve
  forvalues i = 1/100 {
	gen  testweight`i' = max(0,rnormal(1.25,.05))
	local testweightlist `testweightlist' testweight`i'
  }
  bys  high after: gen rownum=_n
  xtile gid=rownum, nq(10)
  svyset gid, bsrweight(`testweightlist') vce(bootstrap)
  svy : cic all y high after
restore


// direct comparision of vce() options
est restore a
cic
set seed 1
cic all y high after, vce(boot, reps(`Nreps') sepercentile)
cic all y high after, vce(delta)
cic all y high after, vce(none)

// test various combinations of estimators/vce/graphs with graphs
cic all    y high after ,  at(25 50 75) vce(bootstrap, reps(5))
cicgraph, name(all)
cic dci    y high after ,  at(25 50 75) vce(bootstrap, reps(5))
cicgraph, name(dci)
cic continuous y high after ,  at(25 50 75) vce(bootstrap, reps(5))
cic bounds     y high after ,  at(25 50 75) vce(bootstrap, reps(5)) did
cicgraph, name(bounds) eq(qdid dci_lower_bnd)
cic continuous y high after ,  at(25 50 75) vce(none)
cicgraph,name(novce) 



* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* The following code can be used to test the program using another
* "fake" data set
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
sysuse nlsw88, clear
set seed 1
gen TREAT1 = uniform() < .5
replace wage = wage + TREAT1
gen POST1 = uniform() < .5
replace wage = wage - POST1

// bootstrap the sample conditional on Ngt for g; t = 0; 1
cic all wage TREAT1 POST1 i.occupation,  at(10(10)90 99.5)



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
cic all d treat post,  at(10(10)90) vce(b)

cicgraph , name(r1)
cicgraph , ci(ci_percentile) name(r2)




* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* The following code tests various sub-functions
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

/* * * * *  BEGIN MATA BLOCK * * * * */
version 11.2
mata:
mata set matastrict on
mata set matafavor speed
mata set matalnum on /* drop this later */

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

for (i=1; i<=100; i++) {
	bs_draw = bs_draw_wgt((1\2),(3\4),(5\6),(7\8),N00, N01, N10, N11,  cumsum00, cumsum01, cumsum10, cumsum11, popsize00, popsize01, popsize10, popsize11)
	if (i==1) bs_draw
	if (i==1) avg = bs_draw'
	else      avg = (avg \ bs_draw')
}
meanvariance(avg)'
colmin(avg)'
colmax(avg)'

cumdfinv((1::10), .94          )
cumdfinv((1::10), .94,          J(10,1,1))
cumdfinv((1::10), .9           )
cumdfinv((1::10), .9 ,          J(10,1,1))

cumdfinv((1::9) , .94          ,(2\J(8,1,1)))
cumdfinv((1::9) , .84          ,(J(8,1,1)\2))
cumdfinv((1::9) , .9           ,(2\J(8,1,1)))
cumdfinv((1::9) , .8           ,(J(8,1,1)\2))

// test fden w/ this vector
Y= ( 1.11 \ 1.1 \ 1.2 \ 1.3 \ 5 \ 5.1 \ 5.2 \ 10 \ 10 \ 10 )

// the following three should all equal 0.072278
fden(1.1,Y)                           // no weights
fden(1.1,Y, J(10,1,1) )               // weight = 1
fden(1.1,Y[1..8], (J(7,1,1) \ 3) )    // weight = 1, except handle three "10" at bottom

// the following two should equal  0.076557
fden(1.1,(Y\Y))                       // stack Y on top of itself
fden(1.1,Y, J(10,1,2))                // weight = 2



mata describe

end
/* * * * *  END OF MATA BLOCK * * * * */


* test cic_vce_parse
cic_vce_parse, vce(boot, reps(1000) saving(myfile.dta, replace) sepercent)
return list
cap nois cic_vce_parse, vce(none, reps(25))
return list
cic_vce_parse, vce(boot, reps(25) mse accel(myvector) saving("c:\temp\test.dta", replace))
return list

log close

