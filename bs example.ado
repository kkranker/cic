* Imagine we would like to bootstrap the standard errors of a command using a bootstrap routine.

  * I created a previous post demonstrating how to write a bootstrap command.  This is a similar post however the bootstrap is much faster than the previous one.
 
  * First let's generate some data.
 
  clear
  set obs 1000
  gen x1 = rnormal()
  gen x2 = 2*rnormal()
  gen u = 5*rnormal()
 
  gen y = 5 + x1 + x2 + u
 
  local dependent_var x1 x2
 
  local command_bs reg y `dependent_var'
 
  * First let's see how the base command works directly.
  `command_bs'
 
  * As a matter of comparison this is the built in bootstrap command.
  bs, rep(100): `command_bs'
 
  * The following code is yet another user written bootstrap alternative.
 
  * I wrote this
 
  * Specify the number of bootstrap iterations
  local bs = 100
 
  * Save the number of observations to be drawn
  local N_obs = _N
 
  * Number of terms plus one for the constant of coefficients to be saved
  local ncol = wordcount("`dependent_var'")+1
 
  mata theta = J(`bs',`ncol',0)
  forv i = 1(1)`bs' {
    * Preserve the initial form of the data
preserve

    * Draw the indicators of the resample
mata draw = ceil(runiform(`N_obs',1):*`N_obs')

* Create an empty matrix to hold the number of items to expand
mata: expander = J(`N_obs',1, 0)

* Count the number of items per observation to generate.
    qui mata: for (i=1 ; i <= `N_obs'; i++) expander[i]=sum(draw:==i) ; "draws complete"

* Pull the expander value into stata
getmata expander=expander

* Drop unnessessary data
qui drop if expander == 0

* Expand data, expand==1 does nothing
    qui expand expander

* Run a regression
    qui `command_bs'

* Send to mata the matrix of results
mata theta[`i',] = st_matrix("e(b)")
   
* Configure the visual display
di _c "."
    if (int(`i'/10)==`i'/10) di _c "|"
    if (int(`i'/50)==`i'/50) di " " `i'

restore
  }
 
  * The estimates of the coefficients have been saved into a theta matrix.
  mata theta
 
  * Now let's calculate the standard deviations.
  mata bs_var = variance(theta)
  mata bs_var
  mata bs_se = diag(bs_var):^.5
  mata bs_se