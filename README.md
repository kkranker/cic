# cic: Stata implementation of the Athey and Imbens (2006) Changes-in-Changes model


# Description

This Stata command, `cic`, implements the changes-in-changes (CIC) model proposed by Athey and Imbens (2006).
The command estimates the average and quantile treatment effects of a treatment in settings where
repeated cross sections of individuals are observed in a treatment group and a control group, before and after the treatment.
The CIC model relaxes several assumptions of the standard linear difference-in-differences model.
Both the continuous CIC model and discrete CIC model (with and without a conditional independence assumption)
are included in the cic command, as are treatment effects from
standard linear difference-in-differences model and a quantile difference-in-differences model.
By implementing the CIC estimator alongside these other two pre-existing estimators,
the `cic` command can illustrate how the effect the treatment varies across a variety of assumptions.


# Background

I wrote the code for this Stata implementation for one of my dissertation projects (Kranker 2011, 2016).
The code started as a simple "port" of [Athey and Imbens' Matlab code](https://athey.people.stanford.edu/research).
Then I made some changes to the code to speed it up in various ways and extend the methods.

While was rummaging through some old files in June 2019,
I ran across this code, thought it might be helpful to others, and decided to post the package online.
I see Blaise Melley also has [his own version of the model for Stata](https://sites.google.com/site/blaisemelly/home/computer-programs/cic_stata).
Let's hope this release prevents yet another person from doing this work!


# Author

Keith Kranker

Based on Matlab code by Susan Athey & Guido W. Imbens,
available at https://athey.people.stanford.edu/research



# Installation

To install from Github, type this from your Stata command line:

```stata
. net from https://raw.githubusercontent.com/kkranker/cic/master/
```

(coming soon) To install from SSC, type this from your Stata command line:

```stata
. net describe cic
```

The Stata help file (cic.sthlp) provides additional documentation and examples.


# Remarks

The use of difference-in-differences (DID) methods is widespread in program evaluation and empirical
economics (Imbens and Wooldridge 2009). DID methods involve comparing outcomes between two groups across two time periods,
where only one of the two groups are exposed to the intervention in one of the periods.
The DID estimator calculates the difference in outcomes between the treatment and comparison groups after the intervention began,
minus the difference in outcomes between the treatment and comparison groups before the intervention began.
Or, equivalently, the DID estimator can be seen as the change in outcomes for the treatment group
before and after the intervention, minus the change in outcomes for the comparison group over the same time period.
It is straightforward to generalize this basic two group, two period DID model has been generalized in various ways—for example,
to adjust for observed covariates, include more than two groups, or include more than two time periods.

Athey and Imbens (2006) proposed a changes-in-changes (CIC) model which generalizes the DID model by relaxing relaxes several assumptions.
(Thus the standard DID model is nested in the CIC model as a special case.)
The CIC model estimates the entire distribution of outcomes under the counterfactual,
allowing one to calculate average treatment effects or estimate effects at specific quantiles.
This Stata command, `cic`,  implements the CIC estimator from Athey and Imbens (2006).
`cic` is written in Mata with an effort to maximize parallel computing; in tests (not shown), I found `cic` estimated
the model more quickly that the Matlab code previously distributed by the Athey and Imbens.
The cic command also offers several previously unavailable features (e.g., to allow for covariates).
In addition, you can use Stata's `bootstrap:` prefix, which offers more flexibility for
computing bootstrapped standard errors (e.g., strata, blocks).

Imbens and Wooldridge (2009) provide a nice, short overview of the method.
Athey and Imbens (2006) explain their methods (with proofs) in a fairly long and complicated article.
The appendix is also quite helpful.

The Stata help file (cic.sthlp) provides additional documentation and examples.


# Referemces

* Athey, Susan and Guido W. Imbens. "Identification and Inference in Nonlinear Difference-in-Differences Models." *Econometrica*, vol. 74, no. 2, March 2006, pp. 431-497. (http://dx.doi.org/10.1111/j.1468-0262.2006.00668.x)

* Imbens, Guido W. and Jeffery M. Wooldridge. “Recent Developments in the Econometrics of Program Evaluation.” *Journal of Economic Literature*, vol. 47, no. 1, 2009, pp. 5–86. (http://dx.doi.org/10.1257/jel.47.1.5)

* Kranker, Keith. “The Effect of Disease Management Programs on Medicaid Expenditures.” Doctoral dissertation. College Park, MD: University of Maryland, 2011. (http://hdl.handle.net/1903/12101)

* Kranker, Keith. “Effects of Medicaid Disease Management Programs on Medical Expenditures: Evidence from a Natural Experiment in Georgia.” *Journal of Health Economics*, vol. 46, March 2016, pp. 52-69. (http://dx.doi.org/10.1016/j.jhealeco.2016.01.008)
