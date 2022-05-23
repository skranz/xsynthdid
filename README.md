Author: Sebastian Kranz, Ulm University


### Overview

[Arkhangelsky et al. (2021)](https://www.aeaweb.org/articles?id=10.1257/aer.20190159) introduce the novel synthetic difference-in-differences (SDID) method that is implemented in the R package [synthdid](https://github.com/synth-inference/synthdid).

This package has the function `adjust.outcome.for.x` that can be used to adjust the outcome for exogenous time-varying covariates. For more background see my vignette [Synthetic Difference-in-Differences with Time-Varying Covariates (2021)](https://github.com/skranz/xsynthdid/blob/main/paper/synthdid_with_covariates.pdf).

If you use the package for research, you could cite the vignette as:

Kranz, Sebastian (2021), "Synthetic Difference-in-Differences with Time-Varying Covariates", mimeo 

### Installation

You can install the package from my [r-universe repository](https://skranz.r-universe.dev/ui#builds) by using the following code:


```r
options(repos = c(skranz = 'https://skranz.r-universe.dev',
    CRAN = 'https://cloud.r-project.org'))
install.packages("xsynthdid")
```

### Usage example

We first simulate a panel data set with a true treatment effect of `tau=50` using the function `xsdid.mc` in our package.

```r
set.seed(1)
library(xsynthdid)
dat = xsdid.mc(N=30, T=30,tau=50,return.data = TRUE)
head(dat)
```

```
##   i t   group treat exp treat_exp time.trend        eps     x.load   xc.load
## 1 1 1 control     0   0         0  0.1439891  0.1828379 0.07844164 0.1473573
## 2 1 2 control     0   0         0 -0.2414066  0.1526865 0.07844164 0.1473573
## 3 1 3 control     0   0         0  0.1752374  0.2351750 0.07844164 0.1473573
## 4 1 4 control     0   0         0  1.8588536 -0.7058448 0.07844164 0.1473573
## 5 1 5 control     0   0         0  3.3896764  0.2107316 0.07844164 0.1473573
## 6 1 6 control     0   0         0  2.9570783 -0.2301411 0.07844164 0.1473573
##             x         xc          y
## 1  0.12850395 -0.8351821  13.807491
## 2  0.06349820 -1.2433237  -8.742734
## 3  0.01129829 -1.2056680   9.561960
## 4 -0.04246972 -1.0943951  90.113348
## 5 -0.20068988 -1.2214046 159.660058
## 6 -0.15537278 -1.0809002 139.855136
```

We first perform SDID estimation without any adjustment for the exogenous covariate `x`:

```r
pm = panel.matrices(dat, unit="i",time = "t",outcome = "y",treatment = "treat_exp")
sdid = synthdid_estimate(Y=pm$Y,N0=pm$N0,T0 = pm$T0) # inconsistent
sdid
```

```
## synthdid: 39.760 +- 3.015. Effective N0/N0 = 14.5/15~1.0. Effective T0/T0 = 1.0/15~0.1. N1,T1 = 15,15.
```

We see how the estimator is a bit off the true causal effect.

We now adjust the outcomes `y` using the `adjust.outcome.for.x` function and use the adjusted outcome for the SDID estimation.

```r
dat$y.adj = adjust.outcome.for.x(dat,unit="i",time = "t",
  outcome = "y",treatment = "treat_exp", x="x")
pm = panel.matrices(dat, unit="i",time = "t",outcome = "y.adj",treatment = "treat_exp")
sdid.adj = synthdid_estimate(Y=pm$Y,N0=pm$N0,T0 = pm$T0)
sdid.adj
```

```
## synthdid: 49.398 +- 0.252. Effective N0/N0 = 15.0/15~1.0. Effective T0/T0 = 6.9/15~0.5. N1,T1 = 15,15.
```

Our estimate is now pretty close to the true causal effect.

### Standard errors

Note that the shown standard errors do not account for the initial adjustment for time varying covariates. Experimentally, I have included the function `xsdid_se_bootstrap` that allows to compute the standard errors using a clustered bootstrap approach similar to Algorithm 2 in [Arkhangelsky et al. (2021)](https://www.aeaweb.org/articles?id=10.1257/aer.20190159).

To speed up the building of the README, here is an example with just `B=10` bootstrap replications. Of course, you should set `B` to a larger number like 500.


```r
xsdid_se_bootstrap(dat,B=10, unit="i",time = "t",outcome = "y",treatment = "treat_exp",x = "x")
```

```
## $se
## [1] 0.1832767
## 
## $boot.tau
##  [1] 49.49238 49.59243 49.33778 49.64696 49.22328 49.31334 49.31244 49.35124
##  [9] 49.56474 49.88914
```




