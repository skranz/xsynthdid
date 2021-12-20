---
title: "xsynthdid"
author: Sebastian Kranz, Ulm University
output: 
  html_document: 
    keep_md: yes
---



### Overview

[Arkhangelsky et al. (2021)](https://www.aeaweb.org/articles?id=10.1257/aer.20190159) introduce the novel synthetic difference-in-differences (SDID) method that is implemented in the R package [synthdid](https://github.com/synth-inference/synthdid).

This package has the function `adjust.outcome.for.x` that can be used to adjust the outcome for exogenous time-varying covariates. For more background see my vignette [Synthetic Difference-in-Differences with Time-Varying Covariates (2021)]().

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
sdid.adj = synthdid_estimate(Y=pm$Y,N0=pm$N0,T0 = pm$T0) # inconsistent
sdid.adj
```

```
## synthdid: 49.398 +- 0.252. Effective N0/N0 = 15.0/15~1.0. Effective T0/T0 = 6.9/15~0.5. N1,T1 = 15,15.
```

Our estimate is now pretty close to the true causal effect.
