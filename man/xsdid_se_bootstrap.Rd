% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap.R
\name{xsdid_se_bootstrap}
\alias{xsdid_se_bootstrap}
\title{Compute bootstrap standard errors for SDID with time varying covariates}
\usage{
xsdid_se_bootstrap(
  panel,
  unit = 1,
  time = 2,
  outcome = 3,
  treatment = 4,
  x,
  B = 100,
  num.cores = 1
)
}
\arguments{
\item{panel}{A data frame with columns consisting of units, time, outcome, and treatment indicator. It should be a balanced panel and not contain any NA.}

\item{unit}{The column number/name corresponding to the unit identifier. Default is 1.}

\item{time}{The column number/name corresponding to the time identifier. Default is 2.}

\item{outcome}{The column number/name corresponding to the outcome identifier. Default is 3.}

\item{treatment}{The column number/name corresponding to the treatment status. Default is 4.}

\item{x}{The column numbers/names of all additional control variables}

\item{B}{number of bootstrap samples (Default 100).}

\item{num.cores}{if larger than 1 use the parallel package to perform the bootstrap estimation on multiple cores. (More than 1 core does not work on Windows.)}
}
\value{
As list with 2 elements:
  - se is the estimated standard error of the treatment effect tau.
  - boot.tau is a vector containing the estimated treatment effects tau for each bootstrap replication.
}
\description{
Each clustered bootstrap sample is a sample of units with replacement.
It is guranteed that at least 2 treatment group and at least 2 control group units are selected.
(Ideally at least 1 of each should suffice but just a single observation in
 a group can lead to downstream errors when calling synthdid_estimate. This issue needs to be
 explored further)
}
\details{
In each bootstrap sample, we use the function adjust.outcome.for.x to adjust for the covariates.
}
\examples{
\dontrun{
 dat = xsdid.mc(N=20, T=20,return.data = TRUE)
 xsdid_estimate(dat, unit="i",time = "t",outcome = "y",treatment = "treat_exp",x = "x")
 xsdid_se_bootstrap(dat, unit="i",time = "t",outcome = "y",treatment = "treat_exp",x = "x", B=100)$se
}
}
