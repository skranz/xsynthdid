% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjust_y.R
\name{adjust.outcome.for.x}
\alias{adjust.outcome.for.x}
\title{Adjust outcome for using synthdid with additional covariates x}
\usage{
adjust.outcome.for.x(
  panel,
  unit = 1,
  time = 2,
  outcome = 3,
  treatment = 4,
  x,
  x.rows = NULL,
  add.mean = FALSE
)
}
\arguments{
\item{panel}{A data frame with columns consisting of units, time, outcome, and treatment indicator. It should be a balanced panel and not contain any NA.}

\item{unit}{The column number/name corresponding to the unit identifier. Default is 1.}

\item{time}{The column number/name corresponding to the time identifier. Default is 2.}

\item{outcome}{The column number/name corresponding to the outcome identifier. Default is 3.}

\item{treatment}{The column number/name corresponding to the treatment status. Default is 4.}

\item{x}{The column numbers/names of all additional control variables}

\item{x.rows}{To estimate the effect of x on y, we use by default all rows in which no treatment takes places. The argument x.rows allows to specify these rows manually. E.g. you could only take all rows before the first treatment starts (i.e. compared to the default exclude the rows of the control group during the treatment period).}

\item{add.mean}{Shall the average of the effect of x on the outcome over all observatios be added to the outcome? If TRUE, the mean of the adjusted outcome is the same as the mean of the original outcome. Default is FALSE.}
}
\value{
The computed adjusted values for y. They can be used as outcome column in a subsequent call to panel.matrices in order to then estimate synthdid adjusted for the control variables that where specified here in x.
}
\description{
Implements a way to adjust the outcome varible for additional covariates x.
The adjusted outcome can be used for a synthdid estimate that accounts for the covariates.
More precisely, you can use the name of the new adjusted outcome (by default "y.adj") as argument outcome in your call to panel.matrices that prepares the data for synthdid_estimate.
See the Vignete for an example.
}
