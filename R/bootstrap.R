example.bootstrap = function() {
  dat = xsdid.mc(N=20, T=20,return.data = TRUE)
  xsdid_estimate(dat, unit="i",time = "t",outcome = "y",treatment = "treat_exp",x = "x")
  xsdid_se_bootstrap(dat, unit="i",time = "t",outcome = "y",treatment = "treat_exp",x = "x", B=100)$se
}

#' A shortcut function to directly perform SDID estimation with covariates from a panel data set
#'
#' Calls adjust.outcome.for.y, panel.matrices and then synthdid_estimate. Has similar syntax to
#' xsdid_se_bootstrap.
#'
#' @param panel A data frame with columns consisting of units, time, outcome, and treatment indicator. It should be a balanced panel and not contain any NA.
#' @param unit The column number/name corresponding to the unit identifier. Default is 1.
#' @param time The column number/name corresponding to the time identifier. Default is 2.
#' @param outcome The column number/name corresponding to the outcome identifier. Default is 3.
#' @param treatment The column number/name corresponding to the treatment status. Default is 4.
#' @param x The column numbers/names of all additional control variables
#' @param x.rows To estimate the effect of x on y, we use by default all rows in which no treatment takes places. The argument x.rows allows to specify these rows manually. E.g. you could only take all rows before the first treatment starts (i.e. compared to the default exclude the rows of the control group during the treatment period).
#' @return The resulting object from the call to synthdid_estimate.
#' @examples
#'\dontrun{
#'  dat = xsdid.mc(N=20, T=20,return.data = TRUE)
#'  xsdid_estimate(dat, unit="i",time = "t",outcome = "y",treatment = "treat_exp",x = "x")
#'  xsdid_se_bootstrap(dat, unit="i",time = "t",outcome = "y",treatment = "treat_exp",x = "x", B=100)$se
#'}
#' @export
xsdid_estimate = function (panel, unit = 1, time = 2, outcome = 3, treatment = 4,x,x.rows=NULL) {
  y.adj = "...y.adj"
  dat[[y.adj]] = adjust.outcome.for.x(dat, unit=unit,time =time,outcome = outcome,treatment = treatment,x = x, xrows=xrows)
  pm = panel.matrices(dat, unit=unit,time =time,outcome = y.adj,treatment = treatment)
  sdid = synthdid_estimate(Y=pm$Y,N0=pm$N0,T0 = pm$T0)
  sdid
}

#' Compute bootstrap standard errors for SDID with time varying covariates
#'
#' Each clustered bootstrap sample is a sample of units with replacement.
#' It is guranteed that at least 2 treatment group and at least 2 control group units are selected.
#' (Ideally at least 1 of each should suffice but just a single observation in
#'  a group can lead to downstream errors when calling synthdid_estimate. This issue needs to be
#'  explored further)
#'
#' In each bootstrap sample, we use the function adjust.outcome.for.x to adjust for the covariates.
#'
#' @param panel A data frame with columns consisting of units, time, outcome, and treatment indicator. It should be a balanced panel and not contain any NA.
#' @param unit The column number/name corresponding to the unit identifier. Default is 1.
#' @param time The column number/name corresponding to the time identifier. Default is 2.
#' @param outcome The column number/name corresponding to the outcome identifier. Default is 3.
#' @param treatment The column number/name corresponding to the treatment status. Default is 4.
#' @param x The column numbers/names of all additional control variables
#' @param B number of bootstrap samples (Default 100).
#' @param num.cores if larger than 1 use the parallel package to perform the bootstrap estimation on multiple cores. (More than 1 core does not work on Windows.)
#' @return As list with 2 elements:
#'   - se is the estimated standard error of the treatment effect tau.
#'   - boot.tau is a vector containing the estimated treatment effects tau for each bootstrap replication.
#' @examples
#'\dontrun{
#'  dat = xsdid.mc(N=20, T=20,return.data = TRUE)
#'  xsdid_estimate(dat, unit="i",time = "t",outcome = "y",treatment = "treat_exp",x = "x")
#'  xsdid_se_bootstrap(dat, unit="i",time = "t",outcome = "y",treatment = "treat_exp",x = "x", B=100)$se
#'}
#' @export
xsdid_se_bootstrap = function (panel, unit = 1, time = 2, outcome = 3, treatment = 4,x, B=100, num.cores = 1) {
  restorepoint::restore.point("xsdid_se_bootstrap")
  units = unique(panel[[unit]])
  panel$..unit = panel[[unit]]
  panel$..treatment = as.integer(panel[[treatment]])
  tr.dat = panel %>%
    group_by(..unit) %>%
    summarize(treated = any(..treatment == 1))

  units = tr.dat$..unit
  co.units = tr.dat$..unit[!tr.dat$treated]
  tr.units = tr.dat$..unit[tr.dat$treated]

  if (length(tr.units)==1) {
    stop("You have only a single unit in the treatment group. But the bootstrap approach is only valid if you have more than one treated unit.")
  }
#
#   if (is.null(sdid)) {
#     # Repeat basic estimation
#     dat$..y.adj = adjust.outcome.for.x(dat, unit=unit,time =time,outcome = outcome,treatment = treatment,x = x)
#     pm = panel.matrices(dat, unit=unit,time =time,outcome = "..y.adj",treatment = treatment)
#     sdid = synthdid_estimate(Y=pm$Y,N0=pm$N0,T0 = pm$T0) # inconsistent
#   }
#   str(sdid)
#   # Get weights from baseline estimation as initial values for
#   # the bootstrap
#   weights = attr(sdid, "weights")
#   lambda = weights$lambda
#   omega = weights$omega


  one.bootstrap = function(...) {

    bdat = draw_bootstrap_sample(panel, co.units, tr.units,units = units)
    bdat$..y.adj = adjust.outcome.for.x(bdat, unit="..bunit",time =time,outcome = outcome,treatment = treatment,x = x)
    pm = panel.matrices(bdat, unit="..bunit",time =time,outcome = "..y.adj",treatment = treatment)
    restorepoint::restore.point("jslkfjkdfjlkdjf")
    bsdid = synthdid_estimate(Y=pm$Y,N0=pm$N0,T0 = pm$T0)
    as.numeric(bsdid)
  }

  if (num.cores > 1) {
    library(parallel)
    tau = unlist(mclapply(1:B, one.bootstrap, mc.cores=num.cores))
  } else {
    tau = unlist(lapply(1:B, one.bootstrap))
  }
  list(se = sd(tau)*(B-1)/B, boot.tau = tau)
}

# draw a valid bootstrap sample
draw_bootstrap_sample = function(dat, co.units, tr.units, units = c(co.units, tr.units)) {
  N = length(units)
  # in the moment there seems a bug later in synthdid_estimate
  # sometimes if we only have a single unit in a group
  # thus we pick at least 2

  if (length(tr.units)<2 | length(co.units)<2) {
    stop("In the moment we need at least 2 treatment group and 2 control group units.")
  }

  tr1 = sample(tr.units,2)
  co1 = sample(co.units,2)



  i.boot = c(tr1,co1,sample(units,N-4,replace=TRUE))
  bdat = left_join(data.frame(..unit=i.boot, ..bunit = 1:NROW(i.boot)), dat, by="..unit") %>%
    mutate(..treatgroup = 1L*(..unit %in% tr1)) %>%
    arrange(..treatgroup, ..bunit)
  bdat
}

