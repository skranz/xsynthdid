#' Adjust outcome for using synthdid with additional covariates x
#'
#' Implements a way to adjust the outcome varible for additional covariates x.
#' The adjusted outcome can be used for a synthdid estimate that accounts for the covariates.
#' More precisely, you can use the name of the new adjusted outcome (by default "y.adj") as argument outcome in your call to panel.matrices that prepares the data for synthdid_estimate.
#' See the Vignete for an example.
#'
#' @param panel A data frame with columns consisting of units, time, outcome, and treatment indicator. It should be a balanced panel and not contain any NA.
#' @param unit The column number/name corresponding to the unit identifier. Default is 1.
#' @param time The column number/name corresponding to the time identifier. Default is 2.
#' @param outcome The column number/name corresponding to the outcome identifier. Default is 3.
#' @param treatment The column number/name corresponding to the treatment status. Default is 4.
#' @param x The column numbers/names of all additional control variables
#' @param x.rows To estimate the effect of x on y, we use by default all rows in which no treatment takes places. The argument x.rows allows to specify these rows manually. E.g. you could only take all rows before the first treatment starts (i.e. compared to the default exclude the rows of the control group during the treatment period).
#' @param add.mean Shall the average of the effect of x on the outcome over all observatios be added to the outcome? If TRUE, the mean of the adjusted outcome is the same as the mean of the original outcome. Default is FALSE.
#' @return The computed adjusted values for y. They can be used as outcome column in a subsequent call to panel.matrices in order to then estimate synthdid adjusted for the control variables that where specified here in x.
#'
#'
#' @export

adjust.outcome.for.x = function (panel, unit = 1, time = 2, outcome = 3, treatment = 4,x, x.rows=NULL, add.mean=FALSE) {
  #restorepoint::restore.point("adjust.outcome.for.x")
  cols = colnames(panel)
  names(cols) = cols
  unit = cols[unit]
  time = cols[time]
  outcome = yvar= cols[outcome]
  treatment = cols[treatment]
  if (missing(x)){
    stop("Please specify the names or column numbers of the time varying control variables in the argument x.")
  }
  xvars = org.xvars = cols[x]

  # Create a model matrix so that the function also works
  # if x contains categorical variables
  xform = as.formula(paste0("~",paste0(xvars, collapse="+")))
  xmat = model.matrix(xform, panel)[,-1, drop=FALSE]

  org.panel = panel
  if (!setequal(colnames(xmat), xvars)) {
    panel = cbind(panel[,setdiff(colnames(panel),xvars)], xmat)
    xvars = colnames(xmat)
  }

  if (is.null(x.rows)) {
    x.rows =as.integer(panel[[treatment]]) == 0
  }
  if (is.character(x.rows)) {
    x.rows = panel[[x.rows]]
  }

  xpanel = panel[x.rows,]

  if (n_distinct(xpanel[[treatment]]) >1) {
    warning("Your choice of x.rows also contains obvserations where treatment takes place. You should estimate the effect of x on y using observations wthout treatment.")
  }

  # 1. Estimate relationship between x and y from pre-experimental periods
  form = as.formula(paste0(outcome, "~", paste0(xvars,collapse=" + "), "|", unit, " + ", time))
  reg = feols(form, data = xpanel)

  # 2. Subtract x effect from y to computed adjusted y
  beta = coef(reg)[seq_along(xvars)]
  xeff = as.numeric(xmat %*% beta)

  y.adj = org.panel[[outcome]] - xeff + add.mean*mean(xeff)
  if (NROW(y.adj) != NROW(panel)) {
    stop("Your sample contained NA in the relevant variables")
  }
  y.adj
}
