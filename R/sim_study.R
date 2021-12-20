
#' Useful for monte-Carlo study about different variants of Synth-DID and DID
#'
#' For details look at the code. For an illustration see the vignette.
#' @export
xsdid.mc = function(
    N=100, T=200, tau=50,
    beta.tt = 50, beta.x=50, beta.xc=0, beta.treat=50,
    x.tt.load = 0, xc.x.load = 0,
    x.treat.load.mean = 3/4,x.contr.load.mean = 1/4, x.load.alpha.beta=2,
    xc.treat.load.mean = 3/4,xc.contr.load.mean = 1/4, xc.load.alpha.beta=2,
    ar = 0.97, ar.tt=ar, ar.x=ar, ar.xc = ar,
    return.data = FALSE)
{
  library(dplyr)

  dat = expand.grid(i = 1:N, t = 1:T)

  # Simulate a common AR(1) time trend
  # will be used to create time varying x
  # that affects treatment and control group differently
  tt = as.numeric(arima.sim(n=T,list(ar = ar.tt)))

  tt.x = x.tt.load*tt+as.numeric(arima.sim(n=T,list(ar = ar.x)))
  tt.xc = xc.x.load*tt.x + as.numeric(arima.sim(n=T,list(ar = ar.xc)))

  #plot(time.trend)
  x.load.treat = rbeta(N,x.treat.load.mean*x.load.alpha.beta,(1-x.treat.load.mean)*x.load.alpha.beta)

  x.load.contr = rbeta(N,x.contr.load.mean*x.load.alpha.beta,(1-x.contr.load.mean)*x.load.alpha.beta)

  xc.load.treat = rbeta(N,xc.treat.load.mean*xc.load.alpha.beta,(1-x.treat.load.mean)*xc.load.alpha.beta)

  xc.load.contr = rbeta(N,xc.contr.load.mean*xc.load.alpha.beta,(1-x.contr.load.mean)*xc.load.alpha.beta)


  dat = mutate(dat,
    group = ifelse(i > N/2,"treat","control"),
    treat = 1L*(group == "treat"),
    exp = 1L*(t > T/2),
    treat_exp = exp*treat,
    time.trend = tt[t],
    eps = rnorm(n()),
    x.load = ifelse(treat==1, x.load.treat[i], x.load.contr[i]),
    xc.load = ifelse(treat==1, xc.load.treat[i], xc.load.contr[i]),
    x = x.load * tt.x[t],
    xc = xc.load * tt.xc[t],
    y = beta.tt * time.trend + treat*beta.treat + treat_exp*tau + beta.x*x + beta.xc*xc+ eps
  ) %>% arrange(treat_exp, treat, exp, i, t)

  if (return.data) return(dat)

  # Vanilla DID with FE
  did = feols(y~treat_exp|i+t, data=dat)
  did.beta = coef(did); names(did.beta) = NULL
  didx = feols(y~treat_exp+x|i+t, data=dat)
  didx.beta = coef(didx); names(didx.beta) = NULL

  library(xsynthdid)
  dat$y.adj = adjust.outcome.for.x(dat, unit="i",time = "t",outcome = "y",treatment = "treat_exp",x="x")
  dat$y.res = resid(lm(y~x, data=dat))

  pm = panel.matrices(dat, unit="i",time = "t",outcome = "y",treatment = "treat_exp")
  sdid = synthdid_estimate(Y=pm$Y,N0=pm$N0,T0 = pm$T0)

  pm = panel.matrices(dat, unit="i",time = "t",outcome = "y.adj",treatment = "treat_exp")
  sdidx = synthdid_estimate(Y=pm$Y,N0=pm$N0,T0 = pm$T0)

  pm = panel.matrices(dat, unit="i",time = "t",outcome = "y.res",treatment = "treat_exp")
  sdid.res = synthdid_estimate(Y=pm$Y,N0=pm$N0,T0 = pm$T0)


  sc = sc_estimate(Y=pm$Y,N0=pm$N0,T0 = pm$T0)
  sc

  data.frame(sdidx = as.numeric(sdidx), sdid = as.numeric(sdid), sdid.res = as.numeric(sdid.res), did=did.beta[1], didx = didx.beta[1], sc=as.numeric(sc),N=N, T=T, tau=tau, beta.treat=beta.treat,
  beta.x=beta.x, beta.tt = beta.tt, beta.xc = beta.xc,
  x.tt.load = x.tt.load, xc.x.load = xc.x.load
  )
}

