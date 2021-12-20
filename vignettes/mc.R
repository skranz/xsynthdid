
# Runs simulation for chapter 3 in Vignette
run.mc = function() {
  library(dplyr)
  library(xsynthdid)

  grid = expand.grid(N=c(20,50,200),T=c(20,50,200), beta.xc=c(0,50),
                     x.tt.load = c(0,1), xc.x.load = c(0)) %>%
    filter(N == T) %>%
    filter(beta.xc == 0 | (x.tt.load == 0))


  repl = 1000
  i = 1
  res = bind_rows(lapply(1:NROW(grid), function(i) {
    #if (i > 2) return(NULL)
    cat("\n i=",i)
    par = as.list(grid[i,])
    inner = bind_rows(lapply(1:repl,function(r) {
      do.call(xsdid.mc,par)
    }))
    inner$scenario = i
    inner
  }))

  saveRDS(res,"xsdid_mc.Rds" )

}
