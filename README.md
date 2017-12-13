# antivir2017
Supplemental simulation code for the paper [link to bioarix] 

```R
source("./code/helpers.R")
source("./code/ModelandParams.R")
source("./code/genNet.R")
source("./code/networkSim.R")
pop <- 10^4
```
## Network generation
```R
# -------------------------------------------------------------------------
# (not run), run once and load from RData instead
# -------------------------------------------------------------------------
# M <- genNet(10000, 123, nAllunique, pAll)
M = readRDS("./data/M.Rds")
hist(degree(M$g), breaks=50, main="Network degree distribution", xlab="Number of contact", col="slategrey", sub="POLYMOD data (Mossong et al. (2008) Plos Med)")
# savePNG("./fig/degree")
```

![](./fig/degree.png)

## Parallel without cluster version

```R
runcon <- function(namex="control", pvaci=0, ssize=1:1000) {
  dir.create(namex)
  ncores <- detectCores()
  loopF  <- function(x) {
    out  <- antivir(ni0=100, pvac=0)
    dput(out, paste0(namex, "/", pvaci*100, "_", x) )
  }
  system.time(mclapply(ssize, loopF, mc.cores=ncores))
}
```
## Scenarios

```R
# -------------------------------------------------------------------------
# Long computation process
# -------------------------------------------------------------------------
# Control no vaccination --------------------------------------------------
  runcon()
# 1 weeks .9 coverage -----------------------------------------------------
  totalDrug <- 7*300*.9*10000
  # Generating D dynamics
  x1  <- c(D = 0); .t  <- seq(0, 365, .01)
  o1  <- ode(x1, .t, drug, pars, events=list(data=ij(0, .5, 7, 150)))
  # Put into a light form, change for every scenario
  fD  <- approxfun(o1[, "time"], o1[, "D"], rule=2)
  runcon("1week90", 0.9)
# 2 weeks .45 coverage ----------------------------------------------------
  (totalDrug/(14*300))/10000
  o1  <- ode(x1, .t, drug, pars, events=list(data=ij(0, .5, 14, 150)))
  fD  <- approxfun(o1[, "time"], o1[, "D"], rule=2)
  runcon("2week45", 0.45, 1:100)
# 3 weeks .30 coverage ----------------------------------------------------
  (totalDrug/(21*300))/10000
  o1  <- ode(x1, .t, drug, pars, events=list(data=ij(0, .5, 21, 150)))
  fD  <- approxfun(o1[, "time"], o1[, "D"], rule=2)
  runcon("3week30", 0.3, 1:100)
# 4 weeks .225 coverage ---------------------------------------------------
  (totalDrug/(28*300))/10000
  o1  <- ode(x1, .t, drug, pars, events=list(data=ij(0, .5, 28, 150)))
  fD  <- approxfun(o1[, "time"], o1[, "D"], rule=2)
  runcon("4week225", .225, 1:100)
# 5 weeks .18 coverage ----------------------------------------------------
  (totalDrug/(35*300))/10000
  o1  <- ode(x1, .t, drug, pars, events=list(data=ij(0, .5, 35, 150)))
  fD  <- approxfun(o1[, "time"], o1[, "D"], rule=2)
  runcon("5week18", .18, 1:100)
# 6 weeks .15 coverage ----------------------------------------------------
  (totalDrug/(42*300))/10000
  o1  <- ode(x1, .t, drug, pars, events=list(data=ij(0, .5, 42, 150)))
  fD  <- approxfun(o1[, "time"], o1[, "D"], rule=2)
  runcon("6week15", 0.15, 1:100)
closeAllConnections()
```