# antivir2017
Supplemental simulation code for the paper [link to bioarix] 

# Python code for evaluating neuraminidase inhibitor dynamics

[Link to Python Notebook???]

# R code for epidemic simulations

## Default settings 

Model
```R
atv <- function(t, x, pars) {
    with(as.list(c(pars, x)), {
        n   <- 1
        Kr  <- EC50
        fR  <- 1/(1+(D/Kr)^n)

        rightV <- ifelse(t >= 1 & V < 1, - c*V, p * fR * I - c*V)
        
        dT = -b * T * V
        dJ = b * T * V - k * J
        dI = k * J - delta * I
        dV = rightV
        dP = theta * I - a * P
        dD = -g * D
        list(c(dT,dJ,dI,dV,dP,dD))
    })
}
```


Parameters
```R
pars <- c(b     = 0.0674, k    = 3.684, delta = 1.364,
          p     = 40356,  EC50 = 42.3,  c     = 8,
          theta = 2.75,   a    = 0.498, g     = 3.26)
```

## Network model's properties

Degree distribution

```R
hist(degree(dget("M")$g), breaks=50, main="Network degree distribution", xlab="Number of contact", col="slategrey", sub="POLYMOD data (Mossong et al. (2008) Plos Med)")
```

 ![](./fig/degree.png)

## Parallel without cluster computer

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

## Simulation scenarios

```R
# Control no vaccination --------------------------------------------------
  runcon()
# 5 days regimen ----------------------------------------------------------
  # Generating D dynamics
  x1  <- c(D = 0); .t  <- seq(0, 365, .01)
  o1  <- ode(x1, .t, drug, pars, events=list(data=ij(0, .5, 5, 150)))
  # Put into a light form, changes for every scenario
  fD  <- approxfun(o1[, "time"], o1[, "D"], rule=2)
# 1 weeks .9 coverage -----------------------------------------------------
  # Generating D dynamics
  x1  <- c(D = 0); .t  <- seq(0, 365, .01)
  totalDrug <- 7*300*.9*10000
  o1  <- ode(x1, .t, drug, pars, events=list(data=ij(0, .5, 7, 150)))
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
