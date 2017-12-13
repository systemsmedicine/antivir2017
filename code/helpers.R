# Customizations
# options(error=recover) # for debug
# -------------------------------------------------------------------------
libs <- c("mgcv", "foreach", "doMC", "igraph", "Matrix", "deSolve", "compiler", "parallel", "devtools", "txtplot")
sapply(libs, library, character.only=TRUE)
# devtools::load_all("~/smidR") # personal package
par(pch=20, lwd=1.2, font.main=1, bty="n")
solarized <- c("#002b36", "#dc322f", "#b58900", "#268bd2", "#859900", "#6c71c4", "#d33682", "#2aa198")
palette(solarized)
# Treatment regimen generation
# -------------------------------------------------------------------------
ij <- function(t0=0, freq=1, duration=5, x=75, name="D") {
	maxt <- ifelse(is.integer(freq), t0+(duration-1), t0+(duration-1)+freq)
	time <- round(seq(t0, maxt, by=freq), 1)
	n    <- length(time)
	o    <- data.frame(var=rep(name, n), time=time, value=rep(x, n), method=rep(2, n))	
	return(o)
}

# Semi steady state of drug
# -------------------------------------------------------------------------
cesar.hi <- function(D = 1, day = 1, tau = 1) D/(1 - exp(-log(2)/day*tau))
l <- function(x) log10(x)
# AUC calculation
# -------------------------------------------------------------------------
AUC <- function(x, y, maxX=length(x)) {
  x <- x[0:maxX]
  y <- y[0:maxX]
  return(sum(diff(x) * (head(y,-1)+tail(y,-1)))/2)
}
find50   <- function(x, y=0, z="V") {
	return(AUC(x[which(x[, z] >= y), "time"], x[which(x[, z] >= y), z]))
}

# Utilities functions for deSolve
# -------------------------------------------------------------------------
vtol     <- 1e-10
rootfun  <- function (t, x, pars) {
	with (as.list(x), {
			return( max(0, V - vtol) )
	  })
}
eventfun <- function(t, x, pars) {
  with (as.list(x),{
    V <- 0
    return(c(T, J, I, V, P, D))
  })
}
orderBy <- function(.data,index) as.data.frame(.data[order(.data[, index]), ])

# Processing output functions
# -------------------------------------------------------------------------

# Getting infected and recovery dynamics
processIR <- function(folder="./control") {
	# Getting infected and recovery dynamics
  IR <- lapply(seq_along(list.files(folder)), function(x) {
    out <- dget(paste0(folder,"/", list.files(folder)[x]))
    rbind(out$SIR$I/pop, out$SIR$R/pop)
  })
  return(IR)
}

# Getting reproduction number
process <- function(folder="./control", full=TRUE) {
	# Getting reproduction number
  R0 <- lapply(seq_along(list.files(folder)), function(x) {
    out <- dget(paste0(folder,"/", list.files(folder)[x]))
    R0net(out, All=full)
  })
  return(R0)
}
