# Generate infection dynamics
# -------------------------------------------------------------------------
getV <- function(drug=0, A = 0, times=seq(0, 21, .01)) {
  # A: time (how long has been passed before the infection)
  if (drug==0) return(fNorm)
    else {
      temp <- ode(xe, times, etv, c(A=A))[, "V"]
      return(approxfun(times, temp, rule=2))
    }
}
# Example getV(1, 5)

# Epidemic with antivir simulations
# -------------------------------------------------------------------------
antivir <- function(net=M, ni0=1, Vc=1.35, maxV=2642.228, tmax=365, tStep=1, maxi=21, pvac=0.3, tVac=0) {
  require(Matrix); require(deSolve)
  `%notin%` <- function(x, y) is.na(match(x, y, nomatch=NA_integer_))
    age   <- net$age #TODO generalize for any given age distribution
    Mx    <- net$g  # the graph of mine
    pop   <- length(age)
    t.inf <- seq(0, maxi, .01)
    t.epi <- seq(0, tmax, tStep)
    vlist <- vector("list", pop)
    SIR   <- Matrix(0, length(t.epi), 4) #todo optional storing
    colnames(SIR) <- c("S", "I", "R", "In")
    S <- Matrix(0, pop, 4)
    colnames(S) <- c("Tr", "Ti", "s", "Tv")
  Fate <- function(who) {
    myV   <- vlist[[who]](t.inf)
    tVmax <- t.inf[which.max(myV)]
    tR    <- t.inf[myV < Vc & t.inf >= tVmax][1]
    return( c(tR) )
  }
  updateI <- function(who, tinf) {
    S[who, "Ti"] <<- tinf # record time of infection
    S[who, "s"]  <<- 1    # change status
    S[who, "Tr"] <<- sapply(who, Fate)
  }
  updateStatus <- function(day) {
    SIR[day, "S"] <<- length(which(S[,"s"]==0))
    SIR[day, "I"] <<- length(which(S[,"s"]==1))
    SIR[day, "R"] <<- length(which(S[,"s"]==2))
  }
  neighbors <- function(mat, r) return(which(mat[r, ]==1))
  findNaive <- function(mat, i) {
    nbi  <- neighbors(mat, i)           # Find all neighbor(s)
    nbi  <- nbi[which(S[,"s"][nbi]==0)] # who is susceptible
    return(nbi)
  }
  ids <- 1:pop
    i0 <- sample(ids, ni0) # choose infected and update status
    nn <- length(neighbors(Mx, i0))
    message("\r", length(i0), " case infected, total ", nn, " contacts!\n")
    invisible(lapply(i0, function(x) vlist[[x]] <<- getV(FALSE)))
    updateI(i0, tinf=0)
    SIR[1, "In"] <- length(i0) # incidences as of now
    if (pvac!=0) {
      vacID <- sample(ids[-i0], round(pop*pvac)) # Got drug IDs
      S[vacID, "Tv"] <- 1 # start drug time
    }
  pId <- vector('character', pop)
  # Main loop -----------------------------------------------------------
  for (ti in 2:length(t.epi)) {
    now <- t.epi[ti]
    updateStatus(ti) # tracking numbers
    cat('\r', "time:", now, "S I R In", as.vector(SIR[ti-1, ]))
    flush.console()
    # Infection processes -----------------------------------------------  
    infId <- which(S[,"s"] %in% c(1) )  # Who infected/death-not-buried
    if ( length(infId) < 1 ) {
      SIR[(ti+1):length(t.epi), "S"] <- SIR[ti, "S"] # current S
      SIR[(ti+1):length(t.epi), "R"] <- SIR[ti, "R"] # current R
      break  # stop if no one infected
    } 
    infDuration <- now - S[infId, "Ti"]  # infected duration till now
    currentV <- mapply(do.call, vlist[infId], lapply(infDuration, list) )
    # Transmission potential: 
    pTrans <- pmin(1, currentV/maxV) 
    # If log, can be <0
    # pTrans <- ifelse (l(currentV) < 0, 0, l(currentV)/l(maxV))
    # pTrans <- pmin(1, pTrans)                                           
    # --------------------------------------------------------------------
    nbi <- lapply(infId, findNaive, mat=Mx)  # Find susceptible neighbor
    # Age weighting here, temporatily disable **************************
    # ageW <- lapply(nbi, function(x) fAge(age[x]))  # Weights age
    # pW <- Map("*", ageW, sapply(pTrans, list) ) # Weights age
    # p <- Map(rbinom, sapply(nbi, length), size=1, pW) # try to infect
    # --------------------------------------------------------------------
    p <- Map(rbinom, sapply(nbi, length), size=1, pTrans) # try to infect
    P <- unlist(p) == 1  # P[c(2,4)] <- TRUE # for debugging
    if ( any(P) ) {
      allnb         <- unlist(nbi)
      Yes           <- allnb[P]  # infected ids
      SIR[ti, "In"] <- length(Yes)  # incidences as of now
      tinf          <- now + runif(length(Yes), -.5, .5) # diffuse tinf
      v             <- S[Yes, "Tv"]
      vlist[ Yes ]  <- Map(getV, drug=v, A=(tinf-1)) # because start day 1
      updateI(Yes, tinf)
    } else {
      SIR[ti, "In"] <- 0  # incidences as of now
    }
    # Update the recover ------------------------------------------------  
    infId             <- which(S[,"s"] == 1 )  # Who infected
    infDuration       <- now - S[infId, "Ti"]
    recoverConditions <- infDuration >= S[infId, "Tr"]
    if (any(recoverConditions)) {
      recoverID         <- infId[recoverConditions]
      S[recoverID, "s"] <- 2
    }
  }
  # Storing for output ----------------------------------------------------
  SIR <- as.data.frame( cbind(t.epi, as.matrix(SIR)) )
  colnames(SIR)[1]  <- "time"
  class(SIR)        <- "smidSIR"
  return(list(SIR=SIR, aux=S[, "s"], nn=nn)) 
}
antivir <- compiler::cmpfun(antivir)

# calculate R0 from sims outputs
# -------------------------------------------------------------------------
R0net <- function(out = sim, net=M, All=FALSE) {
    if (class(out$SIR) != "smidSIR") stop("Incompatible")
    net <- net$g
    n  <- igraph::delete.vertices(net, which(out$aux == 0))
    n  <- igraph::as.directed(n, mode = c("arbitrary"))
    dg <- igraph::degree(n,  mode='out')
    if (All==TRUE) r0 <- dg
      else r0 <- mean(dg)
    return(r0)
}

# Plot SIR and In from output
# -------------------------------------------------------------------------
plotSIR <- function(x=out, startTime=2, endTime=length(x$aux),...) {
  # startTime = 2 due to at the beginning there were already some infected
  par(pch=20, lwd=1.2, font.main=1, bty="n")
  with(x$SIR, {
    lineplot(time[startTime:endTime], S[startTime:endTime], ylim=c(0, length(x$aux)),...)
    lines(time[startTime:endTime], I[startTime:endTime], col=2)
    lines(time[startTime:endTime], R[startTime:endTime], col=3)
    lines(time[startTime:endTime], In[startTime:endTime], col=4)
  })
} 
