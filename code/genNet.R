# POLYMOD Data
# Table 1 POLYMOD
agecont <- orderBy(data.frame(
  age=c(10.21, 14.81, 18.22, 17.58, 13.57, 14.14, 13.83, 12.30, 9.21, 6.89), 
  grp=c("0–4", "5–9", "10–14", "15–19", "20–29", "30–39", "40–49", "50–59", "60–69", "70+"), 
  nmr=1:10), 1) 
agecont$rk <- rank(agecont$age)
# Age grouped by 5-years
agen <- 1:15
ages <- paste(seq(0, 70, 5), c(seq(4, 69, 5), "+"), sep = "-")
# A target age distribution
exampleAge <- dget("./data/exampleAge")
rAge       <- exampleAge$age
pAge       <- exampleAge$p
# Contact distribution
POLYMODcontact <- dget("./data/POLYMODcontact")
nAllunique     <- POLYMODcontact$n
pAll           <- POLYMODcontact$p
# Contact matrix
pages <- dget("./data/contactMatrix")
# Generating network
genNet <- function(N = 100, seed=NULL, freq = X, prob = Y, getGraph=TRUE, n.try=100) {
  # Freq: target contact frequency
  # Prob: target contact prob.
  # =======================================================================
  # Generating network
  # -----------------------------------------------------------------------
  adjmat         <- Matrix::Matrix(0, N, N) # population matrix
  contacts       <- vector("list", N)  # contacts container for each node
  freeId         <- id <- 1:N  # who are free for making contact
  if (!is.null(seed)) set.seed(seed)
  agei           <- sample(rAge, N, replace=1, prob=pAge)  # assign age
  # detailing the nodes in group 70+
  l70            <- length(agei[agei==70])
  agei[agei==70] <- sample(70:80, l70, replace=1)
  ncont          <- sort(sample(freq, N, replace=1, prob=prob), TRUE)
  # new break for 70+ fAge function
  newbreak <- c(0,4,9,14,19,29,39,49,59,69,80) 
  ageid    <- cut(agei, breaks=newbreak, include.lowest=1)
  ageid    <- as.numeric(ageid)
  agerk    <- sapply(ageid, function(x) agecont$rk[x==agecont$nmr] )
  agei     <- rev(orderBy(cbind(agei, agerk), 2)[, 1])
  ageid[is.na(ageid)] <- 15
  freeId   <- id[which(ncont>0)]  # some have no contact?
  # =======================================================================
  # Utils functions
  # -----------------------------------------------------------------------
    fillMat <- function(mat, vec, i) {
      cUP <- vec[vec > i]  # check to fill row
      cLO <- vec[vec < i]  # check to fill col
      if (length(cUP)>0) mat[i, cUP] <- 1 
      if (length(cLO)>0) mat[cLO, i] <- 1
      return(mat)
    }
    pickAge <- function(i) sample(rev(agen), ncont[i], prob=pages[, ageid[i]], replace=1)
    pickId  <- function(k) sample(freeId[ageid[freeId]==apick[k]], npick[k])
    toGraph <- function(mat, toSym=TRUE) {
      # Convert upper.tri adjacent matrix to graph
      if (toSym) mat[lower.tri(mat)] = t(as.matrix(mat))[lower.tri(mat)]
      igraph::graph_from_adjacency_matrix( mat, mode="undirected", diag = FALSE )
    }
    toMat <- function(mat) {
      # Convert upper.tri adjacent matrix to full matrix
      mat[lower.tri(mat)] = t(as.matrix(mat))[lower.tri(mat)]
      return(mat)
    }
  # =========================================================================
  # Main loop
  # -----------------------------------------------------------------------
  for (i in 1:N) {
    cat('\r', "Node:", i)
    flush.console()
    if (!i %in% freeId) next
    # Checking if there is enough peoples to contact
    freeA <- rle(sort(ageid[freeId]))  # current free nodes age
    nfree <- freeA$lengths  # number of free nodes by ages
    afree <- freeA$values # ages that have free nodes
    isEnough <- noWay <- FALSE
    i.try    <- 0
    while(!isEnough) {  # trying to find matched age partners
      i.try <- i.try + 1
      if ( i.try > n.try ) isEnough <- noWay <- TRUE # force stop
      conti <- rle(sort(pickAge(i)))  # pick age of i's contacts
      npick <- conti$lengths  # number of contact by ages
      apick <- conti$values  # age need to pick
      if (!all(apick %in% afree)) next  # if no free age, try again
      temp  <- NULL  # if there are free age, check number if enough
      for (j in 1:length(npick) ) {
        temp <- c(temp, npick[j] <= nfree[which(afree==apick[j])])
      }
      if ( all(temp) ) isEnough <- TRUE  # stop normal
    }
    if (noWay) next  # leave this guy alone
    for (k in 1:length(npick)) {  # pick id from free nodes
      contacts[[i]] <- c(contacts[[i]], pickId(k)) 
    }
    adjmat  <- fillMat(adjmat, contacts[[i]], i)  # mark connections
    ncont[i] <- 0  # update me, no more friend
    ncont[contacts[[i]]] <- ncont[contacts[[i]]] - 1  # update friends
    freeId  <- id[which(ncont>0)]
  }
  if (getGraph) return(list(g=toGraph(adjmat), age=agei, ageid=ageid, ct=contacts))
    else return(list(g=toMat(adjmat), age=agei, ageid=ageid, ct=contacts))
}
genNet <- compiler::cmpfun(genNet)
# Plotting networks
plotNet <- function(M, vs=log(degree(M$g)+2)*2, vc=M$ageid, ly=layout_with_fr, ...) plot(M$g, vertex.size=vs, vertex.label=NA, edge.arrow.mode=0, layout=ly, vertex.color=vc, edge.color="gray90",...)