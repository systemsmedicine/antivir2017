# Main model
# -------------------------------------------------------------------------
atv <- function(t, x, pars) {
    with(as.list(c(pars, x)), {
        n   <- 1 # w   <- 75 e   <- .95 tau <- .5
        Kr  <- EC50
        fR  <- 1/(1+(D/Kr)^n)

        rightV <- ifelse(t >= 1 & V < 1, - c*V, p * fR * I - c*V)
        
        dT = -b * T * V
        dJ = b * T * V - k * J
        dI = k * J - delta * I
        dV = rightV
        dP = theta * I - a * P
        dD = -g * D
        list(c(dT, dJ, dI, dV, dP, dD))
    })
}

# Model for epidemic simulation
# -------------------------------------------------------------------------
etv <- function(t, x, pars) {
    with(as.list(c(pars, x)), {
        D <- fD(t+A)  # shifting the drug dynamic!!!
        fR  <- 1/(1+(D/42.3))
        rightV <- ifelse(t >= 1 & V < 1, - 8*V, 40356 * fR * I - 8*V)
        dT = -0.0674 * T * V
        dJ = 0.0674 * T * V - 3.684 * J
        dI = 3.684 * J - 1.364 * I
        dV = rightV
        dP = 2.75 * I - 0.498 * P
        list(c(dT, dJ, dI, dV, dP))
    })
}

# Drug dynamics
# -------------------------------------------------------------------------
drug <- function(t, x, pars) {
    with(as.list(c(pars, x)), {
        dD = -g * D
        list(c(dD))
    })
}

# EC50 to efficacy
# -------------------------------------------------------------------------
# D0 <- EC50 * ( epsilon/(1 - epsilon) ) * (1 - exp(-lambda * tau))
# D0 / ( EC50 * (1 - exp(-lambda * tau)) ) = epsilon/(1 - epsilon)
# epsilon  = x / (1 + x )
ef <- function(D0 = 1, EC50 = 42, lambda = 3.26, tau = 1) {
    x <- D0 / ( EC50 * (1 - exp(-lambda * tau)) )
    return( x / (1+x) )
}

# Parameters
# -------------------------------------------------------------------------
pars <- c(b = 0.0674, k = 3.684, delta = 1.364, p = 40356, EC50 = 42.3, c = 8, theta = 2.75, a = 0.498, g = 3.26)

# Model settings
# -------------------------------------------------------------------------
x0 <- c(T = 1, J = 0, I = 0, V = 7.5e-6, P = 0, D = 0)
.t <- seq(0, 21, .1)

# Model without any interventions
# -------------------------------------------------------------------------
o <- ode(x0, .t, atv, pars)
fNorm <- approxfun(.t, o[, "V"], rule=2)
# plot(o)
