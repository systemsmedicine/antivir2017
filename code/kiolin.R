kiolin <- function (x, range = 1.5, h = NULL, method="normal", ylim = NULL, names = NULL, horizontal = FALSE, col = "magenta", border = "black", lty = 1, lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, at, add = FALSE, wex = 1, drawRect = TRUE, xlab="xLab", ylab="yLab", main="") {
    datas <- x
    n <- length(datas)
    if (missing(at)) at <- 1:n
    upper <- lower <- q1 <- q3 <- med <- vector("numeric", n)
    base <- height <- vector("list", n)
    baserange <- c(Inf, -Inf)
    args <- list(display = "none", method=method)
    if (!(is.null(h))) args <- c(args, h = h)
    for (i in 1:n) {
        data <- datas[[i]]
        data.min <- min(data)
        data.max <- max(data)
        q1[i] <- quantile(data, 0.25)
        q3[i] <- quantile(data, 0.75)
        med[i] <- median(data)
        iqd <- q3[i] - q1[i]
        upper[i] <- min(q3[i] + range * iqd, data.max)
        lower[i] <- max(q1[i] - range * iqd, data.min)
        est.xlim <- c(min(lower[i], data.min), max(upper[i], 
            data.max))
        smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
            args))
        hscale <- 0.4/max(smout$estimate) * wex
        base[[i]] <- smout$eval.points
        height[[i]] <- smout$estimate * hscale
        t <- range(base[[i]])
        baserange[1] <- min(baserange[1], t[1])
        baserange[2] <- max(baserange[2], t[2])
    }
    if (!add) {
        xlim <- if (n == 1) 
            at + c(-0.5, 0.5)
        else range(at) + min(diff(at))/2 * c(-1, 1)
        if (is.null(ylim)) {
            ylim <- baserange
        }
    }
    if (is.null(names)) label <- 1:n
        else label <- names
    boxwidth <- 0.05 * wex
    plot.new()
    plot.window(xlim = xlim + c(-0.5, 0.5), ylim = ylim)
    Kaxis(2, las=2)
    Kaxis(1, label = label, at = seq_along(label))
    title(main=main, xlab=xlab, ylab=ylab)
    Kgrid()
    # box()
    for (i in 1:n) {
        polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
            c(base[[i]], rev(base[[i]])), col = col[i], border = border[i],
            lty = lty, lwd = lwd)
        if (drawRect) {
            lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
              lty = lty)
            rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
              q3[i], col = rectCol[i])
            points(at[i], med[i], pch = pchMed, col = colMed[i])
        }
    }
    invisible(list(upper = upper, lower = lower, median = med, 
        q1 = q1, q3 = q3))
}
