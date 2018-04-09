shaded_contour  <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1,
    length.out = ncol(z)), z, xlim = x[c(1,length(x))],
    ylim = y[c(1,length(y))], zlim = range(z, finite = TRUE),
    levels = pretty(zlim, nlevels), nlevels = 20, col = rbfun(length(levels) - 1), 
    xaxs='i', yaxs='i', xlab, ylab, las=1, 
    borders=TRUE, l.col=grey(0.7), l.lwd=1,
    block.na=0, na.col='lightblue', add=FALSE, ...){

    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                if (is.null(dim(x)))
                  stop("argument must be matrix-like")
                z <- x
                x <- seq.int(0, 1, length.out = nrow(z))
            }
            if (missing(xlab))
                xlab <- ""
            if (missing(ylab))
                ylab <- ""
        }
        else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
        xn <- deparse(substitute(x))
        if (missing(xlab))
            xlab <- paste(xn, "x", sep = "$")
        if (missing(ylab))
            ylab <- paste(xn, "y", sep = "$")
        y <- x$y
        x <- x$x
    }
    else {
        if (missing(xlab))
            xlab <- if (missing(x))
                ""
            else deparse(substitute(x))
        if (missing(ylab))
            ylab <- if (missing(y))
                ""
            else deparse(substitute(y))
    }
    
    if (all(diff(x) <= 0)){
        x   <- rev(x)
        z   <- z[nrow(z):1,]    
    }
    if (all(diff(y) <= 0)){
        y   <- rev(y)
        z   <- z[,ncol(z):1]
    }

    par(las = las)
    
    if (block.na == 2){
        ylim[which.min(ylim)] <- max(min(ylim), -84.5)
    }
    
    if (!add) plot(NA, NA, xlim = xlim, ylim = ylim, type = "n", xaxs = xaxs,
        yaxs = yaxs, xlab = xlab, ylab = ylab, ...)
    if (!is.matrix(z) || nrow(z) <= 1L || ncol(z) <= 1L)
        stop("no proper 'z' matrix specified")
    if (any(is.na(z)) & block.na){
        xx      <- rep(x, length(y))
        yy      <- rep(y, each=length(x))
        mask    <- !is.na(z)
        zz      <- array(interpp(xx[mask], yy[mask], z[mask], xo=xx, yo=yy)$z, dim(z))
        if (!is.double(zz))
            storage.mode(zz) <- "double"
    } else {
        if (!is.double(z))
            storage.mode(z) <- "double"
    }
    .Internal(filledcontour(as.double(x), as.double(y), 
        if (any(is.na(z)) & block.na) zz else z, as.double(levels),
        col = col))
    if (borders){    
         contour(x,y,z, lev=levels, drawlabels=F, add=TRUE, col=l.col, lwd=l.lwd)
         contour(x,y,z, lev=0, drawlabels=F, add=TRUE, col=l.col, lwd=2*l.lwd)
    }
    if (any(is.na(z)) & block.na == 1){
        image(x,y,mask*1, lev=seq(-0.5,1.5), col=c(na.col, NA), add=T)
    } else if (block.na == 2){
        map(add=T, fill=TRUE, col=na.col)
    }
    box()
    
    out <- list(col=col, lev=levels)
    class(out) <- "plotmap"
    
    invisible(out)
}
