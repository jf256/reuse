`plot_colourbar` <-
function (x, ...){
    UseMethod("plot_colourbar")
}

`plot_colourbar.default` <-
function(levs, cols, side=1, ylab="", labels=NULL,
    xlab="", nn=1, center=F, cex.axis=1, sea.col=NULL,...){
    # plots a colour bar given the levs and cols
    # centers the axis labelling instead of at the 
    # boundary when center==TRUE
    # if sea.col is set, the colourbar is prolonged by one
    # colour with the uncentered label sea underneath

    sea.add     <- 0
    if (!is.null(sea.col)) {
        cols    <- c(sea.col, cols)
        sea.add <- 1
    }
    ncols       <- length(cols)
    lev.arr     <- array(1:ncols, c(1, ncols))
    if (side %in% c(1,3)){
        lev.arr <- t(lev.arr)
    }    
    
    if (side %in% c(1,3)){
        image(1:ncols, 1, lev.arr, axes=F, breaks=1:(ncols+1)-0.5, col=cols,
            ylab=ylab, xlab=xlab, ...)
        abline(v=1:(ncols-1) + 0.5)
    } else {
        image(1, 1:ncols, lev.arr, axes=F, breaks=1:(ncols+1)-0.5, col=cols,
            ylab=ylab, xlab=xlab, ...)
        abline(h=1:(ncols-1) + 0.5)
    }
    if (center){
        at.lev  <- seq(1+sea.add,ncols,nn)
        if (is.null(labels)){
            labels <- (levs[1:(ncols-sea.add)] + levs[2:(ncols+1-sea.add)])/2
        }
        axis(side, at=at.lev, labels=labels, las=1, tick=F, cex.axis=cex.axis)
    } else {
        at.lev  <- seq(2+sea.add,ncols,nn)
        if (is.null(labels)){
            labels  <- levs[at.lev - sea.add]
        }
        axis(side, at=at.lev-0.5,labels=labels, las=1, cex.axis=cex.axis) 
    }
    if (!is.null(sea.col)){
        axis(side, at=1, labels="sea", las=1, tick=F, cex.axis=cex.axis)
    }
    box()
}


`plot_colourbar.plotmap` <-
function(x, incl.units=T, side=1, cex.axis=1, center=F, labels=NULL, ...){
    if (!is.null(x$flag_values)){
        center  <- TRUE
        labels  <- x$flag_values
        if (!is.null(x$flag_meanings)){
            labels  <- x$flag_meanings
        }
    }
        
    plot_colourbar(x$lev, x$col, sea.col=x$sea.col, side=side, 
        cex.axis=cex.axis, center=center, labels=labels, ...)

    # plot unit-string in axis
    if (incl.units){
        if (!is.null(x$longname)) {
            if (!is.null(x$units)){
                unit.string <- paste(x$longname," (",x$units, ")", sep="")
            } else {
                unit.string <- x$longname
            }
            nlev    <- length(x$lev)            
            axis(side, at=nlev/2, labels=unit.string, 
                tick=F, las=1, line=2, cex.axis=cex.axis)
        } else {
            if (!is.null(x$units)){
                axis(side, at=length(x$lev)-0.5, labels=x$units, 
                    tick=F, las=1, cex.axis=cex.axis)
            }
        }
    }
}

