`set_position`  <- function(x, y, size=0.1, bg=NULL){
    # sets the display for a plot centered at x,y
    .old_par    <<- par(no.readonly=TRUE)
    lims        <- par("usr") 
    flims       <- par("plt")

    # expand size if necessary
    if (length(size) %in% c(1,2)){
        size    <- c(-0.5,0.5,-0.5,0.5)*rep(size, each=2)
    } 
    
    if ( x > lims[1] & x < lims[2] & y > lims[3] & y < lims[4]){
        # change to plot coordinates
        figx    <- (x-lims[1])/diff(lims[1:2])*diff(flims[1:2]) + flims[1]
        figy    <- (y-lims[3])/diff(lims[3:4])*diff(flims[3:4]) + flims[3]
        
        # compute the bounds
        bounds  <- size*diff(flims)[c(1,1,3,3)] + rep(c(figx, figy), each=2)        
        
        # set the background to the background colour
        if (!is.null(bg)){
            b.xy    <- c(x,y,x,y) + size[c(1,3,2,4)]*diff(lims)[c(1,3,1,3)]
            rect(b.xy[1], b.xy[2], b.xy[3], b.xy[4], border=NA, col=bg)
        }
                
        par(plt=bounds, new=T)

    } else {
        stop("Desired position not within plot region")
    }
}

`reset_position` <- function(){
    par(.old_par)
}


