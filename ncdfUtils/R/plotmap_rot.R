`plotmap_rot` <-
function(file, file.small=NULL, sponge=8, sponge.small=15, 
    varname = NULL, lsm.file=NULL, lsm.file.small=NULL, 
    col=NULL, levels=NULL, sea.col=NULL, rivers=T, 
    cities=T, label=TRUE, minpop=NULL, ncities=10, city.pch=19,
    alt.contour=F, alt.lev=NULL,
    grid=TRUE, grid.txt=TRUE, grid.lty=2,
    i.time=1, i.lev=1, map.lwd=2, 
    cex.axis=1, cex.lab=1, cex.main=1, cex.txt=1,
    main="", xlab="", ylab="", add=FALSE,
    colourplot=TRUE, hires=FALSE, interior=FALSE, alt.poli=TRUE,
    nlongrid=10, nlatgrid=5, lon.ind, lat.ind){
    
    # load hi-resolution map-data only when hires is set
    # in order to avoid complications on a system without
    # the mapdata package
    if (hires) library("mapdata")
    
    # read in the data from file
    nc <- open.ncdf(file)

    # read in the coordinates of the rotated pole
    pollon  <- as.numeric(att.get.ncdf(nc, nc$var$rotated_pole$id, 
        "grid_north_pole_longitude")$value)
    pollat  <- as.numeric(att.get.ncdf(nc, nc$var$rotated_pole$id, 
        "grid_north_pole_latitude")$value)
    polgam  <- att.get.ncdf(nc, nc$var$rotated_pole$id, 
        "north_pole_grid_longitude")
    if (polgam$hasatt){
        polgam <- polgam$value
    } else {
        polgam <- 0
    }
    rlon    <- nc$dim$rlon$vals
    rlat    <- nc$dim$rlat$vals
    nlon    <- length(rlon)
    nlat    <- length(rlat)

    # find the variable name (if not set)
    if (is.null(varname)){
        noread  <- c("lon", "lat")
        var.n   <- setdiff(names(nc$var), noread)
        if (any(var.n == "HSURF")){
            varname <- "HSURF"
        } else {
            dims    <- lapply(nc$var[names(nc$var) %in% var.n], function(x) x$size)
            which.v <- sapply(dims, function(x) all(x[1:2] == c(nlon, nlat)))
            varname <- names(which.v)[which.max(which.v)]
        }
    }

    data    <- get.var.ncdf(nc, varname)
    if (length(dim(data)) == 3){
        data    <- data[,,i.time]
    } else if (length(dim(data)) == 4){
        data    <- data[,,i.lev,i.time]
    }
        
    # read in unrotated coordinates
    if (any(names(nc$var) %in% c("lon", "lat"))){
        lon <- get.var.ncdf(nc, "lon")
        lat <- get.var.ncdf(nc, "lat")
    } else {
        tmp <- rot2geo(pollon, pollat, rep(rlon, nlat), rep(rlat, each=nlon), polgam)
            lon <- array(tmp$x, c(nlon,nlat))
            lat <- array(tmp$y, c(nlon,nlat))
    }
    
    
    # store information on variables
    longname    <- nc$var[[varname]]$longname
    units       <- nc$var[[varname]]$units
    flag_values <- att.get.ncdf(nc, nc$var[[varname]]$id,
        "flag_values")
    if (flag_values$hasatt){
        flag_values <- flag_values$value
    } else {
        flag_values <- NULL
    }
    flag_meanings <- att.get.ncdf(nc, nc$var[[varname]]$id,
        "flag_meanings")
    if (flag_meanings$hasatt){
        flag_meanings <- unlist(strsplit(flag_meanings$value, " "))
    } else {
        flag_meanings <- NULL
    }
    
    close.ncdf(nc)

    # mask out the sponge zone
    if (!is.na(sponge) & sponge > 0) {
        data[c(1:sponge,(nlon-sponge+1):nlon),] <- NA
        data[,c(1:sponge,(nlat-sponge+1):nlat)] <- NA
    }

    # if HSURF is plotted, the lsm.file is the same as file
    if (varname == "HSURF"){
        lsm.file    <- file
        if (!is.null(file.small)){
            lsm.file.small  <- file.small
        }
    }

    # read in the land sea mask
    if (is.null(lsm.file)){
        lsm         <- array(TRUE, dim(data))
        alt.contour <- FALSE
    } else {
        nc.lsm  <- open.ncdf(lsm.file)
        lsm     <- get.var.ncdf(nc.lsm, "FR_LAND")
        lsm     <- lsm > 0.5
        if (any(names(nc.lsm$var) %in% c("lon", "lat"))){
            lon <- get.var.ncdf(nc.lsm, "lon")
            lat <- get.var.ncdf(nc.lsm, "lat")
        } else {
            tmp <- rot2geo(pollon, pollat, rep(rlon, nlat), rep(rlat, each=nlon), polgam)
            lon <- array(tmp$x, c(nlon,nlat))
            lat <- array(tmp$y, c(nlon,nlat))
        }
        if (alt.contour & any(names(nc.lsm$var) == "FR_LAND")){
            alt <- get.var.ncdf(nc.lsm, "FR_LAND")
        }
        close.ncdf(nc.lsm)
    }
    if (lon[1,1] > lon[nrow(lon),1] & lon[1,1] < lon[2,1]){
        lon[lon < 0] <- lon[lon < 0] + 360
    }
    
    # read in the data from the nested region
    if (!is.null(file.small)){
        nc.small    <- open.ncdf(file.small)
        data.small  <- get.var.ncdf(nc.small, varname)
        rlon.small  <- nc.small$dim$rlon$vals
        rlat.small  <- nc.small$dim$rlat$vals
        nlon.small  <- length(rlon.small)
        nlat.small  <- length(rlat.small)
        close.ncdf(nc.small)

        if (length(dim(data.small)) == 3){
            data.small    <- data.small[,,i.time]
        } else if (length(dim(data.small)) == 4){
            data.small    <- data.small[,,i.lev,i.time]
        }

        # mask out the sponge zone
        if (!is.na(sponge.small) & sponge.small > 0) {
            data.small[c(1:sponge.small,(nlon.small-sponge.small+1):nlon.small),] <- NA
            data.small[,c(1:sponge.small,(nlat.small-sponge.small+1):nlat.small)] <- NA
        }

        # read in the land sea mask
        if (is.null(lsm.file.small)){
            lsm.small     <- array(TRUE, dim(data.small))
        } else {
            nc.lsm.small  <- open.ncdf(lsm.file.small)
            lsm.small     <- get.var.ncdf(nc.lsm.small, "FR_LAND")
            lsm.small     <- lsm.small > 0.5
            if (alt.contour & any(names(nc.lsm.small$var) == "FR_LAND")){
                alt.small <- get.var.ncdf(nc.lsm, "FR_LAND")
            }
            close.ncdf(nc.lsm.small)
        }
    }
     
    # set levels
    if (is.null(levels)){
        if (varname == "HSURF"){
            levs <- c(-200,0,100,200,500,1000,1500,2000,3000,10000)
        } else {
            if (is.null(flag_values)){
                if (exists("data.small")){
                    levs    <- pretty(c(data, data.small), 20)
                } else {
                    levs    <- pretty(data, 20)
                }
            } else {
                levs    <- approx(seq(along=flag_values), flag_values, 
                    0.5 + 0:length(flag_values), yleft=min(flag_values)-diff(flag_values[1:2]),
                    yright=max(flag_values)+diff(flag_values[length(flag_values) - 0:1]))$y
            }
        }
    } else {
        levs        <- levels
    }

    # set the colours and levels   
    ncols <- length(levs)-1
    if (is.null(col)){
        if (colourplot){
            if (varname == "HSURF"){
                colours	<- .colseq(length(levs)-1, .hsurf, smooth=0)
                sea.col <- .water
            } else if (varname == "SOILTYP"){
                colours <- .soil[flag_values+1]
            } else if (varname %in% c("TOT_PREC", "precip", "pr")){
                colours <- .colseq(length(levs)-1, .gpcc, smooth=0)
            } else {
                colours     <- rbfun(ncols)
            }
        } else {
            colours <- grey((ncols+1):1/(ncols+1))[2:(ncols+1)]
            sea.col <- "white"
        }
    } else {
        colours <- rep(col, length.out=ncols)
    }

    if (hires){
        worlddb <- "worldHires"
        if (any(lon > 180)) worlddb <- "worldHires2"
    } else {
        worlddb <- "world"
        if (any(lon > 180)) worlddb <- "world2"
    }
    if (alt.poli & interior){
      data(polibound)
      world <- polibound
      if (any(lon > 180)){
          world$x[world$x < 0] <- world$x[world$x < 0] + 360
      }
      for (add.name in c(".*Lake.*", ".*Sea.*")){
          world.add <- try(map(worlddb, region=add.name, plot=F, xlim=range(lon), ylim=range(lat)), silent=TRUE)
          if (class(world.add) != 'try-error' & length(world.add) > 0){
              world   <- list(x=c(world$x, NA, world.add$x),
                              y=c(world$y, NA, world.add$y))
        }
      }
       
    } else {
        world       <- map(worlddb,interior=interior, plot=F, 
                         xlim=range(lon), ylim=range(lat))
    }
    if (!interior){
        # remove Lesotho and add the Lakes and Seas
        for (add.name in c(".*Lake.*", ".*Sea.*", ".*Island.*")){
            world.add <- try(map(worlddb, region=add.name, plot=F, 
                xlim=range(lon), ylim=range(lat)), silent=TRUE)
            if (class(world.add) != 'try-error' & length(world.add) > 0){
                world   <- list(x=c(world$x, NA, world.add$x),
                    y=c(world$y, NA, world.add$y), 
                    names=c(world$names, world.add$names))
            }
        }
        world.remove<- map(worlddb, region="Lesotho", plot=F)
        ind.i       <- which(world$x %in% world.remove$x & world$y %in% world.remove$y)
        world$x[ind.i] <- NA
        world$y[ind.i] <- NA
    }
    world.rot   <- geo2rot(pollon,pollat,world$x, world$y, polgam)

    if (rivers){
        riv.dat     <- map('rivers', plot=F)
        rivers.rot  <- geo2rot(pollon,pollat, riv.dat$x, riv.dat$y, polgam)
    } 

    data.tmp	    <- data
    data.tmp[!lsm]  <- NA
    image(rlon, rlat, data.tmp, breaks=levs, add=add,
        col=colours, axes=F, xlab=xlab, ylab=ylab, main=main,
        cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main)

    if (any(!lsm) & !is.null(sea.col)){
        ## sea points with sea.col
        data.tmp	<- data
        data.tmp[lsm]   <- NA
        image(rlon, rlat, data.tmp, breaks=c(-1e10,1e10), 
            col=sea.col, add=T, axes=F, xlab="", ylab="")
    }    

    if (exists("alt")){
        if (is.null(alt.lev)) alt.lev <- pretty(alt, 10)
        contour(rlon, rlat, alt, lev=alt.lev, drawlabels=F, add=T)
    }

    if (!is.null(file.small)){

        # make sponge zone transparent if dev == pdf, otherwise white
        if (!is.null(names(dev.cur())) & names(dev.cur()) == "pdf"){
            rect(min(rlon.small), min(rlat.small), max(rlon.small), max(rlat.small), 
                border=1, lwd=1, col=rgb(1,1,1,0.5))
        } else {
            rect(min(rlon.small), min(rlat.small), max(rlon.small), max(rlat.small), 
                border=1, lwd=1, col="white")
        }    

        data.small.tmp	<- data.small

        data.small.tmp[!lsm.small] <- NA
        image(rlon.small, rlat.small, data.small.tmp, breaks=levs, 
            col=colours,  add=T, axes=F, xlab="", ylab="")

        if (any(!lsm.small) & !is.null(sea.col)){
            data.small.tmp	        <- data.small
            data.small.tmp[lsm.small]   <- NA
            image(rlon.small, rlat.small, data.small.tmp, breaks=c(-1e10,1e10), 
                col=sea.col, add=T, axes=F, xlab="", ylab="")
        }

        if (exists("alt.small")){
        if (is.null(alt.lev)) alt.lev <- pretty(alt.small, 10)
            contour(rlon.small, rlat.small, alt.small, lev=alt.lev, drawlabels=F, add=T)
        }
    }

    if (rivers){
        riv.col <- sea.col
        if (is.null(riv.col)) riv.col <- "white"
        lines(rivers.rot$x, rivers.rot$y, col=riv.col)
    }

    if (!exists("alt")){
        ##lon-lat-lines
        lines(world.rot$x, world.rot$y, lwd=map.lwd)
    }

    if (cities){
        # some cities
        data(world.cities)
        # compute rotated coordinates
        coords <- geo2rot(pollon,pollat,world.cities$long,world.cities$lat,polgam)
        world.cities$rlon <- coords$x
        world.cities$rlat <- coords$y
        
        # select the cities within the plot area
        region.cities   <- world.cities[world.cities$rlon > min(rlon) &
            world.cities$rlat > min(rlat) &
            world.cities$rlon < max(rlon) &
            world.cities$rlat < max(rlat),]
        
        # further select the cities according to minpop
        # or ncities
        if (is.null(minpop)){
            region.cities   <- region.cities[order(region.cities$pop, 
                decreasing=TRUE),]
            region.cities   <- region.cities[1:ncities,]
        } else {
            region.cities   <- region.cities[region.cities$pop > minpop,]            
        }
        # their point 
        points(region.cities$rlon,region.cities$rlat, pch=city.pch)
        # their label
        if (label){
            text(region.cities$rlon, region.cities$rlat,
                labels=region.cities$name, offset=0.5, pos=3, cex=cex.txt)
        }
    }

    
    box(lwd=1)
    
    # check whether lon and lat are present, else set grid to FALSE
    if (!(exists("lon") & exists("lat"))){
        grid    <- FALSE
    }
    
    if (grid){

        if (missing(lon.ind)){
            lon.ind	<- pretty(lon,nlongrid)
        }
        if (missing(lat.ind)){
            lat.ind <- pretty(lat,nlatgrid)
        }

        contour(rlon, rlat, lon, levels=lon.ind, 
            lty=grid.lty, drawlabels=F, axes=F, add=T)
        contour(rlon, rlat, lat, levels=lat.ind, 
            lty=grid.lty, drawlabels=F, axes=F, add=T)

        if (grid.txt){
            lon.ind2 <- lon.ind
            if (any(lon > 180)) lon.ind2[lon.ind > 180] <- lon.ind[lon.ind > 180] - 360
            lon.txt <- paste(abs(lon.ind2), '*degree', c('~ W', '', '~ E')[sign(lon.ind2) + 2])
            lab.w <- strwidth(parse(text=lon.txt), cex=cex.axis)
            lon.i <- apply(as.matrix(lon.ind), 1, function(x) if (x > min(lon[,1]) & x < max(lon[,1])) which.min((lon[,1]-x)**2) else NA)
            lon.at <- geo2rot(pollon, pollat, lon[lon.i,1], lat[lon.i,1], polgam)$x 
            for (i in (min(which(!is.na(lon.at)))+1):max(which(!is.na(lon.at)))){
                lo.i <- max(which(!is.na(lon.at[1:(i-1)])))
                dist <- lon.at[i] - lon.at[lo.i]
                if (dist < 0.6*(lab.w[i] + lab.w[lo.i])) lon.at[i] <- NA
            }
            axis(1, at=lon.at, labels=parse(text=lon.txt), tick=F, line=-0.5, cex.axis=cex.axis)
            lon.i <- apply(as.matrix(lon.ind), 1, function(x) if (x > min(lon[,ncol(lon)]) & x < max(lon[,ncol(lon)])) which.min((lon[,ncol(lon)]-x)**2) else NA)
            lon.at <- geo2rot(pollon, pollat, lon[lon.i,ncol(lon)], lat[lon.i,ncol(lon)], polgam)$x 
            for (i in (min(which(!is.na(lon.at)))+1):max(which(!is.na(lon.at)))){
                lo.i <- max(which(!is.na(lon.at[1:(i-1)])))
                dist <- lon.at[i] - lon.at[lo.i]
                if (dist < 0.6*(lab.w[i] + lab.w[lo.i])) lon.at[i] <- NA
            }
            axis(3, at=lon.at, labels=parse(text=lon.txt), tick=F, line=-0.5, cex.axis=cex.axis)
            
            lat.txt <- paste(abs(lat.ind), '*degree', c('~ S', '', '~ N')[sign(lat.ind) + 2])
            lab.w <- strheight(parse(text=lat.txt), cex=cex.axis)
            lat.i <- apply(as.matrix(lat.ind), 1, function(x) if (x > min(lat[,1]) & x < max(lat[1,])) which.min((lat[1,]-x)**2) else NA)
            lat.at <- geo2rot(pollon, pollat, lon[1,lat.i], lat[1,lat.i], polgam)$y
            for (i in (min(which(!is.na(lat.at)))+1):max(which(!is.na(lat.at)))){
                lo.i <- max(which(!is.na(lat.at[1:(i-1)])))
                dist <- lat.at[i] - lat.at[lo.i]
                if (dist < (lab.w[i] + lab.w[lo.i])) lat.at[i] <- NA
            }
            axis(2, at=lat.at, labels=parse(text=lat.txt), tick=F, line=-0.5, cex.axis=cex.axis, las=1)
            lat.i <- apply(as.matrix(lat.ind), 1, function(x) if (x > min(lat[nrow(lat),]) & x < max(lat[nrow(lat),])) which.min((lat[nrow(lat),]-x)**2) else NA)
            lat.at <- geo2rot(pollon, pollat, lon[nrow(lat),lat.i], lat[nrow(lat),lat.i], polgam)$y
            for (i in (min(which(!is.na(lat.at)))+1):max(which(!is.na(lat.at)))){
                lo.i <- max(which(!is.na(lat.at[1:(i-1)])))
                dist <- lat.at[i] - lat.at[lo.i]
                if (dist < (lab.w[i] + lab.w[lo.i])) lat.at[i] <- NA
            }
            axis(4, at=lat.at, labels=parse(text=lat.txt), tick=F, line=-0.5, cex.axis=cex.axis, las=1)
        }
    }


    ## pollon, pollat
    out         <- list(pollon=pollon, pollat=pollat, polgam=polgam, 
        col=colours, lev=levs, sea.col=sea.col, flag_values=flag_values,
        flag_meanings=flag_meanings, longname=longname, units=units)
    class(out)  <- "plotmap"
    invisible(out)

}

