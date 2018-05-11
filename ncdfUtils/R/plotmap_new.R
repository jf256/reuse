read_netcdf <- function(file, varname=NULL, lsm=TRUE, ...){
    # open the netcdf file
    nc <- open.ncdf(file)
    varnames <- names(nc$var)

    # get the dimensions
    rotated <- any(names(nc$dim) %in% c('rlon', 'rlat')) | any(names(nc$var) == 'rotated_pole')
    dimnames <- names(nc$dim)
    x <- nc$dim[[min(grep('lon', dimnames, ignore.case=T))]]$vals
    y <- nc$dim[[min(grep('lat', dimnames, ignore.case=T))]]$vals
    if (rotated){
        lon <- get.var.ncdf(nc, varnames[min(grep('lon', varnames, ignore.case=T))])
        lat <- get.var.ncdf(nc, varnames[min(grep('lat', varnames, ignore.case=T))])
    } else {
        lon <- array(rep(x, length(y)), c(length(x), length(y)))
        lat <- array(rep(y, each=length(x)), c(length(x), length(y)))
    }
    nlon <- length(x)
    nlat <- length(y)
        
    # if missing variable select the first in the netcdf file
    if (is.null(varname)){
        noread  <- c("lon", "lat", "lev", "rot")
        var.no <- NULL
        for (nr in noread) var.no <- c(var.no,grep(nr, varnames, ignore.case=TRUE))
        var.n <- varnames[setdiff(seq(along=varnames), var.no)]
        if (any(var.n == "HSURF")){
            varname <- "HSURF"
        } else {
            dims    <- lapply(nc$var[names(nc$var) %in% var.n], function(x) x$size)
            which.v <- sapply(dims, function(x) all(x[1:2] == c(nlon, nlat)))
            varname <- names(which.v)[which.max(which.v)]
        }
    }
    
    # read in projection data
    if (rotated){
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
    } else {
        pollon <- 0
        pollat <- 90
        polgam <- 0
    }
    
    .rotated_grid <<- list(pollon=pollon, pollat=pollat, polgam=polgam)
        
    # get the data
    data        <- get.var.ncdf(nc, varname)
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
    
    # get the land-sea mask if needed
    if (is.character(lsm)) {
        if (file.exists(lsm)) {
            lsm.nc <- open.ncdf(lsm) 
            lsm.grep <- pmatch(tolower(names(lsm.nc$var)), c('lsm', 'slm', 'fr_land'))
            if (any(!is.na(lsm.grep))) {
                lsm.name <- names(lsm.nc$var)[min(which(!is.na(lsm.grep)), na.rm=T)]
                lsm <- get.var.ncdf(lsm.nc, lsm.name)
            }
            close.ncdf(lsm.nc)
        }
    } else if (lsm) {
        lsm.grep <- pmatch(tolower(varnames), c('lsm', 'slm', 'fr_land'))
        if (any(!is.na(lsm.grep))){
            lsm.name <- varnames[min(which(!is.na(lsm.grep)), na.rm=T)]
            lsm <- get.var.ncdf(nc, lsm.name)
        }
    } 

    if (is.matrix(lsm)) {
        lsm <- lsm*1 > mean(range(lsm))
    } else {
        lsm <- NULL
    }
    
    if (length(lsm) > prod(dim(data)[1:2])) lsm <- array(lsm, dim(data)[1:2])
    
    if (all(diff(y) < 0)) {
        y <- rev(y)
        lon <- lon[,ncol(lon):1]
        lat <- lat[,ncol(lat):1]
        data <- data[,ncol(data):1]
        lsm <- lsm[,ncol(lsm):1]
    }
    
    # write out the data in a list
    out <- list(data=data, x=x, y=y, lon=lon, lat=lat, nlon=nlon, nlat=nlat,
        rotated=rotated, pollon=pollon, pollat=pollat, polgam=polgam,
        varname=varname, longname=longname, units=units, flag_values=flag_values,
        flag_meanings=flag_meanings, lsm=lsm)
    lev.exists <- grep('lev', names(nc$dim)[nc$var[[varname]]$dimis], ignore.case=T)
    if (length(lev.exists) > 0){
        z.i <- grep('lev', dimnames, ignore.case=T)
        if (length(z.i) > 0){
            z <- nc$dim[[min(z.i)]]$vals
            out[['z']] <- z
        }
        lev.i <- grep('lev', varnames, ignore.case=T)
        if (length(lev.i) > 0){
            lev <- get.var.ncdf(nc, varnames[min(lev.i)])
            out[['lev']] <- lev
        }
    }

    close.ncdf(nc)

    invisible(out)
}

# --------------------------------------------------------
dataproject <- function(inlist, projection=NULL, parameters=NULL, orientation=NULL){
    require(mapproj)
    attach(inlist, pos=2)
    on.exit(detach(inlist, pos=2))
    # plots everything in native model grid unless projection is specified
    if (is.null(projection)){
        if (!(pollon == .rotated_grid$pollon &
            pollat == .rotated_grid$pollat &
            polgam == .rotated_grid$polgam)){
                xy <- geo2rot(.rotated_grid$pollon, rotated_grid$pollat, lon, lat, rotated_grid$polgam)             
                inlist$x <- xy$x
                inlist$y <- xy$y
            }
    } else {
        xy <- mapproject(lon, lat, projection=projection, parameters=parameters, orientation=orientation)
        inlist$x <- xy$x
        inlist$y <- xy$y
    }
    inlist$projection <- projection
    inlist$parameters <- parameters
    inlist$orientation <- orientation
    invisible(inlist)
}

# --------------------------------------------------------
plot_poly <- function(inlist, ti=1, sponge=8, lev=NULL, symmetric=FALSE, levtype='normal',
    col=NULL, sea.col='lightblue', colourplot=TRUE, b.lwd=2, b.col=1,
    axes=TRUE, xlab='', ylab='', lons=NULL, lats=NULL,
    grid=TRUE, grid.lty=3, grid.lwd=1, tick=F, add=FALSE, contour=FALSE, border=FALSE, transparency=TRUE, ...) {
    attach(inlist, pos=2)
    on.exit(detach(inlist, pos=2))
    if (length(dim(data)) == 3) data <- data[,,ti]
    lsm <- inlist$lsm
    levtypes <- c('normal', 'log', 'alt')
    if (is.null(lev)){
        if (is.null(col)) {
            nlev <- 12/(symmetric + 1)
        } else {
            nlev <- (length(col) + 1)/(symmetric + 1)
        }
        levdata <- data
        if (symmetric) levdata <- abs(data)
        if (pmatch(levtype, levtypes) == 1) {
            lev <- pretty(levdata, nlev)
        } else if (pmatch(levtype, levtypes) == 2) {
            if (!symmetric & diff(sign(range(data, na.rm=T))) == 2) frac <- 2 else frac <- 1
            lev <- exp(pretty(log(abs(levdata)), nlev/frac))
            lev2 <- c(0,lev)
            for (i in seq(along=lev)) lev[i] <- round(lev[i], 1-round(log(diff(lev2[i+0:1]))/log(10)))
            if (!symmetric) {
                lev <- c(rev(-lev[2:length(lev)]), lev[2:length(lev)])
                lev <- lev[(min(which(lev > min(data, na.rm=T)))-1):(max(which(lev < max(data, na.rm=T)))+1)]
            }
        } else if (pmatch(levtype, levtypes) == 3){
            lev <- c(0,50,100,200,300,500,700,1000,1500,2000,3000,5000,10000)
        }
        if (symmetric) lev <- sort(unique(c(-lev, 0, lev)))
    }

    ncols <- length(lev)-1
    if (is.null(col)){
        if (colourplot){
            if (varname == "HSURF" | pmatch(levtype, levtypes) == 3){
                colours	<- .colseq(length(lev)-1, .hsurf, smooth=0)
                sea.col <- .water
            } else if (varname == "SOILTYP"){
                colours <- .soil[flag_values+1]
            } else if (varname %in% c("TOT_PREC", "precip", "pr")){
                colours <- .colseq(length(lev)-1, .gpcc, smooth=0)
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
    
    if (sponge > 0){
        s.ind <- array(seq(along=data), dim(data))[(sponge+1):(nrow(data)-sponge), (sponge+1):(ncol(data)-sponge)]
        s.data <- data
        s.data[s.ind] <- NA
        s.data[-s.ind] <- 1
        if (transparency){
            s.col <- rgb(1,1,1,0.2)
        } else {
            s.col <- rgb(1,1,1)
        }
        s.lev <- c(0,2)
    }
    
    if (length(x) == length(y) & length(x) == length(data)){
        print('using poly_image')
        poly_image(x, y, data, lsm, levels=lev, col=colours, sea.col=sea.col, xlab=xlab, ylab=ylab, axes=F, add=add, b.col=b.col, b.lwd=b.lwd, ...)
        if (sponge > 0){
            poly_image(x,y,s.data, levels=s.lev, col=s.col, xlab='', ylab='', axes=F, add=TRUE, b.lwd=b.lwd, b.col=b.col)
        }
    } else if (is.matrix(data) & length(x) == nrow(data) & length(y) == ncol(data)){
        if (contour){
            if (!is.null(lsm)){
                data.tmp <- data
                data.tmp[!lsm] <- NA
                shaded_contour(x,y,data.tmp, lev=lev, col=colours, na.col=sea.col, xlab='', ylab='', axes=F, border=border, add=add, ...)
            } else {
                shaded_contour(x, y, data, lev=lev, col=colours, xlab=xlab, ylab=ylab, axes=F, border=border, add=add, ...)        
            }
        } else {
            if (!is.null(lsm)){
                data.tmp <- data
                data.tmp[lsm] <- NA
                image(x,y,data.tmp,breaks=c(-1,1)*1e10, col=sea.col, xlab='', ylab='', axes=F, add=add, ...)
                data.tmp <- data
                data.tmp[!lsm] <- NA
                image(x,y,data.tmp, breaks=lev, col=colours, xlab=xlab, ylab=ylab, axes=F, add=TRUE, ...)
            } else {
                image(x, y, data, breaks=lev, col=colours, xlab=xlab, ylab=ylab, axes=F, add=add, ...)        
            }
        }
        if (sponge > 0){
            image(x,y,s.data, breaks=s.lev, col=s.col, axes=F, add=TRUE)
        }
    }
    if (!add) box()

    out         <- list(pollon=pollon, pollat=pollat, polgam=polgam, 
        col=colours, lev=lev, sea.col=sea.col, flag_values=flag_values,
        flag_meanings=flag_meanings, longname=longname, units=units, varname=varname)
    class(out)  <- "plotmap"

    invisible(out)
}

# --------------------------------------------------------
add_map <- function(inlist, hires=FALSE, interior=F, alt.poli=TRUE, 
    rivers=TRUE, riv.lwd=1, riv.col='lightblue',
    cities=FALSE, minpop=NULL, ncities=100, city.pch=19, 
    label=FALSE, cex.txt=0.9, ...) {
    attach(inlist, pos=2)
    on.exit(detach(inlist, pos=2))
    if (hires){
        require(mapdata)
        worlddb <- "worldHires"
        if (any(lon > 180)) worlddb <- "world2Hires"
    } else {
        worlddb <- "world"
        if (any(lon > 180)) worlddb <- "world2"
    }
    if (alt.poli & interior){
      data(polibound)
      world <- polibound
      if (any(lon > 180)){
          world$x[!is.na(world$x) & world$x < 0] <- world$x[!is.na(world$x) & world$x < 0] + 360
          world$x[abs(diff(world$x)) > 180] <- NA
      }
    } else {
        world       <- map(worlddb,interior=interior, plot=F, 
                         xlim=range(lon) + c(-0.1,0.1)*diff(sort(range(lon))),
                         ylim=range(lat) + c(-0.1, 0.1)*diff(sort(range(lat))))
    }
    for (add.name in c(".*Lake.*", ".*Sea.*")){
        world.add <- try(map(worlddb, region=add.name, plot=F, xlim=range(lon), ylim=range(lat)), silent=TRUE)
        if (class(world.add) != 'try-error' & length(world.add) > 0){
            world   <- list(x=c(world$x, NA, world.add$x),
                            y=c(world$y, NA, world.add$y))
      }
    }
    if (rivers) {
        riv.dat     <- map('rivers', plot=F)   
        if (any(lon > 180)) riv.dat$x[!is.na(riv.dat$x) & riv.dat$x < 0] <- riv.dat$x[!is.na(riv.dat$x) & riv.dat$x < 0] + 360 
    }
    if (cities) {
        data(world.cities)
        if (any(lon > 180)) world.cities$long[world.cities$long < 0] <- world.cities$long[world.cities$long < 0] + 360
    }
    
    if (!is.null(inlist$proj)) {
        world <- mapproject(world$x, world$y, proj=projection, par=parameters, ori=orientation) 
        if (rivers) riv.dat <- mapproject(riv.dat$x, riv.dat$y, proj=projection, par=parameters, ori=orientation)
        if (cities) {
            city.dat <- mapproject(world.cities$long, world.cities$lat, proj=projection, par=parameters, ori=orientation)
            world.cities$long <- city.dat$x
            world.cities$lat <- city.dat$y
        }
    } else if (rotated) {
        world <- geo2rot(pollon, pollat, world$x, world$y, polgam=polgam)
        if (rivers) riv.dat <- geo2rot(pollon, pollat, riv.dat$x, riv.dat$y, polgam=polgam)
        if (cities) {
            city.dat <- geo2rot(pollon, pollat, world.cities$long, world.cities$lat, polgam=golgam)
            world.cities$long <- city.dat$x
            world.cities$lat <- city.dat$y
        }
    }
    lines(world, ...)
    if (rivers) lines(riv.dat, lwd=riv.lwd, col=riv.col)
    if (cities){
        # select the cities within the plot area
        region.cities   <- world.cities[world.cities$long > min(x) &
            world.cities$lat > min(y) &
            world.cities$lon < max(x) &
            world.cities$lat < max(y),]
        
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
        points(region.cities$lon,region.cities$lat, pch=city.pch, cex=cex.txt)
        # their label
        if (label){
            text(region.cities$lon, region.cities$lat,
                labels=region.cities$name, offset=0.5, pos=3, cex=cex.txt)
        }
    }

}

# --------------------------------------------------------
add_grid <- function(inlist, axes=TRUE, lats=NULL, lons=NULL, grid.lty=3, grid.lwd=1) {
    attach(inlist, pos=2)
    on.exit(detach(inlist, pos=2))
    # find out if plot over date-line
    dlon <- apply(lon, 2, diff)
    if (any(abs(dlon) > 180)) lon[lon < 0] <- lon[lon < 0] + 360
    if (is.null(lons)) lons <- pretty(lon, 10)
    if (is.null(lats)) lats <- pretty(lat, 10)
    if (length(x) == nrow(lon) & length(y) == ncol(lon)){
        contour(x, y, lon, lev=lons, lty=grid.lty, lwd=grid.lwd, drawlabels=F, add=T)
        contour(x, y, lat, lev=lats, lty=grid.lty, lwd=grid.lwd, drawlabels=F, add=T)
    } else if (!is.null(projection)){
        dx <- median(diff(lons))
        dy <- median(diff(lats))
        map.grid(c(0,360, -90, 90), nx=360/dx + 1, ny=180/dy + 1, col=1, lty=grid.lty, lwd=grid.lwd, labels=FALSE)
    }
    if (axes) {
        # don't write labels in the first 5 percent of the axis range
        lon.range <- par("usr")[1:2] + c(0.05, - 0.05) * diff(par("usr")[1:2])
        lat.range <- par("usr")[3:4] + c(0.02, - 0.02) * diff(par("usr")[3:4])
        xx <- array(rep(x, length.out=length(lon)), dim(lon))
        lon1 <- spline(x=lon[,1], y=xx[,1], xout=lons)
        lon.i <- lon1$y > lon.range[1] & lon1$y < lon.range[2]
        lolalabels(lapply(lon1, function(x) x[lon.i]), side=1, lola='lon')
        lon1 <- spline(x=lon[,ncol(lon)], y=xx[, ncol(lon)], xout=lons)
        lon.i <- lon1$y > lon.range[1] & lon1$y < lon.range[2]
        lolalabels(lapply(lon1, function(x) x[lon.i]), side=3, lola='lon')
        yy <- t(array(rep(y, length.out=length(lat)), rev(dim(lat))))
        lat1 <- spline(x=lat[1,], y=yy[1,], xout=lats)
        lat.i <- lat1$y > lat.range[1] & lat1$y < lat.range[2]
        lolalabels(lapply(lat1, function(x) x[lat.i]), side=2, lola='lat')
        lat1 <- spline(x=lat[nrow(lat),], y=yy[nrow(lat),], xout=lats)
        lat.i <- lat1$y > lat.range[1] & lat1$y < lat.range[2]
        lolalabels(lapply(lat1, function(x) x[lat.i]), side=4, lola='lat')
    }
}

lolalabels <- function(loni, side=1, lola='lon', tick=F, cex.axis=par('cex.axis')){
    lon.ind <- loni$x
    lon.ind[lon.ind > 180] <- lon.ind[lon.ind > 180] - 360
    lon.idx <- c('~ W', '', '~ E')
    if (lola == 'lat') lon.idx <- c('~ S', '', '~ N')
    lon.txt <- paste(abs(lon.ind), '*degree', lon.idx[sign(lon.ind) + 2])
    lab.w <- strwidth(parse(text=lon.txt), cex=cex.axis)
    lon.i <- apply(as.matrix(lon.ind), 1, function(x) if (x > min(lon[,1]) & x < max(lon[,1])) which.min((lon[,1]-x)**2) else NA)
    lon.at <- loni$y
    for (i in (min(which(!is.na(lon.at)))+1):max(which(!is.na(lon.at)))){
        lo.i <- max(which(!is.na(lon.at[1:(i-1)])))
        dist <- lon.at[i] - lon.at[lo.i]
        if (dist < 0.6*(lab.w[i] + lab.w[lo.i])) lon.at[i] <- NA
    }
    axis(side, at=lon.at, labels=parse(text=lon.txt), tick=tick, line=-0.5, cex.axis=cex.axis)
}


# --------------------------------------------------------
edges   <- function(x){
    edges   <- c(x - 0.5*diff(x)[c(1,seq(diff(x)))], 
        x[length(x)] + 0.5*diff(x)[length(x) - 1])
    edges
}

# --------------------------------------------------------
edges_grid  <- function(lon,lat){
    if (!is.matrix(lon) | length(lon) != length(lat)){
        # assume it is the grid coordinates in 1-d
        nlon    <- length(lon)
        nlat    <- length(lat)
        lon <- array(lon, c(nlon,nlat))
        lat <- t(array(lat, c(nlat, nlon)))
    }
    elon    <- t(apply(apply(lon, 2, edges), 1, edges))                            
    elat    <- t(apply(apply(lat, 2, edges), 1, edges))                            
    # check if  -180 to 180 ...
    ediff   <- abs(elon[1,] - elon[nrow(elon),])
    ediff2  <- abs(elon[1,] - elon[2,])
    if (mean(ediff/ediff2, na.rm=T) < 0.2){
        eint    <- 0.5*(elon[1,] + elon[nrow(elon),])
        elon[1,]<- eint
        elon[nrow(elon),]<- eint
        eint    <- 0.5*(elat[1,] + elat[nrow(elat),])
        elat[1,]<- eint
        elat[nrow(elat),] <- eint
    }        
    egrid <- list()
    N   <- nrow(elon)
    M   <- ncol(elon)
    egrid$x <- rbind(as.vector(elon[1:(N-1), 1:(M-1)]), 
                    as.vector(elon[1:(N-1), 2:M]),
                    as.vector(elon[2:N, 2:M]), 
                    as.vector(elon[2:N, 1:(M-1)]),
                    NA)
    egrid$y <- rbind(as.vector(elat[1:(N-1), 1:(M-1)]), 
                    as.vector(elat[1:(N-1), 2:M]),
                    as.vector(elat[2:N, 2:M]), 
                    as.vector(elat[2:N, 1:(M-1)]),
                    NA)
    return(egrid)
}

# --------------------------------------------------------
poly_image <- function (x, y, z, lsm=NULL, levels = pretty(z, 20), 
    col = rbfun(length(levels) - 1), sea.col=NA,
    xlim = range(x, na.rm=T), ylim = range(y, na.rm=T), 
    add = FALSE, b.col=1, b.lwd=3, ...) 
{
    if (length(x) != length(y) | length(x) != length(z) | !is.matrix(z)){
        stop('Dimensions in poly_image do not match')
    }
    x <- array(x, dim(z))
    y <- array(y, dim(z))
    egrid <- edges_grid(x, y)
    # find out about colours
    zcol <- apply(z, 1:2, function(x) if (!is.na(x)) sum(levels < x) else NA)
    zcol[zcol == 0] <- NA
    zcol[zcol == length(levels)] <- NA
    cols <- c(col, sea.col)
    if (!is.null(lsm)) zcol[!lsm] <- length(cols)
    
    if (!add) {
        plot(xlim, ylim, type = "n", ...)
    }
    emask <- apply(!is.na(egrid$x), 2, sum) > 2
    polygon(egrid$x[,emask], egrid$y[,emask], border = cols[zcol[emask]], col = cols[zcol[emask]], lend=3)
    if (b.lwd > 0){
        enew <- lapply(egrid, function(x) array(x, c(5, dim(z))))
        border <- cbind(
            c(apply(enew$x[,1,], 2, min, na.rm=T),
              apply(enew$x[,,dim(enew$x)[3]], 2, mean, na.rm=T),
              rev(apply(enew$x[,dim(enew$x)[2],], 2, max, na.rm=T)),
              rev(apply(enew$x[,,1], 2, mean, na.rm=T))),
            c(apply(enew$y[,1,], 2, mean, na.rm=T),
              apply(enew$y[,,dim(enew$y)[3]], 2, max, na.rm=T),
              rev(apply(enew$y[,dim(enew$y)[2],], 2, mean, na.rm=T)),
              rev(apply(enew$y[,,1], 2, min, na.rm=T))))
        lines(border, col=b.col, lwd=b.lwd)
    }
}



# --------------------------------------------------------
`oldplotmap_rot` <-
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
    
    filelarge <- read_netcdf(file, varname)
    if (is.null(varname)) varname <- names(filelarge)[1]
    attach(filelarge)
    
    # select time and levels if appropriate
    if (length(dim(data)) == 3){
        data    <- data[,,i.time]
    } else if (length(dim(data)) == 4){
        data    <- data[,,i.lev,i.time]
    }

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

