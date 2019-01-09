`frac_in_polygon` <- 
function(lon, lat, poly.x, poly.y, multiply=4){

    # blow up longitudes
    nlon    <- length(lon)
    lon.alt <- c(lon[1] - 0.5*diff(lon[1:2]), lon, 
        lon[nlon] + 0.5*diff(lon[nlon - 1:0]))
    lon.ind <- c(0,seq(0.5/nlon, 1-0.5/nlon,length=nlon),1)
    ind.out <- seq(0.5/nlon/multiply, 1- 0.5/nlon/multiply, length=nlon*multiply)
    lon.new <- approx(lon.ind, lon.alt, ind.out)$y
    
    # blow up latitudes
    nlat    <- length(lat)
    lat.alt <- c(lat[1] - 0.5*diff(lat[1:2]), lat, 
        lat[nlat] + 0.5*diff(lat[nlat - 1:0]))
    lat.ind <- c(0,seq(0.5/nlat, 1-0.5/nlat,length=nlat),1)
    ind.out <- seq(0.5/nlat/multiply, 1- 0.5/nlat/multiply, length=nlat*multiply)
    lat.new <- approx(lat.ind, lat.alt, ind.out)$y
    
    # find points in polygon
    outgrid <- point.in.polygon(rep(lon.new, length(lat.new)),
        rep(lat.new, each=length(lon.new)), poly.x, poly.y)
    if (multiply == 1){
        outgrid[outgrid > 1] <- NA
    } else {
        outgrid[outgrid > 1] <- 0.5
    }
    
    out.tmp <- array(outgrid, c(multiply, nlon, multiply, nlat))
    fracout <- apply(out.tmp, c(2,4), mean, na.rm=T)
    
    fracout
}
