.k2p <- function(v){
    if (length(v) == 3){
        v   <- t(t(v))
    }
    
    lon     <- atan(v[2,]/v[1,])/pi*180
    indi    <- v[1,] < 0 & v[2,] < 0 & !is.na(v[1,]) & !is.na(v[2,])
    lon[indi] <- lon[indi] - 180
    indi    <- v[1,] < 0 & v[2,] >= 0 & !is.na(v[1,]) & !is.na(v[2,])
    lon[indi] <- lon[indi] + 180
    
    lat     <- asin(v[3,])/pi*180
    if (length(lat) == 1){
        p   <- c(lon,lat)
    } else {
        p   <- cbind(lon,lat)
    }
        
    p
}

##########################################################

.p2k <- function(lon,lat){

    lonr <- lon/180*pi
    latr <- lat/180*pi

    if (length(lonr) == 1){
        v <- c(cos(lonr)*cos(latr),sin(lonr)*cos(latr),sin(latr)) 
    } else {
        v <- rbind(cos(lonr)*cos(latr),sin(lonr)*cos(latr),sin(latr))
    }
    v
}

##########################################################

.phirot2phi  <- function(phirot, rlarot, polphi, pollam, polgam=0){
    # Description:
    #   This function converts phi from one rotated system to phi in another
    #   system. If the optional argument polgam is present, the other system
    #   can also be a rotated one, where polgam is the angle between the two
    #   north poles.
    #   If polgam is not present, the other system is the real geographical
    #   system.
    #
    # Method:
    #   Transformation formulas for converting between these two systems.
    #
    
    # convert to -180 to 180
    ind         <- rlarot > 180 & !is.na(rlarot)
    rlarot[ind]    <- rlarot[ind] - 360
    
    # convert to -pi to pi
    polphi      <- polphi / 180 * pi
    pollam      <- pollam / 180 * pi
    polgam      <- polgam / 180 * pi
    phirot      <- phirot / 180 * pi
    rlarot      <- rlarot / 180 * pi
    
    zsinpol     <- sin(polphi)
    zcospol     <- cos(polphi)
 
    # adjust for switched north pole
    if (polgam != 0){
        zarg    <- zsinpol * sin(phirot) + zcospol * cos(phirot) *
            (cos(rlarot)*cos(polgam) - sin(polgam) * sin(rlarot))
    } else {
        zarg    <- zcospol * cos(phirot) * cos(rlarot) + zsinpol * sin(phirot)
    }
    
    out         <- asin(zarg) * 180 / pi
}


##########################################################

.phi2phirot  <- function(phi, rla, polphi, pollam){

    #------------------------------------------------------------------------------
    # Description:
    #   This routine converts phi from the real geographical system to phi
    #   in the rotated system.
    #
    # Method:
    #   Transformation formulas for converting between these two systems.
    #
    #------------------------------------------------------------------------------

    zsinpol     <- sin(polphi / 180 * pi)
    zcospol     <- cos(polphi / 180 * pi)
    zlampol     <- pollam / 180 * pi
    zphi        <- phi / 180 * pi
    
    ind         <- rla > 180 & !is.na(rla)
    rla[ind]    <- rla[ind] - 360
    zrla        <- rla / 180 * pi

    zarg1       <- sin(zphi) * zsinpol
    zarg2       <- cos(zphi) * zcospol * cos(zrla - zlampol)

    phi2phirot  <- asin(zarg1 + zarg2) * 180 / pi
    
    phi2phirot

}

##########################################################

.rlarot2rla  <- function(phirot, rlarot, polphi, pollam, polgam=0){

    # Description:
    #   This function converts lambda from one rotated system to lambda in another
    #   system. If the optional argument polgam is present, the other system
    #   can also be a rotated one, where polgam is the angle between the two
    #   north poles.
    #   If polgam is not present, the other system is the real geographical
    #   system.
    #
    # Method:
    #   Transformation formulas for converting between these two systems.
    #
    
    zpir18      <- pi/180
    zrpi18      <- 180/pi
    
    zsinpol     <- sin(zpir18 * polphi)
    zcospol     <- cos(zpir18 * polphi)

    zlampol     <- zpir18 * pollam
    zphis       <- zpir18 * phirot
    
    # convert to -180 to 180
    ind         <- rlarot > 180 & !is.na(rlarot)
    rlarot[ind] <- rlarot[ind] - 360
    zrlas       <- zpir18 * rlarot


    if (polgam != 0) {
        zgam    <- zpir18 * polgam
        zarg1   <- sin(zlampol) * (- zsinpol*cos(zphis) * (cos(zrlas)*cos(zgam) - 
            sin(zrlas)*sin(zgam)) + zcospol * sin(zphis)) - cos(zlampol)*cos(zphis) * 
            (sin(zrlas)*cos(zgam) + cos(zrlas)*sin(zgam))
        zarg2   <- cos(zlampol) * (- zsinpol*cos(zphis) * (cos(zrlas)*cos(zgam) - 
            sin(zrlas)*sin(zgam)) + zcospol * sin(zphis)) + sin(zlampol)*cos(zphis) *
            (sin(zrlas)*cos(zgam) + cos(zrlas)*sin(zgam))
    } else {
        zarg1   <- sin(zlampol) * (-zsinpol * cos(zrlas) * cos(zphis) + zcospol * 
            sin(zphis)) - cos(zlampol) * sin(zrlas) * cos(zphis)
        zarg2   <- cos(zlampol) * (-zsinpol * cos(zrlas) * cos(zphis) + zcospol *
            sin(zphis)) + sin(zlampol) * sin(zrlas) * cos(zphis)
    }
    
    zarg2[zarg2 == 0] <- 1.0e-20

    rlarot2rla  <- zrpi18 * atan2(zarg1, zarg2)
    rlarot2rla
}

##########################################################

.rla2rlarot  <- function(phi, rla, polphi, pollam, polgam=0) { 

    #------------------------------------------------------------------------------
    #
    # Description:
    #   This routine converts lambda from the real geographical system to lambda 
    #   in the rotated system.
    #
    # Method:
    #   Transformation formulas for converting between these two systems.
    #
    #------------------------------------------------------------------------------

    zpir18      <- pi/180
    zrpi18      <- 180/pi
    
    zsinpol     <- sin(zpir18 * polphi)
    zcospol     <- cos(zpir18 * polphi)
    zlampol     <- zpir18 * pollam
    zphi        <- zpir18 * phi

    # convert everything to -180 to 180
    ind         <- rla > 180 & !is.na(rla)
    rla[ind]    <- rla[ind] - 360
    zrla        <- zpir18 * rla

    zarg1       <- - sin(zrla-zlampol) * cos(zphi)
    zarg2       <- - zsinpol * cos(zphi) * cos(zrla-zlampol) + 
        zcospol * sin(zphi)

    zarg2[zarg2 == 0] <- 1.0e-20

    rla2rlarot  <- zrpi18 * atan2(zarg1,zarg2)

    if (polgam != 0) {
        rla2rlarot   <- polgam + rla2rlarot
    }
    ind         <- rla2rlarot > 180 & !is.na(rla2rlarot)
    rla2rlarot[ind] <- rla2rlarot[ind] - 360
    
    rla2rlarot
}


##########################################################

geo2rot <- function(pollon, pollat, lon, lat, polgam=0){

    rlon    <- .rla2rlarot(lat, lon, pollat, pollon, polgam)
    rlat    <- .phi2phirot(lat, lon, pollat, pollon)
    
    out     <- list(x=rlon, y=rlat)
    out
}

##########################################################

rot2geo <- function(pollon, pollat, rlon, rlat, polgam=0){

    lon     <- .rlarot2rla(rlat, rlon, pollat, pollon, polgam)
    lat     <- .phirot2phi(rlat, rlon, pollat, pollon, polgam)
    
    out     <- list(x=lon, y=lat)
    out
}










        
    
    
