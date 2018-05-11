\name{ncdfUtils-package}
\alias{ncdfUtils-package}
\alias{ncdfUtils}
\docType{package}
\title{
Package for plotting and analysing NetCDF data
}
\description{
  This package provides several functions to plot NetCDF data on rotated or regular lonlat grids, as well as a function to read in the time dimension of NetCDF files and a function to compute grid point fractions within a given polygon.
}
\details{
\tabular{ll}{
Package: \tab ncdfUtils\cr
Type: \tab Package\cr
Version: \tab 0.4-9\cr
Date: \tab 2010-05-07\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
  The time dimesion of NetCDF files can be read with:
  
  \code{nc <- open.ncdf("<path-to-netcdf-file>")}
  
  \code{times <- ncdf_times(nc)}
  
  Data on rotated grids can be plotted with:
  
  \code{plotmap_rot("<path-to-netcdf-file>")}
  
  Fractional area within a given polygon can be computed via:
  
  \code{frac_in_polygon(lon, lat, poly.lon, poly.lat)}
  
  Plots within already existing plots are setup via:
  
  \code{set_position(x,y,size=0.2)}
  
}
\author{
  Maintainer: Jonas Bhend <jonas.bhend@env.ethz.ch>
}
\keyword{ package }
\seealso{
  \code{\link[<ncdf>:<ncdf>-package]{<ncdf>}}, 
  \code{\link[<sp>:<sp>-package]{<sp>}}, 
  \code{\link[<maps>:<maps>-package]{<maps>}}, 
  \code{\link[<mapdata>:<mapdata>-package]{<mapdata>}}, 
  \code{\link[<udunits>:<udunits>-package]{<udunits>}} 
}
