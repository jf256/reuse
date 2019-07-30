\name{plotmap_rot}
\alias{plotmap_rot}
\alias{plotmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot NetCDF data on rotated or regular lonlat grids}
\description{
  This function reads in NetCDF data and produces a lon-lat plot of the specified variable.
}
\usage{
plotmap_rot(file, file.small = NULL, sponge = 8, 
    sponge.small = 15, varname = NULL, 
    lsm.file = NULL, lsm.file.small = NULL, col = NULL, 
    levels = NULL, sea.col = NULL, rivers = T, cities = T, 
    label = TRUE, minpop = NULL, ncities = 10, city.pch = 19, 
    alt.contour = F, alt.lev = NULL, grid = TRUE, 
    grid.txt = TRUE, grid.lty = 2, i.time = 1, i.lev = 1, 
    map.lwd = 2, cex.axis = 1, cex.lab = 1, cex.main = 1, 
    cex.txt = 1, main = "", xlab = "", ylab = "", 
    colourplot = TRUE, hires = FALSE, interior = FALSE,
    alt.poli=TRUE, add=TRUE, nlongrid=10, nlatgrid=5,
    lon.ind, lat.ind)
    
plotmap(file, file.small = NULL, sponge = 8, 
    sponge.small = 15, varname = NULL, lsm.file = NULL, 
    lsm.file.small = NULL, col = NULL, levels = NULL, 
    sea.col = NULL, alt.contour = F, alt.lev = NULL, 
    grid = TRUE, grid.txt = TRUE, grid.lty = 2, i.time = 1, 
    i.lev = 1, map.lwd = 2, cex.axis = 1, cex.lab = 1, 
    cex.main = 1, main = "", xlab = "", ylab = "", 
    colourplot = TRUE, hires = FALSE, interior = FALSE,
    alt.poli=TRUE, add=TRUE, nlongrid=10, nlatgrid=5,
    lon.ind, lat.ind)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{Path to the NetCDF file containing the data to plot}
  \item{file.small}{... for the nest in the large region}
  \item{sponge}{Number of grid cells to mask out at boundary}
  \item{sponge.small}{... for the nest in the large region}
  \item{varname}{Variable to plot (defaults to \code{HSURF} for orography or any variable with at least lon and lat dimensions)}
  \item{lsm.file}{Additional file containing land-sea mask and/or altitude information in \code{FR\_LAND} and \code{HSURF}}
  \item{lsm.file.small}{... for the nest in the large region}
  \item{col}{Colours to use}
  \item{levels}{Levels to use}
  \item{sea.col}{Colour for masked out water grid points}
  \item{rivers}{Switch for automatically added rivers}
  \item{cities}{Switch for automatically added cities}
  \item{label}{Switch for labels of cities}
  \item{minpop}{Minimal population for cities to get drawn}
  \item{ncities}{Number of cities to draw}
  \item{city.pch}{Symbol for city points}
  \item{cex.txt}{Character expansion for city labels}
  \item{alt.contour}{Switch for usage of the coastline as in file or lsm.file}
  \item{alt.lev}{Level at which to draw coastline (if alt.contour == TRUE) in units of "\code{FR\_LAND}"}
  \item{grid}{Switch for inclusion of grid lines}
  \item{grid.txt}{Switch for inclusion of grid axis labels at the border of the plot}
  \item{grid.lty}{Line type for grid lines}
  \item{i.time}{Time step to plot}
  \item{i.lev}{Level to plot}
  \item{map.lwd}{Line width for the coastline/ country border}
  \item{hires}{Use hi-resolution coastline}
  \item{interior}{Draw political boundaries}
  \item{colourplot}{Switch for colour or grayscale plot}
  \item{alt.poli}{Switch to use alternative political boundaries}
  \item{add}{Should plot be added to existing plot}
  \item{nlongrid, nlatgrid}{Approximate number of longitude and latitude grid lines}
  \item{lon.ind, lat.ind}{Sequence of longitudes (latitudes) at which to draw grid lines. Has to be increasing. Longitudes greater than 180 are converted.}
}
\value{
  A List of class "plotmap" containing the following objects:
  \item{pollon}{Longitude of the rotated pole}
  \item{pollat}{Latitude of the rotated pole}
  \item{polgam}{Switch for mirrored rotation (e.g. 0 or -180 degrees)}
  \item{levs}{The levels used}
  \item{col}{The colours used}
  \item{sea.col}{The colour used for water grid points}
  \item{longname}{The longname attribute of the variable plotted}
  \item{units}{The units of the variable plotted}
}
\author{Jonas Bhend: jonas.bhend@env.ethz.ch}
\seealso{\code{\link{plot_colourbar.plotmap}}}
\examples{
  print("Sorry no examples, too complex...")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aplot }
