\name{frac_in_polygon}
\alias{frac_in_polygon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to compute what fraction of grid cells lies within  the given polygon}
\description{
  This function computes for a regular grid (one vector of longitudes and latitudes), what fraction of each grid cell lies within the given polygon.
}
\usage{
frac_in_polygon(lon, lat, poly.x, poly.y, multiply = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lon}{A vector of grid longitudes}
  \item{lat}{A vector of grid latitudes}
  \item{poly.x}{The longitudes of the edges of the polygon}
  \item{poly.y}{The latitudes of the edges of the polygon}
  \item{multiply}{Factor by which the grid is inflated for computation}
}
\details{
  The accuracy for this method scales with the multiplying factor. After inflation, grid points within the polygon are found with \code{point.in.polygon} from the package \code{sp}.
}
\value{
  Returns a matrix of dimension nlon x nlat of fractional gridbox area within the polygon.
}
\author{Jonas Bhend <email: jonas.bhend@env.ethz.ch>}
\seealso{\code{\link[sp]{point.in.polygon}}}
\examples{
  # set up latitudes and longitudes
  lon       <- seq(-3.5,3.5,0.5)
  lat       <- seq(-3.5,3.5,0.5)
  
  # polygon is a circle with radius 3
  poly.x    <- 3*sin(seq(0,2*pi,length=100))+0.13
  poly.y    <- 3*cos(seq(0,2*pi,length=100))-0.32
  
  outgrid   <- frac_in_polygon(lon,lat,poly.x,poly.y, multiply=5)
  
  # plot the result
  image(lon,lat,outgrid, asp=1)
  lines(poly.x, poly.y, lwd=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }
