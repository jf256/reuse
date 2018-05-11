\name{geo2rot}
\alias{geo2rot}
\alias{rot2geo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to convert }
\description{
  Functions to convert rotated latitudes and longitudes to positions on regular grids and vice versa.
}
\usage{
  geo2rot(pollon, pollat, lon, lat, polgam = 0)

  rot2geo(pollon, pollat, rlon, rlat, polgam = 0)
}

\arguments{
  \item{pollon}{Longitude of the rotated pole in degrees}
  \item{pollat}{Latitude of the rotated pole in degrees}
  \item{polgam}{Switch for mirrored rotation (0 or -180 deg.)}
  \item{lon}{A vector of longitudes}
  \item{lat}{A vector of latitudes}
  \item{rlon}{A vector of longitudes on the rotated grid}
  \item{rlat}{A vector of latitudes on the rotated grid}
}
\value{
  Returns a list with two vectors:
  \item{x}{Processed longitudes}
  \item{y}{Processed latitudes}
}
\note{
  These functions are based on the respective fortran code from the COSMO model (numerical weather prediction model of the DWD, Germany).
}
\references{missing reference to COSMO code}
\author{Jonas Bhend: jonas.bhend@env.ethz.ch}
\examples{
    # Location of c(0,0) in rotated coordinates
    rot2geo(-170, 32.5, 0, 0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
\keyword{ utilities }% __ONLY ONE__ keyword per line
