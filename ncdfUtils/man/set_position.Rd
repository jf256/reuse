\name{set_position}
\alias{set_position}
\alias{reset_position}
\title{Set the position of a mini-plot within the active plot region and reset to mother plot.}
\description{
  Sets the position of a mini plot within the active plot region at user coordinates x and y.
}
\usage{
set_position(x, y, size = 0.1, bg = NULL)
}
\usage{
reset_position()
}
\arguments{
  \item{x}{x coordinate value of center of mini plot}
  \item{y}{y coordinate value of center of mini plot }
  \item{size}{size relative to the existing plot}
  \item{bg}{Background colour for the new plot}
}
\details{
  The size can be given either as a single value 
  or a pair for centered mini plots, or as a vector
  of length 4 for arbitrary positioned plots. The 
  sequence of the arguments is c(xleft, xright,
  ybottom, ytop) (e.g. size = c(0,0.2,0,0.2) positions 
  the plot to the upper right of the indicated position).
  
  \code{set_position} writes a hidden variable named \code{.par_old}
  to the parent environment, which is subsequently used by
  \code{reset_position}.
}
\author{Jonas Bhend: <email: jonas.bhend@env.ethz.ch>}
\seealso{\code{\link{par}}, \code{\link{map}}}
\examples{

    # add random numbers in Northern Germany
    map(xlim=c(-15,35), ylim=c(35,75))
    set_position(10,55,0.2,bg="white")
    plot(rnorm(50), pch=".", xaxt="n", yaxt="n", xlab="", ylab="")
    reset_position()

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aplot }
\keyword{ hplot }
