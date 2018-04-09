\name{shaded_contour}
\alias{shaded_contour}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to contour-plot matrizes with shading.
}
\description{
Unlike image(), shaded_contour() interpolates the matrix before plotting. The 
}
\usage{
shaded_contour(x = seq(0, 1, length.out = nrow(z)), 
    y = seq(0, 1, length.out = ncol(z)), 
    z, xlim = x[c(1, length(x))], ylim = y[c(1, length(y))], 
    zlim = range(z, finite = TRUE), levels = pretty(zlim, nlevels),
     nlevels = 20, col = rbfun(length(levels) - 1), 
    xaxs = "i", yaxs = "i", add = FALSE, xlab, ylab, las = 1, 
    borders = TRUE, l.col = grey(0.7), l.lwd = 1, block.na=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,y}{
locations of grid lines at which the values in z are measured
}
  \item{z}{
a matrix containing the values to be plotted ('NA's are allowed). Note that 'x' as a list with elements x, y, and z can be used instead of 'z' for convenience.
}
  \item{xlim, ylim}{
ranges for the plotted 'x' and 'y' values, defaulting to the ranges of 'x' and 'y'.
}
  \item{zlim}{
the minimum and maximum 'z' values for which colors should be
plotted, defaulting to the range of the finite values of 'z'.
Each of the given colors will be used to color an equispaced
interval of this range. The _midpoints_ of the intervals
cover the range, so that values just outside the range will
be plotted.
}
  \item{levels}{
a set of levels which partition the values of 'z'. Must give one more level than colors.
}
  \item{nlevels}{
approximate number of levels to be used if levels is not set.
}
  \item{col}{
an explicit set of colors to be used in the plot.
}
  \item{las}{
the style of axis labels. See help(par) for further info.
}
  \item{borders}{
logical, whether or not to draw boundary lines in between the shaded areas.
}
  \item{l.col}{
color used for the boundaries if borders is set to TRUE.
}
  \item{l.lwd}{
line width of the boundaries if borders is set to TRUE.
}
  \item{block.na}{
switch whether block style missing values or default (partly interpolated) 
missing values should be plotted.
}
  \item{\dots}{
further arguments to be passed to plot().
}
}
\value{
shaded_contour() returns a list of class plotmap containing the levels and colors used in the plot. This list can directly be passed to plot_colourbar() for the appropriate colorbar.
}
\author{
Jonas Bhend: jonas.bhend@env.ethz.ch
}

\seealso{\code{\link{plot_colourbar}},\code{\link{image}}, \code{\link{filled.contour}}}
\examples{

shaded_contour(volcano)

layout(matrix(1:2, 2, 1), height=c(5,2))
par(mar=c(3,3,1,1))
tmp <- shaded_contour(volcano, col=terrain.colors(21), borders=FALSE)
tmp$units    <- "m"
plot_colourbar(tmp)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
