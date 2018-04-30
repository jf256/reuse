\name{plot_colourbar}
\alias{plot_colourbar}
\alias{plot_colourbar.default}
\alias{plot_colourbar.plotmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to plot a colourbar}
\description{
  Generic function to plot a horizontal or vertical colourbar with axis on the specified side.
}
\usage{
  plot_colourbar(levs, cols, side=1, 
    ylab="",xlab="", nn=1, center=F, cex.axis=1, 
    sea.col=NULL, ...)
    
  plot_colourbar(x.list, incl.units=T,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{levs}{Levels (boundaries)}
  \item{cols}{colours to use}
  \item{side}{side at which the labelling takes place}
  \item{nn}{how many labels are set (every nn'th label is set)}
  \item{center}{switch for labelling in the middle of the colour box)}
  \item{sea.col}{colour for the lowest box, gets a different labelling)}
  \item{x.list}{A list of class "plotmap" containing the levels and colours and additional information as returned by \code{plotmap} or \code{plotmap_rot}}
  \item{incl.units}{allows for the long-name and units attribute of the plotted variable to be added to the axis annotation)}
}
\details{
  There must be one more levels than colours including the colour for water (\code{sea.col}). 
  
  The enumeration of side is consistent with \code{axis()}, e.g. from 1 to 4 clockwise from bottom.
  
  Only every nn'th label is drawn. 
  
  Ticks are omitted if centered labels are drawn.
}
\author{Jonas Bhend: jonas.bhend@env.ethz.ch}

\seealso{\code{\link{plotmap}}, \code{\link{plotmap_rot}}}
\examples{
  plot_colourbar(levs=0:11, col=rbfun(10), 
      xlab="imaginary units", sea.col="lightblue")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
