\name{rbfun}
\alias{rbfun}
\alias{gbfun}
\alias{colfun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to provide colour sequences }
\description{
 The functions provide a red to blue or green to brown colour sequence.
}
\usage{
rbfun(x, start = 0.4, log = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Number of colours}
  \item{start}{starting point (0 is black, 1 is white)}
  \item{log}{logarithmic scaling (more difference in light colours)}
}
\value{
a vector of rgb colours in R colour format
}
\author{Jonas Bhend: jonas.bhend@env.ethz.ch}
\examples{

image(array(1:9, c(1,9)), breaks=0.5 + 0:9, col=rbfun(9), axes=F, xlab="", ylab="")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{color}
