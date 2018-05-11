\name{ncdf_times}
\alias{ncdf_times}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to extract dates and times from NetCDF file opened with open.ncdf().}
\description{
  This function extracts the date and time information from NetCDF files and converts it either to the internal Date format or to the udunits format if udunits is present.
}
\usage{
  ncdf_times(nc, as.Rdate = TRUE, force = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nc}{The connection to the netcdf file as with open.ncdf()}
  \item{as.Rdate}{switch to convert output to an object of class "Date"}
  \item{force}{switch to force conversion to an object of class "Date" even if non-standard calendars are used}
}
\details{
  ...  
}
\value{
  A vector of Dates or Datetimes (e.g. \code{is.Date()} or \code{is.POSIXct()}) depending on the units in the NetCDF file (if units are "days/months/years since ...", then the output is of class "Date", otherwise, it is of class "POSIXct") or a list of udunit dates.
}
\author{Jonas Bhend: jonas.bhend@env.ethz.ch}
\section{Warning}{
  The forced conversion to output of class "Date" will lead to wrong dates if non-standard calendars are used!

  At present, the conversion of times with a reference to time zones other than GMT could be flawed (R objects of class date are not aware of time zones). Please carefully check the times when working with data from different time zones.
}
\seealso{\code{\link[ncdf]{ncdf.open}},\code{\link{as.Date}}, \code{\link{as.POSIXct}}}
\examples{
  # Define the time
  dimtime   <- dim.def.ncdf( "time", "days since 1990-01-01", 0:364 )

  # define a random variable
  vart2m    <- var.def.ncdf("T_2M", "K", dimtime, -1e20, 
  longname="Near-surface air temperature")

  # Create a netCDF file with this variable
  ncnew     <- create.ncdf( "ncdf_times_example.nc", vart2m )

  # Write some values to this variable on disk.
  put.var.ncdf( ncnew, vart2m, rnorm(365, mean=10, sd=2) + 
    12*cos(seq(-pi,pi,length=365)), start=1, count=365 )

  close.ncdf(ncnew)

  # Now open the file and read its data
  nc        <- open.ncdf("ncdf_times_example.nc")
  times     <- ncdf_times(nc)
  print("here is the data in the file:")
  head(times) 
  tail(times)
  t2m       <- get.var.ncdf(nc, "T_2M")
  close.ncdf( nc )
  file.remove("ncdf_times_example.nc")
    
  plot(times, t2m, type="l")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ attribute }
\keyword{ chron }% __ONLY ONE__ keyword per line
