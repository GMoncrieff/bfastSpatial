% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/saveExtent.R
\name{saveExtent}
\alias{saveExtent}
\title{Draw an extent on a raster plot and print command to console}
\usage{
saveExtent(..., digits = 0)
}
\arguments{
\item{...}{Arguments to be passed to \code{\link{drawExtent}}}

\item{digits}{integer indicating the number of decimal places (passed to \code{\link{round}})}
}
\value{
an object of class \code{extent} with the side effect of having the \code{extent} command with resulting parameters printed to the console for reproducability.
}
\description{
Identical to \code{\link{drawExtent}} except that the resulting extent is printed to the console with a formatted \code{extent} command which can be copied and pasted to a script for reproducability.
}
\examples{
\dontrun{
data(tura)
plot(tura, 6)
e <- saveExtent() # draw extent on the plot, and associated extent() command will be printed to console, while extent is saved to workspace
}

}
\seealso{
\code{\link{drawExtent}}
}
\author{
Ben DeVries
}
