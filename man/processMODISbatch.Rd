% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processMODISbatch.R
\name{processMODISbatch}
\alias{processMODISbatch}
\title{Wrapper/batcher to pre-process MODIS data in batch mode}
\usage{
processMODISbatch(x, pattern = NULL, data_SDS, QC_SDS, bit = FALSE,
  QC_val, fill = NULL, outdir, mosaic = TRUE, mc.cores = 1)
}
\arguments{
\item{x}{Character. Directory where the data is located. Or list of file names.}

\item{pattern}{Only useful if x is of length 1. See \link{list.files} for more details}

\item{data_SDS}{See \link{cleanMODIS}.}

\item{QC_SDS}{See \link{cleanMODIS}.}

\item{bit}{See \link{cleanMODIS}.}

\item{QC_val}{\link{cleanMODIS}.}

\item{fill}{\link{cleanMODIS}.}

\item{outdir}{Character. Directory where the output should be written.}

\item{mosaic}{Logical. When working with several tiles, should these be mosaicked or kept as separate output files. Default is \code{TRUE}}

\item{mc.cores}{Numeric. Number of cores to use for processing. See \link{mclapply}}
}
\value{
This function is used for its side effect of producing cleaned modis data, hence written to disk at \code{outdir} location
}
\description{
Processes MODIS data for subsequentent use in time-series analysis. Performs MODIS cleaning and mosaicking operations with parallel support. This function requires rgdal to be configured with HDF4 driver.
}
\seealso{
\code{\link{cleanMODIS}}, and \code{\link{timeStackMODIS}} for stacking.
}
