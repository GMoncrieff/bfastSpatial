#' @title Apply bfmastmonitor on a single pixel with a covariate
#' 
#' @description Apply bfastmonitor (bfm) on a single pixel of known cell index (by supplying a single numeric value for cell), xy coordinates (\code{cell=c(x, y)}), or interactively by clicking on a plot (by setting \code{interactive=TRUE}). Outputs a list with (1) an object of class 'bfastmonitor' and (2) the resulting cell number (useful for follow-up analysis).
#' 
#' @param x RasterBrick with raster time series data.
#' @param start Numeric. Vector of length = 2 representing the start of the monitoring period (in the format c(year, julian day))
#' @param monend Numeric. Optional: the end of the monitoring period (in the format c(year, julian day)), at which point the time series will be trimmed.
#' @param cell Numeric. Can be one of: (1) a numeric of length 1 indicating the raster cell to be observed; (2) a numeric of length 2 representing the (x,y) coordinate of the raster cell to be observed. Can also be omitted, in which case 'interactive' must be set to TRUE (see below)
#' @param f Numeric. Factor by which to rescale values before running \code{bfastmonitor}. Defaults to 1 (no rescaling)
#' @param min.thresh Numeric. Optional: A minimum threshold below which NA's are assigned to data points. NOTE: the threshold is applied \emph{before} rescaling the data by \code{f} (see above)
#' @param sensor Character. Optional: Limit analysis to data from one or more sensors. Can be one or more of \code{c("ETM+", "ETM+ SLC-on", "ETM+ SLC-off", "TM", "OLI")} according to the sensor information returned by \code{\link{getSceneinfo}}
#' @param interactive Logical. Select cell by clicking on an already plotted map? Defaults to \code{FALSE}. If \code{FALSE}, a value must be assigned to \code{cell} (see above).
#' @param plot Logical. Plot the result? Defaults to \code{FALSE}.
#' @param xreg See \code{\link{bfastmonitor}}
#' @param ... Arguments to be passed to \code{\link{bfastmonitor}}
#' 
#' @return A list with the following components: 1) $bfm - an object of class 'bfastmonitor' (see \code{\link{bfastmonitor}}) 2) $cell - the cell index (an integer of length 1). This can be used to run \code{bfmPixel} again on the same pixel (with different parameters) without having to click on a plot again to find the same pixel (in that case, be sure to set interactive=FALSE for subsequent trials!).
#' 
#' @details \code{bfmPixel_xreg} is theoretically designed to work on any generic raster time series, as long as a \code{dates} vector is provided. In the absence of a \code{dates} vector, \code{names(x)} should correspond exactly to respective Landsat scene ID's. In this case, \code{\link{getSceneinfo}} is used to extract a dates vector, and subset by sensor if desired.
#' 
#' @author Ben DeVries
#' 
#' @examples
#' 
#'  \dontrun{
#' # load in time series raster data
#' data(tura)
#' 
#' xreg <- cellStats(tura,stat= 'mean')
#' xreg<-as.numeric(xreg)
#' xreg <- zo::na.approx(xreg)
#' 
#' # run bfm on a pixel of known (x,y) location
#' bfm <- bfmPixel(tura, cell=c(820900, 831340), start=c(2005, 1),xreg=xreg)
#' plot(bfm$bfm)
#' }
#' 
#' @seealso \code{\link{bfastmonitor}}, \code{\link{bfmSpatial}}
#' 
#' @import raster
#' @import bfast
#' @export


bfmPixel_xreg <- function (x, dates=NULL, start, monend=NULL, xreg, cell=NULL, f=1, min.thresh=NULL, sensor=NULL, interactive=FALSE, plot=FALSE, ...) 
{
    
    
    # select cell from the input raster brick x in 1 of 3 ways:
    # 1) interactively (by clicking on an already plotted map)
    # 2) by supplying the cell index as an integer of length=1
    # 3) by supplying a vector of length=2 representing the (x,y) coordinates
    if (interactive) { # condition 1:
        cell <- as.data.frame(click(x, n=1, id=TRUE, cell=TRUE, show=FALSE))$cell
    } else { # conditions 2 and 3:
        cell <- ifelse(length(cell)==2, cellFromXY(x, t(as.matrix(cell))), cell)
    }
    
    # extract pixel time series
    pixelts <- as.vector(x[cell])
    
    
    # optional: apply a threshold (if supplied)
    if (!is.null(min.thresh))
        pixelts[pixelts <= min.thresh] <- NA
    
    # optional: rescale values
    if(f != 1)
        pixelts <- pixelts * f
    
    ##
    pixelts <- cbind(pixelts,xreg)
    
    # convert to a bfast ts object
    pixelts <- bfastts(pixelts, dates, type=c("irregular"))
    
    # optional: trim ts if monend is supplied
    if(!is.null(monend))
        pixelts <- window(pixelts, end=monend)
    
    # run bfm on the pixel time series
    bfm <- bfastmonitor(data=pixelts, start=start, response ~ trend + xreg, ...)
    
    # plot if plot=TRUE
    if(plot)
        plot(bfm)
    
    # return a list with (1) a bfm object, and (2) the cell number (for follow-up)
    return(list(bfm=bfm, cell=cell))
}