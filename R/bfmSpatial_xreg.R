
#' @title Function to run bfastmonitor on any kind of raster brick using a covariate, with parallel support
#' 
#' @description Implements bfastmonitor function, from the bfast package on any kind of rasterBrick object. Time information is provided as an extra object and the time series can be regular or irregular.
#' 
#' @param x rasterBrick or rasterStack object, or file name to a multilayer raster object stored on disk.
#' @param dates A date vector. The number of dates must match the number of layers of x.
#' @param pptype Character. Type of preprocessing to be applied to individual time series vectors. The two options are 'irregular' and '16-days'. See \code{\link{bfastts}} for more details.
#' @param start See \code{\link{bfastmonitor}}
#' @param monend Numeric. Optional: end of the monitoring period in the format c(year, julian day). All raster data after this time will be removed before running \code{bfastmonitor}
#' @param formula See \code{\link{bfastmonitor}}
#' @param order See \code{\link{bfastmonitor}}
#' @param lag See \code{\link{bfastmonitor}}
#' @param slag See \code{\link{bfastmonitor}}
#' @param xreg See \code{\link{bfastmonitor}}
#' @param history See \code{\link{bfastmonitor}}
#' @param type See \code{\link{bfastmonitor}}
#' @param h See \code{\link{bfastmonitor}}
#' @param level See \code{\link{bfastmonitor}}
#' @param mc.cores Numeric. Number of cores to be used for the job.
#' @param returnLayers Character. Result layers to be returned. Can be any combination of \code{c("breakpoint", "magnitude", "error", "history", "r.squared", "adj.r.squared", "coefficients")}. By default, \code{breakpoint}, \code{magnitude} and \code{error} are returned by the function. See \code{details} for more information.
#' @param ... Arguments to be passed to \code{\link{mc.calc}}
#' 
#' @return A rasterBrick with layers depending on what has been supplied to \code{returnLayers}. By default, 3 layers are returned: (1) breakpoint: timing of breakpoints detected for each pixel; (2) magnitude: the median of the residuals within the monitoring period; (3) error: a value of 1 for pixels where an error was encountered by the algorithm (see \code{\link{try}}), and NA where the method was successfully run. See \code{\link{bfastmonitor}} for more information on the other possible layers.
#' 
#' @details
#' \code{bfmSpatial_xreg} applies \code{\link{bfastmonitor}} over a raster time series. For large raster datasets, processing times can be long. Given the number of parameters that can be set, it is recommended to first run \code{\link{bfmPixel}} over some test pixels or \code{bfmSpatial_xreg} over a small test area to gain familiarity with the time series being analyzed and to test several parameters.
#' 
#' Note that there is a difference between the \code{monend} argument included here and the \code{end} argument passed to \code{\link{bfastmonitor}}. Supplying a date in the format \code{c(year, Julian day)} to \code{monend} will result in the time series being trimmed \emph{before} running \code{\link{bfastmonitor}}. While this may seem identical to trimming the resulting \code{bfastmonitor} object per pixel, trimming the time series before running \code{bfastmonitor} will have an impact on the change magnitude layer, which is calculated as the median residual withint the entire monitoring period, whether or not a breakpoint is detected.
#' 
#' \code{returnLayers} can be used to specify which \code{bfasmonitor} results to return. Regardless of which parameters are assigned, the output layers will always follow the order: \code{c("breakpoint", "magnitude", "error", "history", "r.squared", "adj.r.squared", "coefficients")}. This is important if \code{mc.cores} is set to be greater than 1, since this causes the layer names in the output brick to be lost, so it is important to know which layers have been requested and in which order they will be exported. Note that if "coefficients" is included, the output will include the following: "(Intercept)" and any trend and/or harmonic coefficients depending on the values of \code{formula} and \code{order}.
#' 
#' @author Loic Dutrieux and Ben DeVries
#' @import bfast
#' @import parallel
#' @import raster
#' 
#' @seealso \code{\link{bfastmonitor}}, \code{\link{bfmPixel}}
#' 
#' @examples
#' \dontrun{
#' # load tura dataset
#' data(tura)
#' 
#' xreg <- cellStats(tura,stat= 'mean')
#' xreg<-as.numeric(xreg)
#' xreg <- zo::na.approx(xreg)
#' # run BFM over entire time series with a monitoring period starting at the beginning of 2009
#' t1 <- system.time(bfm <- bfmSpatial_xreg(tura, start=c(2009, 1),xreg=xreg))
#' 
#' }
#' 
#' \dontrun{
#' # with multi-core support (see ?mc.calc)
#' t2 <- system.time(bfm <- bfmSpatial_xreg(tura, start=c(2009, 1),xreg=xreg, mc.cores=4))
#' # difference processing time
#' t1 - t2
#' }
#' 
#' # plot the result
#' plot(bfm)
#' 
#' @export
#' 


# Author: Loic Dutrieux
# January 2014

bfmSpatial_xreg <- function(x, dates=NULL, xreg, pptype='irregular', start, monend=NULL,
                            formula = response ~ trend + xreg, order = 3, lag = NULL, slag = NULL,
                            history = c("ROC", "BP", "all"),
                            type = "OLS-MOSUM", h = 0.25, end = 10, level = 0.01, mc.cores=1, returnLayers = c("breakpoint", "magnitude", "error"), ...) {
    
    if(is.character(x)) {
        x <- brick(x)
    }
    
    
    # determine length of coefficient vector
    # = intercept [+ trend] [+ harmoncos*order] [+ harmonsin*order]
    coef_len <- 1 # intercept
    modterms <- attr(terms(formula), "term.labels")
    if("trend" %in% modterms)
        coef_len <- coef_len + 1
    if("harmon" %in% modterms)
        coef_len <- coef_len + (order * 2) # sin and cos terms
    
    fun <- function(x) {

        x <- cbind(x,xreg)
        
        # convert to bfast ts
        ts <- bfastts(x, dates=dates, type=pptype)
        
        #optional: apply window() if monend is supplied
        if(!is.null(monend))
            ts <- window(ts, end=monend)
        
        # run bfastmonitor(), or assign NA if only NA's (ie. if a mask has been applied)
        if(!all(is.na(ts))){
            bfm <- try(bfastmonitor(data=ts, start=start,
                                    formula=formula,
                                    order=order, lag=lag, slag=slag,
                                    history=history,
                                    type=type, h=h,
                                    end=end, level=level), silent=TRUE)
            
            # assign 1 to error and NA to all other fields if an error is encountered
            if(class(bfm) == 'try-error') {
                bkpt <- NA
                magn <- NA
                err <- 1
                history <- NA
                rsq <- NA
                adj_rsq <- NA
                coefficients <- rep(NA, coef_len)
            } else {
                bkpt <- bfm$breakpoint
                magn <- bfm$magnitude
                err <- NA
                history <- bfm$history[2] - bfm$history[1]
                rsq <- summary(bfm$model)$r.squared
                adj_rsq <- summary(bfm$model)$adj.r.squared
                coefficients <- coef(bfm$model)
            }
        } else {
            bkpt <- NA
            magn <- NA
            err <- NA
            history <- NA
            rsq <- NA
            adj_rsq <- NA
            coefficients <- rep(NA, coef_len)
        }
        res <- c(bkpt, magn, err, history, rsq, adj_rsq)
        names(res) <- c("breakpoint", "magnitude", "error", "history", "r.squared", "adj.r.squared")
        res <- res[which(names(res) %in% returnLayers)]
        if("coefficients" %in% returnLayers)
            res <- c(res, coefficients)
        return(res)
    }
    
    out <- mc.calc(x=x, fun=fun, mc.cores=mc.cores, ...)
    
    return(out)
}