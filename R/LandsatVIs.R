# Loic Dutrieux
# April 2014

# Each function needs to return a vector of indices (band numbers) and a function to be passed to overlay


# NDVI --------------------------------------------------------------------

.ndvi <- function() {
    ind <- c('band3','band4')
    fun <- function(x1, x2) {
        ndvi <- 10000 * (x2 - x1)/(x2 + x1)
        return(ndvi)
    }
    return(list(ind=ind,
                fun=fun))
}


# EVI ---------------------------------------------------------------------

.evi <- function() {
    ind <- c('band1','band3','band4')
    fun <- function(x1, x3, x4){ 
        evi <- 10000 * 2.5 * (x4/10000 - x3/10000)/(x4/10000 + 6 * x3/10000 - 7.5 * x1/10000 + 1)
        return(evi)
    }
    return(list(ind=ind,
                fun=fun))
    
}

# SAVI --------------------------------------------------------------------
# Soil-adjusted Vegetation Index:
# SAVI = (1 + L)(B4 - B3)/(B4 + B3 + L)
# where L ~ [0,1] depending on climate, and is often assumed to be 0.5

.savi <- function(L=0.5) {
    ind <- c('band3','band4')
    fun <- function(x1, x2){ 
        savi <- 10000 * (1 + L) * (x2/10000 - x1/10000) / (x2/10000 + x1/10000 + L)
        return(savi)
    }
    return(list(ind=ind,
                fun=fun))
}

# NBR ---------------------------------------------------------------------
# Normalized Burn Ratio:
# NBR = (B4 - B7) / (B4 + B7)

.nbr <- function() {
    ind <- c('band4','band7')
    fun <- function(x1, x2) {
        ndvi <- 10000 * (x1 - x2)/(x1 + x2)
        return(ndvi)
    }
    return(list(ind=ind,
                fun=fun))
}

# Tasseled Cap Components ------------------------------------------------- 

.tcbright <- function(sensor) { 
    ind <- c('band1','band2','band3','band4','band5','band7')  
    # make compatible with getSceneinfo() output 
    if(sensor %in% c("ETM+", "ETM+ SLC-on", "ETM+ SLC-off")) 
        sensor <- 7 
    if(sensor == "TM") 
        sensor <- 5 
    
    if(sensor == 5) { 
        tc_coef <- c(0.2043, 0.4158, 0.5524, 0.5741, 0.3124, 0.2303) 
    } else if (sensor == 7) { 
        tc_coef <- c(0.3561, 0.3972, 0.3904, 0.6966, 0.2286, 0.1596)
    } 
    
    fun <- function(x1, x2, x3, x4, x5, x7) { 
        tcbright <- sum(c(x1, x2, x3, x4, x5, x7) * tc_coef) 
    } 
    
    return(list(ind=ind, 
                fun=fun)) 
} 


.tcgreen <- function(sensor) { 
    ind <- c('band1','band2','band3','band4','band5','band7') 
    # make compatible with getSceneinfo() output 
    if(sensor %in% c("ETM+", "ETM+ SLC-on", "ETM+ SLC-off")) 
        sensor <- 7 
    if(sensor == "TM") 
        sensor <- 5 
    
    if(sensor == 5) { 
        tc_coef <- c(-0.1603, -0.2819, -0.4934,  0.7940, 0.0002, -0.1446) 
    } else if (sensor == 7) { 
        tc_coef <- c(-0.3344, -0.3544, -0.4556,  0.6966, -0.0242,-0.2630)
    } 
    
    fun <- function(x1, x2, x3, x4, x5, x7) { 
        tcbright <- sum(c(x1, x2, x3, x4, x5, x7) * tc_coef) 
    } 
    
    return(list(ind=ind, 
                fun=fun)) 
} 


.tcwet <- function(sensor) { 
    ind <- c('band1','band2','band3','band4','band5','band7') 
    # make compatible with getSceneinfo() output 
    if(sensor %in% c("ETM+", "ETM+ SLC-on", "ETM+ SLC-off")) 
        sensor <- 7 
    if(sensor == "TM") 
        sensor <- 5 
    
    if(sensor == 5) { 
        tc_coef <- c(0.0315,  0.2021,  0.3102,  0.1594, 0.6806, -0.6109) 
    } else if (sensor == 7) { 
        tc_coef <- c(0.2626, 0.2141, 0.0926, 0.0656, -0.7629, -0.5388) 
    } 
    
    fun <- function(x1, x2, x3, x4, x5, x7) { 
        tcbright <- sum(c(x1, x2, x3, x4, x5, x7) * tc_coef) 
    } 
    
    return(list(ind=ind, 
                fun=fun)) 
} 