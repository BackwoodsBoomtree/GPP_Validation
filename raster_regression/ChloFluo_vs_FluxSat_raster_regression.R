library(raster)

cf_file  <- "G:/ChloFluo/product/v02/clipfill/ChloFluo.GPP.v02.1deg.CF80.2019.clipfill.nc"
gpp_file <- "G:/FluxSat/monthly/1deg/GPP_FluxSat_8day_1deg_v2_2019.nc"
y_name   <- "gpp"
x_name   <- "GPP"
out_dir  <- "G:/ChloFluo/comps/fluxsat/raster_regressions"
out_name <- "ChloFluo_vs_FluxSat.v02.1deg.CF80.2019.clipfill"
f_name   <- NA # Filter by value. Example, error, std, or n. If none use NA.
f_thresh <- 30  # Values => will be kept


y <- brick(cf_file, varname = y_name)
x <- brick(gpp_file, varname = x_name)

# # Create single brick from input files
# for (i in 1:(length(gpp_file))) {
#   
#   x_tmp <- brick(gpp_file[i], varname = x_name)
#   
#   if (i == 1) {
#     x <- x_tmp
#   } else {
#     x     <- addLayer(x, x_tmp)
#   }
# }
# 
# # Do temporal aggregation (non-leap year)
# for (i in seq(1, 365, by = 8)) {
#   
#   if (i != 361) {
#     x_mean_tmp <- mean(x[[i : (i + 7)]], na.rm = TRUE)
#   } else {
#     x_mean_tmp <- mean(x[[i : (i +  4)]], na.rm = TRUE)
#   }
#   
#   if (i == 1) {
#     x_8day <- x_mean_tmp
#   } else {
#     x_8day <- addLayer(x_8day, x_mean_tmp)
#   }
# }

# x <- aggregate(x_8day, 20, fun = mean, na.rm = TRUE)

# writeRaster(x, "G:/FluxSat/1deg/GPP_FluxSat_8day_1deg_v2_2019.nc4", overwrite = TRUE, format="CDF", varname="GPP", varunit = "g C m-2 d-1", 
#             longname="Gross Primary Production", xname = "Longitude", yname = "Latitude", zname = "Time")


rastlm <- function(x) {
  full <- length(x)
  half <- full / 2
  
  if (all(is.na(x[1:half])) || all(is.na(x[(half + 1):full]))){ 
    
    return(c(NA,NA,NA,NA,NA,NA))
    
  } else { 
    reg       <- lm(x[1:half] ~ x[(half +1):full])
    s         <- summary(reg)
    r2        <- s$r.squared
    pval      <- s$coefficients[8]
    slope     <- s$coefficients[2]
    intercept <- s$coefficients[1]
    rmse      <- sqrt(mean(s$residuals^2))
    n         <- nobs(reg)
    
    return(c(r2, pval, slope, intercept, rmse, n)) 
  }
}


rast_reg <- function(y, x, out_dir, out_name) {
  
  y <- mask(y, x)
  x <- mask(x, y)
  
  # Filter as needed
  if (!is.na(f_name)) {
    print(paste0("Filtering data using: ", f_name))
    print(paste0("Filter threshold is: ", f_thresh))
    
    f <- brick(in_file, varname = f_name)
    f[f < f_thresh] <- NA
    
    y <- mask(y, f)
    x <- mask(x, f)
    
  }
  
  # Combine bricks into single brick. lm convention is y~x
  yx <- stack(y, x)
  
  beginCluster(12)
  lm.result <- clusterR(yx, calc, args = list(fun = rastlm))
  endCluster()
  
  # Create dir
  if (!dir.exists(out_dir)){
    dir.create(out_dir, recursive = TRUE)
    print("Dir created!")
  } else {
    print("Dir already exists!")
  }
  
  writeRaster(lm.result[[1]], paste0(out_dir, "/", out_name, "_Rsquare.tif"), overwrite = TRUE)
  writeRaster(lm.result[[2]], paste0(out_dir, "/", out_name, "_Pval.tif"), overwrite = TRUE)
  writeRaster(lm.result[[3]], paste0(out_dir, "/", out_name, "_Slope.tif"), overwrite = TRUE)
  writeRaster(lm.result[[4]], paste0(out_dir, "/", out_name, "_Intercept.tif"), overwrite = TRUE)
  writeRaster(lm.result[[5]], paste0(out_dir, "/", out_name, "_RMSE.tif"), overwrite = TRUE)
  writeRaster(lm.result[[6]], paste0(out_dir, "/", out_name, "_Nobs.tif"), overwrite = TRUE)
  
  remove(x, y, yx) # get it out of memory
  
}

rast_reg(y, x, out_dir, out_name)
