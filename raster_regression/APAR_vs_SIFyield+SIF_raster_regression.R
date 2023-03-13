library(terra)
library(lm.beta)

### Multiple regression

y_file     <- "G:/ChloFluo/input/APARchl/1deg/apar.2019.8-day.1deg.nc"
x1_file    <- "G:/ChloFluo/input/yield/1deg/yield.2019.8-day.1deg.nc"
x2_file    <- "G:/ChloFluo/input/SIF/1deg/SIFqc.8day.1deg.CF80.2019.nc"
y_name     <- "aparchl"
x1_name    <- "sif_yield"
x2_name    <- "sif743_qc"
out_dir    <- "G:/ChloFluo/comps/aparchl_vs_yield_sif/raster_regressions"
out_name   <- "APARchl_vs_SIFyield+SIF.v01.1deg.CF80.2019.clipfill"
f_name     <- NA # Filter by value. Example, error, std, or n. If none use NA.
f_thresh   <- 30  # Values => will be kept


rastlm <- function(x) {
  full <- length(x)
  third <- full / 3
  
  if (all(is.na(x[1:third])) || all(is.na(x[(third + 1):(third * 2)])) || all(is.na(x[(third * 2 + 1):full]))) { 
    
    return(c(NA,NA,NA,NA,NA,NA,NA))
    
  } else { 
    reg     <- lm(x[1:third] ~ x[(third + 1):(third * 2)] + x[(third * 2 + 1):full])
    s       <- summary(reg)
    if (is.na(s$coefficients[2]) || is.na(s$coefficients[3])) {
      
      print("cor not good")
      return(c(NA,NA,NA,NA,NA,NA,NA))
      
    } else {
      reg_standard <- lm.beta(reg)
      standard_s   <- summary(reg_standard)
      
      pval_x1          <- s$coefficients[2,4]
      pval_x2          <- s$coefficients[3,4]
      beta_x1          <- s$coefficients[2]
      beta_x2          <- s$coefficients[3]
      beta_standard_x1 <- standard_s$coefficients[2,2]
      beta_standard_x2 <- standard_s$coefficients[3,2]
      n                <- nobs(reg)
      
      return(c(pval_x1, pval_x2, beta_x1, beta_x2, beta_standard_x1, beta_standard_x2, n)) 
    }
  }
}

rast_reg <- function(y_file, x1_file, x2_file, y_name, x1_name, x2_name, out_dir, out_name) {
  
  y  <- rast(y_file, subds = y_name)
  x1 <- rast(x1_file, subds = x1_name)
  x2 <- rast(x2_file, subds = x2_name)
  
  y  <- mask(y, x1)
  y  <- mask(y, x2)
  x1 <- mask(x1, y)
  x1 <- mask(x1, x2)
  x2 <- mask(x2, y)
  x2 <- mask(x2, x1)
  
  # Filter as needed
  if (!is.na(f_name)) {
    print(paste0("Filtering data using: ", f_name))
    print(paste0("Filter threshold is: ", f_thresh))
    
    f <- brick(in_file, varname = f_name)
    f[f < f_thresh] <- NA
    
    y  <- mask(y, x1)
    y  <- mask(y, x2)
    x1 <- mask(x1, y)
    x1 <- mask(x1, x2)
    x2 <- mask(x2, y)
    x2 <- mask(x2, x1)
    
  }
  
  # Combine bricks into single brick. lm convention is y~x
  yx <- c(y, x1, x2)
  
  lm.result <- app(yx, rastlm)
  
  writeRaster(lm.result[[1]], paste0(out_dir, "/", out_name, "_Pval_", x1_name, ".tif"), overwrite = TRUE)
  writeRaster(lm.result[[2]], paste0(out_dir, "/", out_name, "_Pval_", x2_name, ".tif"), overwrite = TRUE)
  writeRaster(lm.result[[3]], paste0(out_dir, "/", out_name, "_Beta_", x1_name, ".tif"), overwrite = TRUE)
  writeRaster(lm.result[[4]], paste0(out_dir, "/", out_name, "_Beta_", x2_name, ".tif"), overwrite = TRUE)
  writeRaster(lm.result[[5]], paste0(out_dir, "/", out_name, "_Beta_Standard_", x1_name, ".tif"), overwrite = TRUE)
  writeRaster(lm.result[[6]], paste0(out_dir, "/", out_name, "_Beta_Standard_", x2_name, ".tif"), overwrite = TRUE)
  writeRaster(lm.result[[7]], paste0(out_dir, "/", out_name, "_Nobs.tif"), overwrite = TRUE)
  
  remove(x1, x2, y, yx) # get it out of memory
  
}

rast_reg(y_file, x1_file, x2_file, y_name, x1_name, x2_name, out_dir, out_name)
