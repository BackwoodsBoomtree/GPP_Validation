library(terra)
library(viridis)
library(rgdal)
library(RColorBrewer)

out_name <- "G:/ChloFluo/comps/tower_comparisons_all.pdf"

cf_file     <- "G:/ChloFluo/product/v02/clipfill/ChloFluo.GPP.v02.1deg.CF80.2019.clipfill.nc"
fcom_file   <- "G:/FluxCom/RS/GPP.RS_V006.FP-ALL.MLM-ALL.METEO-NONE.720_360.8daily.2019.nc"
fsat_files  <- "G:/FluxSat/monthly/1deg/GPP_FluxSat_8day_1deg_v2_2019.nc"
sif_file    <- "G:/ChloFluo/input/SIF/1deg/SIFqc.8day.1deg.veg_only.CF80.2019.nc"
tower_list  <- "G:/ChloFluo/comps/tower-data/Joiner_2020_Sites.csv"
tower_dir   <- "G:/ChloFluo/comps/tower-data/unzipped"

### Can be used for pvalues
round2 = function(x, n, p) {
  posneg = sign(x)
  z <- abs(x)*10^n
  z <- z + 0.5
  z <- trunc(z)
  z <- z/10^n
  z <- z*posneg
  
  if (p == TRUE) {
    if (z < 0.05 && z >= 0.01) {
      z <- "p < 0.05"
    } else if (z < 0.01) {
      z <- "p < 0.01"
    } else {
      z <- paste0("p = ", z)
    }
  }
  return(z)
}

# Monthly aggregation of model and sif data
to_month <- function(vals, year) {
  
  # Number of days for the year
  ndays <- ifelse(((year %% 100 != 0) & (year %%4 ==0)) | (year %% 400==0), 366 , 365)
  
  # DOY each value
  nn <- rep(1:length(vals), each = 8)[1:ndays]
  
  # day of year for each month
  m <- as.integer(format(as.Date(1:ndays, origin=paste0(year-1, "-12-31")), "%m"))
  
  # x describes for each day of the year, which layer to use, and which month it is.
  x <- cbind(layer=nn, month=m)
  
  # Now we can, for each month, determine how much of each layer is in that month
  weights <- table(x[,1], x[,2])
  
  s <- c()
  for (i in 1:12) {
    w <- weights[,i]
    x <- vals[which(w > 0)]
    ww <- w[w > 0] / 8
    s[i] <- weighted.mean(x, ww, na.rm = TRUE)
  }
  
  return(s)
  
}

# Grab tower data and compute annual mean
get_tower_gpp <- function(f) {
  df       <- read.csv(f, header = TRUE)
  
  # Take mean of DT and NT
  gpp_mean <- (df$GPP_NT_VUT_MEAN + df$GPP_DT_VUT_MEAN) / 2
  gpp_mean[gpp_mean == -9999] <- NA
  
  # Annual means
  n_years <- length(gpp_mean) / 12
  
  if (n_years > 1) {
    df_years <- data.frame(matrix(ncol = n_years, nrow = 12))
    
    for (i in 1:n_years) {
      df_years[,i] <- gpp_mean[(i * 12 - 11) : (i * 12)]
    }
    annual_mean <- rowMeans(df_years, na.rm = TRUE)
  } else {
    annual_mean <- gpp_mean
    df_years    <- gpp_mean
  }
  return(list(annual_mean, gpp_mean, df_years))
}

cf     <- rast(cf_file, subds = "gpp")
fcom   <- rast(fcom_file, subds = "GPP")
fsat   <- rast(fsat_files, subds = "GPP")
sif    <- rast(sif_file, subds = "sif_qc")
fcom   <- aggregate(fcom, 2, fun = mean, na.rm = TRUE)
towers <- read.csv(tower_list, header = TRUE)

### Get csv files for each tower
tower_dirs <- list.dirs(tower_dir, full.names = TRUE, recursive = FALSE)
for (i in 1:length(towers$site)) {
  site       <- towers$site[i] # from tower df
  dir_loc    <- grep(site, tower_dirs)
  
  if (length(dir_loc) != 0) {
    site_files    <- list.files(tower_dirs[dir_loc], full.names = TRUE, recursive = FALSE)
    file_loc      <- grep("FULLSET_MM_", site_files)
    file          <- site_files[file_loc]
    towers$data[i] <- file
  } else {
    towers$data[i] <- NA
  }
}
towers <- na.omit(towers)


## Plot variables
x            <- 1:12
x_lab        <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
y_lab_gpp    <- bquote("GPP (g C m"^"-2"*" d"^"-1"*")")
y_lab_sif    <- bquote("SIF"[Daily]*" (mW/m"^"2"*"/sr/nm)")
y_limit_gpp  <- c(0,16)
y_limit_sif  <- c(0,1)

par(mfrow = c(3, 3), oma=c(1,3.5,2.5,3.5))

# for (i in 1:length(towers$sites)) {
for (i in 1:9) {
  
  site_name <- towers$site[i]
  
  tower_data <- get_tower_gpp(towers$data[i])[[1]]
  
  # Get model GPP and SIF data at tower
  tower_loc <- vect(cbind(towers$lon[i], towers$lat[i]), crs = "+proj=longlat +ellps=WGS84")
  
  cf_site   <- extract(cf, tower_loc)
  fcom_site <- extract(fcom, tower_loc)
  fsat_site <- extract(fsat, tower_loc)
  sif_site  <- extract(sif, tower_loc)
  
  cf_site   <- unlist(cf_site, use.names = FALSE)[-c(1)]
  fcom_site <- unlist(fcom_site, use.names = FALSE)[-c(1)]
  fsat_site <- unlist(fsat_site, use.names = FALSE)[-c(1)]
  sif_site  <- unlist(sif_site, use.names = FALSE)[-c(1)]
  
  cf_site_month   <- to_month(cf_site, 2019)
  fcom_site_month <- to_month(fcom_site, 2019)
  fsat_site_month <- to_month(fsat_site, 2019)
  sif_site_month  <- to_month(sif_site, 2019)
  
  # Do not plot if SIF or ChloFluo is missing
  if (all(!is.na(sif_site_month)) || all(!is.na(cf_site_month)) ||
      all(!is.nan(sif_site_month)) || all(!is.nan(cf_site_month))) {
    ## PLOTTING
    
    op <- par(mar = c(2.5,0,0,0))
    
    plot(x, tower_data, col = "black", type = "l", axes = FALSE, lwd = 2, lty = 2, xaxs="i", ylim = y_limit_gpp, ylab = "")
    
    lines(x, fcom_site_month, col = "#FE6100", lwd = 2)
    lines(x, fsat_site_month, col = "#648FFF", lwd = 2)
    lines(x, cf_site_month, col = "#DC267F", lwd = 2)
    
    if (i == 7 || i == 8 || i == 9) {
      axis(1, tck = 0.03, labels = x_lab, at = x, mgp=c(3, 0.2, 0))
    } else {
      axis(1, tck = 0.03, labels = FALSE, at = x, mgp=c(3, 0.2, 0))
    }
    
    if (i == 1 || i == 4 || i == 7) {
      axis(2, tck = 0.03, labels = TRUE, mgp=c(3, 0.2, 0), las = 2)
    } else {
      axis(2, tck = 0.03, labels = FALSE, mgp=c(3, 0.2, 0), las = 2)
    }
    
    par(new = TRUE)
    plot(x, sif_site_month, col = "gray50", type = "l", axes = FALSE, lwd = 2, lty = 2, xaxs="i", ylim = y_limit_sif, ylab = "")
    
    if (i == 3 || i == 6 || i == 9) {
      axis(4, tck = 0.03, labels = TRUE, mgp=c(3, 0.2, 0), las = 2)
    } else {
      axis(4, tck = 0.03, labels = FALSE, mgp=c(3, 0.2, 0), las = 2)
    }
    
    
    mtext(3, text = site_name, line = 0)
    
    if (i == 1) {
      legend("topleft", legend=c("Tower", "ChloFluo", "FluxCom", "FluxSat", "SIF"),
             col=c("black", "#DC267F", "#FE6100", "#648FFF", "gray50"), lty = c(2,1,1,1,2),
             lwd = 1, ncol = 2, y.intersp = 1, cex = 0.75, bty = "n", bg =)
    }
    
    # Regressions
    cf_reg   <- summary(lm(cf_site_month ~ tower_data))
    fcom_reg <- summary(lm(fcom_site_month ~ tower_data))
    fsat_reg <- summary(lm(fsat_site_month ~ tower_data))
    sif_reg  <- summary(lm(sif_site_month ~ tower_data))
    
    cf_r   <- round2(cf_reg$adj.r.squared, 2, FALSE)
    fcom_r <- round2(fcom_reg$adj.r.squared, 2, FALSE)
    fsat_r <- round2(fsat_reg$adj.r.squared, 2, FALSE)
    sif_r  <- round2(sif_reg$adj.r.squared, 2, FALSE)
    
    cf_p   <- round2(cf_reg$coefficients[2,4], 2, TRUE)
    fcom_p <- round2(fcom_reg$coefficients[2,4], 2, TRUE)
    fsat_p <- round2(fsat_reg$coefficients[2,4], 2, TRUE)
    sif_p  <- round2(sif_reg$coefficients[2,4], 2, TRUE)
    
    cf_rmse   <- round2(sqrt(mean(cf_reg$residuals^2)), 2, FALSE)
    fcom_rmse <- round2(sqrt(mean(fcom_reg$residuals^2)), 2, FALSE)
    fsat_rmse <- round2(sqrt(mean(fsat_reg$residuals^2)), 2, FALSE)
    sif_rmse  <- round2(sqrt(mean(sif_reg$residuals^2)), 2, FALSE)
    
    if (i == 4) { # Cor is insignificant for K34 for fcom and fsat
      cf_report   <- paste0("R2 = ", cf_r, "; RMSE = ", cf_rmse)
      fcom_report <- paste0("R2 = ", fcom_r, "*; RMSE = ", fcom_rmse)
      fsat_report <- paste0("R2 = ", fsat_r, "*; RMSE = ", fsat_rmse)
      sif_report  <- paste0("R2 = ", sif_r, "; RMSE = ", sif_rmse)
    } else {
      cf_report   <- paste0("R2 = ", cf_r, "; RMSE = ", cf_rmse)
      fcom_report <- paste0("R2 = ", fcom_r, "; RMSE = ", fcom_rmse)
      fsat_report <- paste0("R2 = ", fsat_r, "; RMSE = ", fsat_rmse)
      sif_report  <- paste0("R2 = ", sif_r, "; RMSE = ", sif_rmse)
    }
    
    legend("topright", legend=c(cf_report, fcom_report, fsat_report, sif_report),
           text.col=c("#DC267F", "#FE6100", "#648FFF", "gray50"), lty = NA,
           lwd = NA, ncol = 1, y.intersp = 1, cex = 0.75, bty = "n", bg = "white")
    
    box()
  }
    
  }

mtext(1, text = "Month", outer = TRUE, line = -0.5)
mtext(2, text = as.expression(y_lab_gpp), outer = TRUE, line = 1.5)
mtext(4, text = as.expression(y_lab_sif), outer = TRUE, line = 2.5)



# 
# ## Plot variables
# x            <- 1:12
# x_lab        <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
# y_lab_gpp    <- bquote("GPP (g C m"^"-2"*" d"^"-1"*")")
# y_lab_sif    <- bquote("SIF"[Daily]*" (mW/m"^"2"*"/sr/nm)")
# y_limit_gpp  <- c(0,16)
# y_limit_sif  <- c(0,1)
# 
# ### PLOT
# cairo_pdf(out_name, width = 7.5, height = 6.25)
# 
# par(mfrow = c(3, 3), oma=c(1,3.5,2.5,3.5))
# 
# for (i in 1:9){
#   
#   op <- par(mar = c(2.5,0,0,0))
#   
#   plot(x, df_tower[,i], col = "black", type = "l", axes = FALSE, lwd = 2, lty = 2, xaxs="i", ylim = y_limit_gpp, ylab = "")
#   
#   lines(x, fcom_site_month, col = "#FE6100", lwd = 2)
#   lines(x, df_fsat[,i], col = "#648FFF", lwd = 2)
#   lines(x, df_cf[,i], col = "#DC267F", lwd = 2)
#   
#   if (i == 7 || i == 8 || i == 9) {
#     axis(1, tck = 0.03, labels = x_lab, at = x, mgp=c(3, 0.2, 0))
#   } else {
#     axis(1, tck = 0.03, labels = FALSE, at = x, mgp=c(3, 0.2, 0))
#   }
#   
#   if (i == 1 || i == 4 || i == 7) {
#     axis(2, tck = 0.03, labels = TRUE, mgp=c(3, 0.2, 0), las = 2)
#   } else {
#     axis(2, tck = 0.03, labels = FALSE, mgp=c(3, 0.2, 0), las = 2)
#   }
#   
#   par(new = TRUE)
#   plot(x, sif_site_month, col = "gray50", type = "l", axes = FALSE, lwd = 2, lty = 2, xaxs="i", ylim = y_limit_sif, ylab = "")
#   
#   if (i == 3 || i == 6 || i == 9) {
#     axis(4, tck = 0.03, labels = TRUE, mgp=c(3, 0.2, 0), las = 2)
#   } else {
#     axis(4, tck = 0.03, labels = FALSE, mgp=c(3, 0.2, 0), las = 2)
#   }
# 
#   
#   mtext(3, text = tower_names[i], line = 0)
#   
#   if (i == 1) {
#       legend("topleft", legend=c("Tower", "ChloFluo", "FluxCom", "FluxSat", "SIF"),
#          col=c("black", "#DC267F", "#FE6100", "#648FFF", "gray50"), lty = c(2,1,1,1,2),
#          lwd = 1, ncol = 2, y.intersp = 1, cex = 0.75, bty = "n", bg =)
#   }
#   
#   # Regressions
#   cf_reg   <- summary(lm(df_cf[,i] ~ df_tower[,i]))
#   fcom_reg <- summary(lm(fcom_site_month ~ df_tower[,i]))
#   fsat_reg <- summary(lm(df_fsat[,i] ~ df_tower[,i]))
#   sif_reg  <- summary(lm(sif_site_month ~ df_tower[,i]))
# 
#   cf_r   <- round2(cf_reg$adj.r.squared, 2, FALSE)
#   fcom_r <- round2(fcom_reg$adj.r.squared, 2, FALSE)
#   fsat_r <- round2(fsat_reg$adj.r.squared, 2, FALSE)
#   sif_r  <- round2(sif_reg$adj.r.squared, 2, FALSE)
#   
#   cf_p   <- round2(cf_reg$coefficients[2,4], 2, TRUE)
#   fcom_p <- round2(fcom_reg$coefficients[2,4], 2, TRUE)
#   fsat_p <- round2(fsat_reg$coefficients[2,4], 2, TRUE)
#   sif_p  <- round2(sif_reg$coefficients[2,4], 2, TRUE)
#   
#   cf_rmse   <- round2(sqrt(mean(cf_reg$residuals^2)), 2, FALSE)
#   fcom_rmse <- round2(sqrt(mean(fcom_reg$residuals^2)), 2, FALSE)
#   fsat_rmse <- round2(sqrt(mean(fsat_reg$residuals^2)), 2, FALSE)
#   sif_rmse  <- round2(sqrt(mean(sif_reg$residuals^2)), 2, FALSE)
#   
#   if (i == 4) { # Cor is insignificant for K34 for fcom and fsat
#     cf_report   <- paste0("R2 = ", cf_r, "; RMSE = ", cf_rmse)
#     fcom_report <- paste0("R2 = ", fcom_r, "*; RMSE = ", fcom_rmse)
#     fsat_report <- paste0("R2 = ", fsat_r, "*; RMSE = ", fsat_rmse)
#     sif_report  <- paste0("R2 = ", sif_r, "; RMSE = ", sif_rmse)
#   } else {
#     cf_report   <- paste0("R2 = ", cf_r, "; RMSE = ", cf_rmse)
#     fcom_report <- paste0("R2 = ", fcom_r, "; RMSE = ", fcom_rmse)
#     fsat_report <- paste0("R2 = ", fsat_r, "; RMSE = ", fsat_rmse)
#     sif_report  <- paste0("R2 = ", sif_r, "; RMSE = ", sif_rmse)
#   }
#   
#   legend("topright", legend=c(cf_report, fcom_report, fsat_report, sif_report),
#          text.col=c("#DC267F", "#FE6100", "#648FFF", "gray50"), lty = NA,
#          lwd = NA, ncol = 1, y.intersp = 1, cex = 0.75, bty = "n", bg = "white")
#   
#   box()
#   
# }
# 
# mtext(1, text = "Month", outer = TRUE, line = -0.5)
# mtext(2, text = as.expression(y_lab_gpp), outer = TRUE, line = 1.5)
# mtext(4, text = as.expression(y_lab_sif), outer = TRUE, line = 2.5)
# 
# dev.off()