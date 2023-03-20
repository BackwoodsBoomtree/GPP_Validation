library(terra)
library(viridis)
library(rgdal)
library(RColorBrewer)

out_name    <- "G:/ChloFluo/comps/towers/gpp/tower_comparisons_all"

cf_file     <- "G:/ChloFluo/product/v02/clipfill/ChloFluo.GPP.v02.1deg.CF80.2019.clipfill.nc"
fcom_file   <- "G:/FluxCom/RS/GPP.RS_V006.FP-ALL.MLM-ALL.METEO-NONE.720_360.8daily.2019.nc"
fsat_files  <- "G:/FluxSat/monthly/1deg/GPP_FluxSat_8day_1deg_v2_2019.nc"
sif_file    <- "G:/ChloFluo/input/SIF/1deg/clipfill/SIFqc.8day.1deg.veg_only.CF80.2019.clipfill.nc"
tower_list  <- "G:/ChloFluo/comps/tower-data/Joiner_2020_Sites.csv"
tower_dir   <- "G:/ChloFluo/comps/tower-data/unzipped"
c3_c4_file  <- "G:/ChloFluo/input/C3C4/ISLSCP/c4_percent_1d.nc"

kick_out_sites <- c("AU-ASM", "AU-Gin", "ES-LgS", "GF-Guy", "IT-Cp2", "IT-Cpz",
                    "NL-Hor", "NL-Loo", "RU-Cok", "US-KS2")

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
  df <- read.csv(f, header = TRUE)
  
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

extract_tower_var <- function(f, loc) {
  site_data <- extract(f, loc)
  site_data <- unlist(site_data, use.names = FALSE)[-c(1)]
  site_data[is.nan(site_data)] <- NA
  return(site_data)
}

cf     <- rast(cf_file, subds = "gpp")
fcom   <- rast(fcom_file, subds = "GPP")
fsat   <- rast(fsat_files, subds = "GPP")
sif    <- rast(sif_file, subds = "sif743_qc")
fcom   <- aggregate(fcom, 2, fun = mean, na.rm = TRUE)
towers <- read.csv(tower_list, header = TRUE)
c3_c4  <- rast(c3_c4_file, subds = "Band1")

### Get csv files for each tower after filtering
for (i in 1:length(kick_out_sites)) {
  towers <- subset(towers, site != kick_out_sites[i])
}

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
towers           <- na.omit(towers)
rownames(towers) <- 1:nrow(towers)

tower_df_len <- length(towers$site)


## Plot variables
x            <- 1:12
x_lab        <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
y_lab_gpp    <- bquote("GPP (g C m"^"-2"*" d"^"-1"*")")
y_lab_sif    <- bquote("SIF"[Daily]*" (mW/m"^"2"*"/sr/nm)")
y_limit_gpp  <- c(0,24)
y_limit_sif  <- c(0,1)

run_n      <- (round(length(towers$site) / 9))
plot_index <- seq(0, (run_n * 9), by = 9)

for (j in 1:run_n) {
  
  # out_file <- paste0(out_name, "_", sprintf("%02d", j), ".pdf")
  # cairo_pdf(out_file, width = 8, height = 6.25)
  out_file <- paste0(out_name, "_", sprintf("%02d", j), ".svg")
  svg(out_file, width = 8, height = 6.25)
  
  par(mfrow = c(3, 3), oma=c(1,3,2.5,3))
  
  # for (i in 1:length(towers$sites)) {
  for (i in 1:9) {
    i_adj <- i + plot_index[j]
    
    if (i_adj <= tower_df_len) {
      site_name  <- paste0(towers$site[i_adj], " (", towers$veg[i_adj], ")")
      tower_file <- towers$data[i_adj]
      
      tower_all_data <- get_tower_gpp(tower_file)
      tower_gpp_data <- tower_all_data[[1]]
      
      # Get model GPP and SIF data at tower
      tower_loc <- vect(cbind(towers$lon[i_adj], towers$lat[i_adj]), crs = "+proj=longlat +ellps=WGS84")
      
      cf_site    <- extract_tower_var(cf, tower_loc)
      fcom_site  <- extract_tower_var(fcom, tower_loc)
      fsat_site  <- extract_tower_var(fsat, tower_loc)
      sif_site   <- extract_tower_var(sif, tower_loc)
      c3_c4_site <- extract(c3_c4, tower_loc)[[2]]
      c3_c4_report <- paste0("C4: ", c3_c4_site, "%")
      
      cf_site_month   <- to_month(cf_site, 2019)
      fcom_site_month <- to_month(fcom_site, 2019)
      fsat_site_month <- to_month(fsat_site, 2019)
      sif_site_month  <- to_month(sif_site, 2019)
      
      ## PLOTTING
      
      op <- par(mar = c(2.5,0.5,0,0.5))
      
      plot(x, tower_gpp_data, col = "black", type = "l", axes = FALSE, lwd = 2, lty = 2, xaxs="i", ylim = y_limit_gpp, ylab = "")
      
      # Plot and run regression only if there is data. Add results to tower df
      if (all(is.na(fcom_site_month)) == FALSE) {
        lines(x, fcom_site_month, col = "#FE6100", lwd = 2)
        fcom_reg    <- summary(lm(fcom_site_month ~ tower_gpp_data))
        fcom_r      <- round2(fcom_reg$adj.r.squared, 2, FALSE)
        fcom_p      <- round2(fcom_reg$coefficients[2,4], 2, FALSE)
        fcom_rmse   <- round2(sqrt(mean(fcom_reg$residuals^2)), 2, FALSE)
        towers$fcom_r[i_adj]    <- fcom_r
        towers$fcom_p[i_adj]    <- fcom_p
        towers$fcom_rmse[i_adj] <- fcom_rmse
        if (fcom_p < 0.05) {
          fcom_report <- paste0("R2 = ", fcom_r, "; RMSE = ", fcom_rmse)
        } else {
          fcom_report <- paste0("*R2 = ", fcom_r, "; RMSE = ", fcom_rmse)
        }
      } else {
        towers$fcom_r[i_adj]    <- NA
        towers$fcom_p[i_adj]    <- NA
        towers$fcom_rmse[i_adj] <- NA
        fsat_report <- NA
      }
      if (all(is.na(fsat_site_month)) == FALSE) {
        lines(x, fsat_site_month, col = "#648FFF", lwd = 2)
        fsat_reg    <- summary(lm(fsat_site_month ~ tower_gpp_data))
        fsat_r      <- round2(fsat_reg$adj.r.squared, 2, FALSE)
        fsat_p      <- round2(fsat_reg$coefficients[2,4], 2, FALSE)
        fsat_rmse   <- round2(sqrt(mean(fsat_reg$residuals^2)), 2, FALSE)
        towers$fsat_r[i_adj]    <- fsat_r
        towers$fsat_p[i_adj]    <- fsat_p
        towers$fsat_rmse[i_adj] <- fsat_rmse
        if (fsat_p < 0.05) {
          fsat_report <- paste0("R2 = ", fsat_r, "; RMSE = ", fsat_rmse)
        } else {
          fsat_report <- paste0("*R2 = ", fsat_r, "; RMSE = ", fsat_rmse)
        }
      } else {
        towers$fsat_r[i_adj]    <- NA
        towers$fsat_p[i_adj]    <- NA
        towers$fsat_rmse[i_adj] <- NA
        fsat_report <- NA
      }
      if (all(is.na(cf_site_month)) == FALSE) {
        lines(x, cf_site_month, col = "#DC267F", lwd = 2)
        cf_reg    <- summary(lm(cf_site_month ~ tower_gpp_data))
        cf_r      <- round2(cf_reg$adj.r.squared, 2, FALSE)
        cf_p      <- round2(cf_reg$coefficients[2,4], 2, FALSE)
        cf_rmse   <- round2(sqrt(mean(cf_reg$residuals^2)), 2, FALSE)
        towers$cf_r[i_adj]    <- cf_r
        towers$cf_p[i_adj]    <- cf_p
        towers$cf_rmse[i_adj] <- cf_rmse
        if (cf_p < 0.05) {
          cf_report <- paste0("R2 = ", cf_r, "; RMSE = ", cf_rmse)
        } else {
          cf_report <- paste0("*R2 = ", cf_r, "; RMSE = ", cf_rmse)
        }
      } else {
        towers$cf_r[i_adj]    <- NA
        towers$cf_p[i_adj]    <- NA
        towers$cf_rmse[i_adj] <- NA
        cf_report <- NA
      }
      
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
      
      # SIF
      if (all(!is.na(sif_site_month)) && all(!is.nan(sif_site_month))) {
        par(new = TRUE)
        plot(x, sif_site_month, col = "gray50", type = "l", axes = FALSE, lwd = 2, lty = 2, xaxs="i", ylim = y_limit_sif, ylab = "")
        
        sif_reg    <- summary(lm(sif_site_month ~ tower_gpp_data))
        sif_r      <- round2(sif_reg$adj.r.squared, 2, FALSE)
        sif_p      <- round2(sif_reg$coefficients[2,4], 2, FALSE)
        sif_rmse   <- round2(sqrt(mean(sif_reg$residuals^2)), 2, FALSE)
        towers$sif_r[i_adj]    <- sif_r
        towers$sif_p[i_adj]    <- sif_p
        towers$sif_rmse[i_adj] <- sif_rmse
        if (sif_p < 0.05) {
          sif_report <- paste0("R2 = ", sif_r, "; RMSE = ", sif_rmse)
        } else {
          sif_report <- paste0("*R2 = ", sif_r, "; RMSE = ", sif_rmse)
        }
        
        if (i == 3 || i == 6 || i == 9) {
          axis(4, tck = 0.03, labels = TRUE, mgp=c(3, 0.2, 0), las = 2)
        } else {
          axis(4, tck = 0.03, labels = FALSE, mgp=c(3, 0.2, 0), las = 2)
        }
      } else {
        towers$sif_r[i_adj]    <- NA
        towers$sif_p[i_adj]    <- NA
        towers$sif_rmse[i_adj] <- NA
        sif_report <- NA
      }
      
      mtext(3, text = site_name, line = 0)
      
      if (i == 1) {
        legend("topleft", legend = c("Tower", "ChloFluo", "FluxCom", "FluxSat", "SIF"),
               col=c("black", "#DC267F", "#FE6100", "#648FFF", "gray50"), lty = c(2,1,1,1,2),
               lwd = 1, ncol = 2, y.intersp = 1, cex = 0.75, bty = "n", bg =)
      }
      
      legend("topright", legend = c(cf_report, fcom_report, fsat_report, sif_report, c3_c4_report),
             text.col=c("#DC267F", "#FE6100", "#648FFF", "gray50", "black"), lty = NA,
             lwd = NA, ncol = 1, y.intersp = 1, cex = 0.75, bty = "n", bg = "white")
      
      box()
      
    } else if (i_adj == (tower_df_len + 1)) {
      k34_details           <- t(data.frame(c("BR-K34", -2.6091, -60.2093, 1999, 2006, "EBF", "Wu et al. (2015)", 
                                              "G:/SIF_comps/figs/Wu_2016/K34_GEP.csv", rep(NA, 12))))
      colnames(k34_details) <- colnames(towers)
      
      towers           <- rbind(towers, k34_details)
      rownames(towers) <- seq(1, nrow(towers))
      
      k34_gpp <- read.csv("G:/SIF_comps/figs/Wu_2016/K34_GEP.csv", header = FALSE)[,2]
      k34     <- vect(cbind(-60.2093, -2.6091), crs = "+proj=longlat +ellps=WGS84")
      
      cf_k34     <- extract(cf, k34)
      fcom_k34   <- extract(fcom, k34)
      fsat_k34   <- extract(fsat, k34)
      sif_k34    <- extract(sif, k34)
      c3_c4_site <- extract(c3_c4, k34)
      c3_c4_report <- paste0("C4: ", c3_c4_site, "%")[[2]]
      
      cf_k34   <- unlist(cf_k34, use.names = FALSE)[-c(1)]
      fcom_k34 <- unlist(fcom_k34, use.names = FALSE)[-c(1)]
      fsat_k34 <- unlist(fsat_k34, use.names = FALSE)[-c(1)]
      sif_k34  <- unlist(sif_k34, use.names = FALSE)[-c(1)]
      
      cf_site_month   <- to_month(cf_k34, 2019)
      fcom_site_month <- to_month(fcom_k34, 2019)
      fsat_site_month <- to_month(fsat_k34, 2019)
      sif_site_month  <- to_month(sif_k34, 2019)
      
      plot(x, k34_gpp, col = "black", type = "l", axes = FALSE, lwd = 2, lty = 2, xaxs="i", ylim = y_limit_gpp, ylab = "")
      
      # Plot and run regression only if there is data. Add results to tower df
      if (all(!is.na(fcom_site_month)) && all(!is.nan(fcom_site_month))) {
        lines(x, fcom_site_month, col = "#FE6100", lwd = 2)
        fcom_reg    <- summary(lm(fcom_site_month ~ k34_gpp))
        fcom_r      <- round2(fcom_reg$adj.r.squared, 2, FALSE)
        fcom_p      <- round2(fcom_reg$coefficients[2,4], 2, FALSE)
        fcom_rmse   <- round2(sqrt(mean(fcom_reg$residuals^2)), 2, FALSE)
        towers$fcom_r[i_adj]    <- fcom_r
        towers$fcom_p[i_adj]    <- fcom_p
        towers$fcom_rmse[i_adj] <- fcom_rmse
        if (fcom_p < 0.05) {
          fcom_report <- paste0("R2 = ", fcom_r, "; RMSE = ", fcom_rmse)
        } else {
          fcom_report <- paste0("*R2 = ", fcom_r, "; RMSE = ", fcom_rmse)
        }
      } else {
        fsat_report <- NA
      }
      if (all(!is.na(fsat_site_month)) && all(!is.nan(fsat_site_month))) {
        lines(x, fsat_site_month, col = "#648FFF", lwd = 2)
        fsat_reg    <- summary(lm(fsat_site_month ~ k34_gpp))
        fsat_r      <- round2(fsat_reg$adj.r.squared, 2, FALSE)
        fsat_p      <- round2(fsat_reg$coefficients[2,4], 2, FALSE)
        fsat_rmse   <- round2(sqrt(mean(fsat_reg$residuals^2)), 2, FALSE)
        towers$fsat_r[i_adj]    <- fsat_r
        towers$fsat_p[i_adj]    <- fsat_p
        towers$fsat_rmse[i_adj] <- fsat_rmse
        if (fsat_p < 0.05) {
          fsat_report <- paste0("R2 = ", fsat_r, "; RMSE = ", fsat_rmse)
        } else {
          fsat_report <- paste0("*R2 = ", fsat_r, "; RMSE = ", fsat_rmse)
        }
      } else {
        fsat_report <- NA
      }
      if (all(!is.na(cf_site_month)) && all(!is.nan(cf_site_month))) {
        lines(x, cf_site_month, col = "#DC267F", lwd = 2)
        cf_reg    <- summary(lm(cf_site_month ~ k34_gpp))
        cf_r      <- round2(cf_reg$adj.r.squared, 2, FALSE)
        cf_p      <- round2(cf_reg$coefficients[2,4], 2, FALSE)
        cf_rmse   <- round2(sqrt(mean(cf_reg$residuals^2)), 2, FALSE)
        towers$cf_r[i_adj]    <- cf_r
        towers$cf_p[i_adj]    <- cf_p
        towers$cf_rmse[i_adj] <- cf_rmse
        if (cf_p < 0.05) {
          cf_report <- paste0("R2 = ", cf_r, "; RMSE = ", cf_rmse)
        } else {
          cf_report <- paste0("*R2 = ", cf_r, "; RMSE = ", cf_rmse)
        }
      } else {
        cf_report <- NA
      }
      
      # SIF
      if (all(!is.na(sif_site_month)) && all(!is.nan(sif_site_month))) {
        par(new = TRUE)
        plot(x, sif_site_month, col = "gray50", type = "l", axes = FALSE, lwd = 2, lty = 2, xaxs="i", ylim = y_limit_sif, ylab = "")
        
        sif_reg    <- summary(lm(sif_site_month ~ k34_gpp))
        sif_r      <- round2(sif_reg$adj.r.squared, 2, FALSE)
        sif_p      <- round2(sif_reg$coefficients[2,4], 2, FALSE)
        sif_rmse   <- round2(sqrt(mean(sif_reg$residuals^2)), 2, FALSE)
        towers$sif_r[i_adj]    <- sif_r
        towers$sif_p[i_adj]    <- sif_p
        towers$sif_rmse[i_adj] <- sif_rmse
        if (sif_p < 0.05) {
          sif_report <- paste0("R2 = ", sif_r, "; RMSE = ", sif_rmse)
        } else {
          sif_report <- paste0("*R2 = ", sif_r, "; RMSE = ", sif_rmse)
        } 
      }
      
      axis(4, tck = 0.03, labels = TRUE, mgp=c(3, 0.2, 0), las = 2)
      axis(1, tck = 0.03, labels = x_lab, at = x, mgp=c(3, 0.2, 0))
      
      mtext(3, text = "BR-K34 (EBF)", line = 0)
      
      legend("topright", legend = c(cf_report, fcom_report, fsat_report, sif_report, c3_c4_report),
             text.col=c("#DC267F", "#FE6100", "#648FFF", "gray50", "black"), lty = NA,
             lwd = NA, ncol = 1, y.intersp = 1, cex = 0.75, bty = "n", bg = "white")

      box()
      
    } else {
      #nothing
    }
  }
  
  mtext(1, text = "Month", outer = TRUE, line = -0.5)
  mtext(2, text = as.expression(y_lab_gpp), outer = TRUE, line = 1)
  mtext(4, text = as.expression(y_lab_sif), outer = TRUE, line = 2)
  
  dev.off()
}


towers_qc <- na.omit(towers)
hist(as.numeric(towers_qc$cf_r))
hist(as.numeric(towers_qc$fcom_r))
hist(as.numeric(towers_qc$fsat_r))
hist(as.numeric(towers_qc$sif_r))

hist(as.numeric(towers_qc$cf_rmse))
hist(as.numeric(towers_qc$fcom_rmse))
hist(as.numeric(towers_qc$fsat_rmse))
hist(as.numeric(towers_qc$sif_rmse))
