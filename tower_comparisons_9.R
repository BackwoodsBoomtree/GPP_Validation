library(terra)
library(viridis)
library(rgdal)
library(RColorBrewer)

out_name <- "G:/ChloFluo/comps/tower_comparisons_9_v3.pdf"

cf_file    <- "G:/ChloFluo/product/v02/clipfill/ChloFluo.GPP.v02.1deg.CF80.2019.clipfill.nc"
fcom_file  <- "G:/FluxCom/RS/GPP.RS_V006.FP-ALL.MLM-ALL.METEO-NONE.720_360.8daily.2019.nc"
fsat_files <- "G:/FluxSat/monthly/1deg/GPP_FluxSat_8day_1deg_v2_2019.nc"
sif_file   <- "G:/ChloFluo/input/SIF/1deg/SIFqc.8day.1deg.veg_only.CF80.2019.nc"
tower_locs <- "G:/ChloFluo/comps/tower-data/Sites_Chosen.csv"

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


cf     <- rast(cf_file, subds = "gpp")
fcom   <- rast(fcom_file, subds = "GPP")
fsat   <- rast(fsat_files, subds = "GPP")
sif    <- rast(sif_file, subds = "sif_qc")
fcom   <- aggregate(fcom, 2, fun = mean, na.rm = TRUE)
towers <- read.csv(tower_locs, header = TRUE)

# Get Tower Data
k34_gpp    <- read.csv("G:/SIF_comps/figs/Wu_2016/K34_GEP.csv", header = FALSE)[,2]
ade_gpp    <- read.csv("G:/ChloFluo/comps/tower-data/FLX_AU-Ade/FLX_AU-Ade_FLUXNET2015_FULLSET_MM_2007-2009_1-4.csv", header = TRUE)
das_gpp    <- read.csv("G:/ChloFluo/comps/tower-data/FLX_AU-DaS/FLX_AU-DaS_FLUXNET2015_FULLSET_MM_2008-2014_2-4.csv", header = TRUE)
dry_gpp    <- read.csv("G:/ChloFluo/comps/tower-data/FLX_AU-Dry/FLX_AU-Dry_FLUXNET2015_FULLSET_MM_2008-2014_2-4.csv", header = TRUE)
lom_gpp    <- read.csv("G:/ChloFluo/comps/tower-data/FLX_FI-Lom/FLX_FI-Lom_FLUXNET2015_FULLSET_MM_2007-2009_1-4.csv", header = TRUE)
fyo_gpp    <- read.csv("G:/ChloFluo/comps/tower-data/FLX_RU-Fyo/FLX_RU-Fyo_FLUXNET2015_FULLSET_MM_1998-2014_2-4.csv", header = TRUE)
ar1_gpp    <- read.csv("G:/ChloFluo/comps/tower-data/FLX_US-AR1/FLX_US-AR1_FLUXNET2015_FULLSET_MM_2009-2012_1-4.csv", header = TRUE)
ar2_gpp    <- read.csv("G:/ChloFluo/comps/tower-data/FLX_US-AR2/FLX_US-AR2_FLUXNET2015_FULLSET_MM_2009-2012_1-4.csv", header = TRUE)
mon_gpp    <- read.csv("G:/ChloFluo/comps/tower-data/FLX_ZM-Mon/FLX_ZM-Mon_FLUXNET2015_FULLSET_MM_2000-2009_2-4.csv", header = TRUE)

ade_gpp_dt    <- ade_gpp$GPP_NT_VUT_MEAN
das_gpp_dt    <- das_gpp$GPP_NT_VUT_MEAN
dry_gpp_dt    <- dry_gpp$GPP_NT_VUT_MEAN
lom_gpp_dt    <- lom_gpp$GPP_NT_VUT_MEAN
fyo_gpp_dt    <- fyo_gpp$GPP_NT_VUT_MEAN
ar1_gpp_dt    <- ar1_gpp$GPP_NT_VUT_MEAN
ar2_gpp_dt    <- ar2_gpp$GPP_NT_VUT_MEAN
mon_gpp_dt    <- mon_gpp$GPP_NT_VUT_MEAN

ade_gpp_nt    <- ade_gpp$GPP_DT_VUT_MEAN
das_gpp_nt    <- das_gpp$GPP_DT_VUT_MEAN
dry_gpp_nt    <- dry_gpp$GPP_DT_VUT_MEAN
lom_gpp_nt    <- lom_gpp$GPP_DT_VUT_MEAN
fyo_gpp_nt    <- fyo_gpp$GPP_DT_VUT_MEAN
ar1_gpp_nt    <- ar1_gpp$GPP_DT_VUT_MEAN
ar2_gpp_nt    <- ar2_gpp$GPP_DT_VUT_MEAN
mon_gpp_nt    <- mon_gpp$GPP_DT_VUT_MEAN

ade_gpp_mean    <- (ade_gpp_nt + ade_gpp_dt) / 2
das_gpp_mean    <- (das_gpp_nt + das_gpp_dt) / 2
dry_gpp_mean    <- (dry_gpp_nt + dry_gpp_dt) / 2
lom_gpp_mean    <- (lom_gpp_nt + lom_gpp_dt) / 2
fyo_gpp_mean    <- (fyo_gpp_nt + fyo_gpp_dt) / 2
ar1_gpp_mean    <- (ar1_gpp_nt + ar1_gpp_dt) / 2
ar2_gpp_mean    <- (ar2_gpp_nt + ar2_gpp_dt) / 2
mon_gpp_mean    <- (mon_gpp_nt + mon_gpp_dt) / 2

# Get annual monthly mean
ade_gpp_mean <- (ade_gpp_mean[1:12] + ade_gpp_mean[13:24] + ade_gpp_mean[25:36]) / 3
das_gpp_mean <- (das_gpp_mean[1:12] + das_gpp_mean[13:24] + das_gpp_mean[25:36] +
                 das_gpp_mean[37:48] + das_gpp_mean[49:60] + das_gpp_mean[61:72] +
                 das_gpp_mean[73:84]) / 7
dry_gpp_mean <- (dry_gpp_mean[1:12] + dry_gpp_mean[13:24] + dry_gpp_mean[25:36] +
                 dry_gpp_mean[37:48] + dry_gpp_mean[49:60] + dry_gpp_mean[61:72] +
                 dry_gpp_mean[73:84]) / 7
lom_gpp_mean <- (lom_gpp_mean[1:12] + lom_gpp_mean[13:24] + lom_gpp_mean[25:36]) / 3
fyo_gpp_mean <- (fyo_gpp_mean[1:12] + fyo_gpp_mean[13:24] + fyo_gpp_mean[25:36] +
                 fyo_gpp_mean[37:48] + fyo_gpp_mean[49:60] + fyo_gpp_mean[61:72] +
                 fyo_gpp_mean[73:84] + fyo_gpp_mean[85:96] + fyo_gpp_mean[97:108] +
                 fyo_gpp_mean[109:120] + fyo_gpp_mean[121:132] + fyo_gpp_mean[133:144] +
                 fyo_gpp_mean[145:156] + fyo_gpp_mean[157:168] + fyo_gpp_mean[169:180] +
                 fyo_gpp_mean[181:192] + fyo_gpp_mean[193:204]) / 17
ar1_gpp_mean <- (ar1_gpp_mean[1:12] + ar1_gpp_mean[13:24] + ar1_gpp_mean[25:36] +
                 ar1_gpp_mean[37:48]) / 4
ar2_gpp_mean <- (ar2_gpp_mean[1:12] + ar2_gpp_mean[13:24] + ar2_gpp_mean[25:36] +
                 ar2_gpp_mean[37:48]) / 4
mon_gpp_mean <- (mon_gpp_mean[1:12] + mon_gpp_mean[85:96] + mon_gpp_mean[97:108] +
                 mon_gpp_mean[109:120]) / 4

## Tower vectors
k34 <- vect(cbind(-60.2093, -2.6091), crs = "+proj=longlat +ellps=WGS84")
ade <- vect(cbind(towers$lon[1], towers$lat[1]), crs = "+proj=longlat +ellps=WGS84")
das <- vect(cbind(towers$lon[2], towers$lat[2]), crs = "+proj=longlat +ellps=WGS84")
dry <- vect(cbind(towers$lon[3], towers$lat[3]), crs = "+proj=longlat +ellps=WGS84")
lom <- vect(cbind(towers$lon[4], towers$lat[4]), crs = "+proj=longlat +ellps=WGS84")
fyo <- vect(cbind(towers$lon[5], towers$lat[5]), crs = "+proj=longlat +ellps=WGS84")
ar1 <- vect(cbind(towers$lon[6], towers$lat[6]), crs = "+proj=longlat +ellps=WGS84")
ar2 <- vect(cbind(towers$lon[7], towers$lat[7]), crs = "+proj=longlat +ellps=WGS84")
mon <- vect(cbind(towers$lon[8], towers$lat[8]), crs = "+proj=longlat +ellps=WGS84")

# Timeseries:
cf_k34   <- extract(cf, k34)
fcom_k34 <- extract(fcom, k34)
fsat_k34 <- extract(fsat, k34)
sif_k34  <- extract(sif, k34)
cf_k34   <- unlist(cf_k34, use.names = FALSE)[-c(1)]
fcom_k34 <- unlist(fcom_k34, use.names = FALSE)[-c(1)]
fsat_k34 <- unlist(fsat_k34, use.names = FALSE)[-c(1)]
sif_k34  <- unlist(sif_k34, use.names = FALSE)[-c(1)]

cf_ade   <- extract(cf, ade)
fcom_ade <- extract(fcom, ade)
fsat_ade <- extract(fsat, ade)
sif_ade  <- extract(sif, ade)
cf_ade   <- unlist(cf_ade, use.names = FALSE)[-c(1)]
fcom_ade <- unlist(fcom_ade, use.names = FALSE)[-c(1)]
fsat_ade <- unlist(fsat_ade, use.names = FALSE)[-c(1)]
sif_ade  <- unlist(sif_ade, use.names = FALSE)[-c(1)]

cf_das   <- extract(cf, das)
fcom_das <- extract(fcom, das)
fsat_das <- extract(fsat, das)
sif_das  <- extract(sif, das)
cf_das   <- unlist(cf_das, use.names = FALSE)[-c(1)]
fcom_das <- unlist(fcom_das, use.names = FALSE)[-c(1)]
fsat_das <- unlist(fsat_das, use.names = FALSE)[-c(1)]
sif_das  <- unlist(sif_das, use.names = FALSE)[-c(1)]

cf_dry   <- extract(cf, dry)
fcom_dry <- extract(fcom, dry)
fsat_dry <- extract(fsat, dry)
sif_dry  <- extract(sif, dry)
cf_dry   <- unlist(cf_dry, use.names = FALSE)[-c(1)]
fcom_dry <- unlist(fcom_dry, use.names = FALSE)[-c(1)]
fsat_dry <- unlist(fsat_dry, use.names = FALSE)[-c(1)]
sif_dry  <- unlist(sif_dry, use.names = FALSE)[-c(1)]

cf_lom   <- extract(cf, lom)
fcom_lom <- extract(fcom, lom)
fsat_lom <- extract(fsat, lom)
sif_lom  <- extract(sif, lom)
cf_lom   <- unlist(cf_lom, use.names = FALSE)[-c(1)]
fcom_lom <- unlist(fcom_lom, use.names = FALSE)[-c(1)]
fsat_lom <- unlist(fsat_lom, use.names = FALSE)[-c(1)]
sif_lom  <- unlist(sif_lom, use.names = FALSE)[-c(1)]

cf_fyo   <- extract(cf, fyo)
fcom_fyo <- extract(fcom, fyo)
fsat_fyo <- extract(fsat, fyo)
sif_fyo  <- extract(sif, fyo)
cf_fyo   <- unlist(cf_fyo, use.names = FALSE)[-c(1)]
fcom_fyo <- unlist(fcom_fyo, use.names = FALSE)[-c(1)]
fsat_fyo <- unlist(fsat_fyo, use.names = FALSE)[-c(1)]
sif_fyo  <- unlist(sif_fyo, use.names = FALSE)[-c(1)]

cf_ar1   <- extract(cf, ar1)
fcom_ar1 <- extract(fcom, ar1)
fsat_ar1 <- extract(fsat, ar1)
sif_ar1  <- extract(sif, ar1)
cf_ar1   <- unlist(cf_ar1, use.names = FALSE)[-c(1)]
fcom_ar1 <- unlist(fcom_ar1, use.names = FALSE)[-c(1)]
fsat_ar1 <- unlist(fsat_ar1, use.names = FALSE)[-c(1)]
sif_ar1  <- unlist(sif_ar1, use.names = FALSE)[-c(1)]

cf_ar2   <- extract(cf, ar2)
fcom_ar2 <- extract(fcom, ar2)
fsat_ar2 <- extract(fsat, ar2)
sif_ar2  <- extract(sif, ar2)
cf_ar2   <- unlist(cf_ar2, use.names = FALSE)[-c(1)]
fcom_ar2 <- unlist(fcom_ar2, use.names = FALSE)[-c(1)]
fsat_ar2 <- unlist(fsat_ar2, use.names = FALSE)[-c(1)]
sif_ar2  <- unlist(sif_ar2, use.names = FALSE)[-c(1)]

cf_mon   <- extract(cf, mon)
fcom_mon <- extract(fcom, mon)
fsat_mon <- extract(fsat, mon)
sif_mon  <- extract(sif, mon)
cf_mon   <- unlist(cf_mon, use.names = FALSE)[-c(1)]
fcom_mon <- unlist(fcom_mon, use.names = FALSE)[-c(1)]
fsat_mon <- unlist(fsat_mon, use.names = FALSE)[-c(1)]
sif_mon  <- unlist(sif_mon, use.names = FALSE)[-c(1)]

# Monthly aggregation
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

# Aggregate to month
cf_k34_month   <- to_month(cf_k34, 2019)
fcom_k34_month <- to_month(fcom_k34, 2019)
fsat_k34_month <- to_month(fsat_k34, 2019)
sif_k34_month  <- to_month(sif_k34, 2019)

cf_ade_month   <- to_month(cf_ade, 2019)
fcom_ade_month <- to_month(fcom_ade, 2019)
fsat_ade_month <- to_month(fsat_ade, 2019)
sif_ade_month  <- to_month(sif_ade, 2019)

cf_das_month   <- to_month(cf_das, 2019)
fcom_das_month <- to_month(fcom_das, 2019)
fsat_das_month <- to_month(fsat_das, 2019)
sif_das_month  <- to_month(sif_das, 2019)

cf_dry_month   <- to_month(cf_dry, 2019)
fcom_dry_month <- to_month(fcom_dry, 2019)
fsat_dry_month <- to_month(fsat_dry, 2019)
sif_dry_month  <- to_month(sif_dry, 2019)

cf_lom_month   <- to_month(cf_lom, 2019)
fcom_lom_month <- to_month(fcom_lom, 2019)
fsat_lom_month <- to_month(fsat_lom, 2019)
sif_lom_month  <- to_month(sif_lom, 2019)

cf_fyo_month   <- to_month(cf_fyo, 2019)
fcom_fyo_month <- to_month(fcom_fyo, 2019)
fsat_fyo_month <- to_month(fsat_fyo, 2019)
sif_fyo_month  <- to_month(sif_fyo, 2019)

cf_ar1_month   <- to_month(cf_ar1, 2019)
fcom_ar1_month <- to_month(fcom_ar1, 2019)
fsat_ar1_month <- to_month(fsat_ar1, 2019)
sif_ar1_month  <- to_month(sif_ar1, 2019)

cf_ar2_month   <- to_month(cf_ar2, 2019)
fcom_ar2_month <- to_month(fcom_ar2, 2019)
fsat_ar2_month <- to_month(fsat_ar2, 2019)
sif_ar2_month  <- to_month(sif_ar2, 2019)

cf_mon_month   <- to_month(cf_mon, 2019)
fcom_mon_month <- to_month(fcom_mon, 2019)
fsat_mon_month <- to_month(fsat_mon, 2019)
sif_mon_month  <- to_month(sif_mon, 2019)

# DFs for plotting
df_tower <- data.frame(ade_gpp_mean, das_gpp_mean, dry_gpp_mean, k34_gpp, lom_gpp_mean,
                       fyo_gpp_mean, ar1_gpp_mean, ar2_gpp_mean, mon_gpp_mean)
df_cf   <- data.frame(cf_ade_month, cf_das_month, cf_dry_month, cf_k34_month, cf_lom_month,
                      cf_fyo_month, cf_ar1_month, cf_ar2_month, cf_mon_month)
df_fcom <- data.frame(fcom_ade_month, fcom_das_month, fcom_dry_month, fcom_k34_month, fcom_lom_month,
                      fcom_fyo_month, fcom_ar1_month, fcom_ar2_month, fcom_mon_month)
df_fsat <- data.frame(fsat_ade_month, fsat_das_month, fsat_dry_month, fsat_k34_month, fsat_lom_month,
                      fsat_fyo_month, fsat_ar1_month, fsat_ar2_month, fsat_mon_month)
df_sif  <- data.frame(sif_ade_month, sif_das_month, sif_dry_month, sif_k34_month, sif_lom_month,
                      sif_fyo_month, sif_ar1_month, sif_ar2_month, sif_mon_month)
# Tower names
tower_names <- towers[,1]
tower_names <- c(tower_names[1:3], "BR-K34", tower_names[4:8])

## Plot variables
x            <- 1:12
x_lab        <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
y_lab_gpp    <- bquote("GPP (g C m"^"-2"*" d"^"-1"*")")
y_lab_sif    <- bquote("SIF"[Daily]*" (mW/m"^"2"*"/sr/nm)")
y_limit_gpp  <- c(0,16)
y_limit_sif  <- c(0,1)

### PLOT
cairo_pdf(out_name, width = 7.5, height = 6.25)

par(mfrow = c(3, 3), oma=c(1,3.5,2.5,3.5))

for (i in 1:9){
  
  op <- par(mar = c(2.5,0,0,0))
  
  plot(x, df_tower[,i], col = "black", type = "l", axes = FALSE, lwd = 2, lty = 2, xaxs="i", ylim = y_limit_gpp, ylab = "")
  
  lines(x, df_fcom[,i], col = "#FE6100", lwd = 2)
  lines(x, df_fsat[,i], col = "#648FFF", lwd = 2)
  lines(x, df_cf[,i], col = "#DC267F", lwd = 2)
  
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
  plot(x, df_sif[,i], col = "gray50", type = "l", axes = FALSE, lwd = 2, lty = 2, xaxs="i", ylim = y_limit_sif, ylab = "")
  
  if (i == 3 || i == 6 || i == 9) {
    axis(4, tck = 0.03, labels = TRUE, mgp=c(3, 0.2, 0), las = 2)
  } else {
    axis(4, tck = 0.03, labels = FALSE, mgp=c(3, 0.2, 0), las = 2)
  }

  
  mtext(3, text = tower_names[i], line = 0)
  
  if (i == 1) {
      legend("topleft", legend=c("Tower", "ChloFluo", "FluxCom", "FluxSat", "SIF"),
         col=c("black", "#DC267F", "#FE6100", "#648FFF", "gray50"), lty = c(2,1,1,1,2),
         lwd = 1, ncol = 2, y.intersp = 1, cex = 0.75, bty = "n", bg =)
  }
  
  # Regressions
  cf_reg   <- summary(lm(df_cf[,i] ~ df_tower[,i]))
  fcom_reg <- summary(lm(df_fcom[,i] ~ df_tower[,i]))
  fsat_reg <- summary(lm(df_fsat[,i] ~ df_tower[,i]))
  sif_reg  <- summary(lm(df_sif[,i] ~ df_tower[,i]))

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

mtext(1, text = "Month", outer = TRUE, line = -0.5)
mtext(2, text = as.expression(y_lab_gpp), outer = TRUE, line = 1.5)
mtext(4, text = as.expression(y_lab_sif), outer = TRUE, line = 2.5)

dev.off()