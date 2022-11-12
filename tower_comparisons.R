library(terra)
library(viridis)
library(rgdal)
library(RColorBrewer)

out_name <- "G:/ChloFluo/comps/tower_comparisons.pdf"

cf_file    <- "G:/ChloFluo/product/v01/1deg/clipfill/ChloFluo.GPP.v01.1deg.CF80.2019.clipfill.nc"
fcom_file  <- "G:/FluxCom/RS/GPP.RS_V006.FP-ALL.MLM-ALL.METEO-NONE.720_360.8daily.2019.nc"
fsat_files <- list.files("G:/FluxSat", full.names = TRUE, pattern = "*.nc")
k34_gep    <- read.csv("G:/SIF_comps/figs/Wu_2016/K34_GEP.csv", header = FALSE)[,2]
k67_gep    <- read.csv("G:/SIF_comps/figs/Wu_2016/K67_GEP.csv", header = FALSE)[,2]
rja_gep    <- read.csv("G:/SIF_comps/figs/Wu_2016/RJA_GEP.csv", header = FALSE)[,2]
cax_gep    <- read.csv("G:/SIF_comps/figs/Wu_2016/CAX_GEP.csv", header = FALSE)[,2]
tch_gpp    <- read.csv("G:/ChloFluo/comps/tower-data/CG-Tch/FLX_CG-Tch_FLUXNET2015_FULLSET_MM_2006-2009_1-4.csv", header = TRUE)
ank_gpp    <- read.csv("G:/ChloFluo/comps/tower-data/GH-Ank/FLX_GH-Ank_FLUXNET2015_FULLSET_MM_2011-2014_1-4.csv", header = TRUE)

# Get Africa GPP tower data
tch_gpp_nt <- tch_gpp$GPP_NT_VUT_MEAN
ank_gpp_nt <- ank_gpp$GPP_NT_VUT_MEAN

tch_gpp_dt <- tch_gpp$GPP_DT_VUT_MEAN
ank_gpp_dt <- ank_gpp$GPP_DT_VUT_MEAN

tch_gpp_mean <- (tch_gpp_nt + tch_gpp_dt) / 2
ank_gpp_mean <- (ank_gpp_nt + ank_gpp_dt) / 2

# Get model GPP data
cf   <- rast(cf_file, subds = "gpp")
fcom <- rast(fcom_file, subds = "GPP")
fsat <- rast(fsat_files, subds = "GPP")

### Aggregate FluxSat

# Do 8-day temporal aggregation (non-leap year)
for (i in seq(1, 365, by = 8)) {
  
  if (i != 361) {
    fsat_mean_tmp <- mean(fsat[[i : (i + 7)]], na.rm = TRUE)
  } else {
    fsat_mean_tmp <- mean(fsat[[i : (i +  4)]], na.rm = TRUE)
  }
  
  if (i == 1) {
    fsat_8day <- fsat_mean_tmp
  } else {
    fsat_8day <- c(fsat_8day, fsat_mean_tmp)
  }
}

fsat <- aggregate(fsat_8day, 20, fun = mean, na.rm = TRUE)

## Tower
k34 <- cbind(-60.2093, -2.6091)
k34 <- vect(k34, crs = "+proj=longlat +ellps=WGS84")

cax <- cbind(-51.4536, -1.7483)
cax <- vect(cax, crs = "+proj=longlat +ellps=WGS84")

k67 <- cbind(-54.959, -2.857)
k67 <- vect(k67, crs = "+proj=longlat +ellps=WGS84")

tch <- cbind(11.6564, -4.2892)
tch <- vect(tch, crs = "+proj=longlat +ellps=WGS84")

ank <- cbind(-2.6942, 5.2685)
ank <- vect(ank, crs = "+proj=longlat +ellps=WGS84")

# Timeseries
cf_k34   <- extract(cf, k34)
fcom_k34 <- extract(fcom, k34)
fsat_k34 <- extract(fsat, k34)

cf_cax   <- extract(cf, cax)
fcom_cax <- extract(fcom, cax)
fsat_cax <- extract(fsat, cax)

cf_k67   <- extract(cf, k67)
fcom_k67 <- extract(fcom, k67)
fsat_k67 <- extract(fsat, k67)

cf_tch   <- extract(cf, tch)
fcom_tch <- extract(fcom, tch)
fsat_tch <- extract(fsat, tch)

cf_ank   <- extract(cf, ank)
fcom_ank <- extract(fcom, ank)
fsat_ank <- extract(fsat, ank)

cf_k34   <- unlist(cf_k34, use.names = FALSE)[-c(1)]
fcom_k34 <- unlist(fcom_k34, use.names = FALSE)[-c(1)]
fsat_k34 <- unlist(fsat_k34, use.names = FALSE)[-c(1)]

cf_cax   <- unlist(cf_cax, use.names = FALSE)[-c(1)]
fcom_cax <- unlist(fcom_cax, use.names = FALSE)[-c(1)]
fsat_cax <- unlist(fsat_cax, use.names = FALSE)[-c(1)]

cf_k67   <- unlist(cf_k67, use.names = FALSE)[-c(1)]
fcom_k67 <- unlist(fcom_k67, use.names = FALSE)[-c(1)]
fsat_k67 <- unlist(fsat_k67, use.names = FALSE)[-c(1)]

cf_tch   <- unlist(cf_tch, use.names = FALSE)[-c(1)]
fcom_tch <- unlist(fcom_tch, use.names = FALSE)[-c(1)]
fsat_tch <- unlist(fsat_tch, use.names = FALSE)[-c(1)]

cf_ank   <- unlist(cf_ank, use.names = FALSE)[-c(1)]
fcom_ank <- unlist(fcom_ank, use.names = FALSE)[-c(1)]
fsat_ank <- unlist(fsat_ank, use.names = FALSE)[-c(1)]

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
    s[i] <- weighted.mean(x, ww)
  }
  
  return(s)
  
}

# Aggregate to month
cf_k34_month   <- to_month(cf_k34, 2019)
fcom_k34_month <- to_month(fcom_k34, 2019)
fsat_k34_month <- to_month(fsat_k34, 2019)

cf_cax_month   <- to_month(cf_cax, 2019)
fcom_cax_month <- to_month(fcom_cax, 2019)
fsat_cax_month <- to_month(fsat_cax, 2019)

cf_k67_month   <- to_month(cf_k67, 2019)
fcom_k67_month <- to_month(fcom_k67, 2019)
fsat_k67_month <- to_month(fsat_k67, 2019)

cf_tch_month   <- to_month(cf_tch, 2019)
fcom_tch_month <- to_month(fcom_tch, 2019)
fsat_tch_month <- to_month(fsat_tch, 2019)

cf_ank_month   <- to_month(cf_ank, 2019)
fcom_ank_month <- to_month(fcom_ank, 2019)
fsat_ank_month <- to_month(fsat_ank, 2019)

tch_gpp_mean_month <- (tch_gpp_mean[1:12] + tch_gpp_mean[13:24] + tch_gpp_mean[25:36] + tch_gpp_mean[37:48]) / 4
ank_gpp_mean_month <- (ank_gpp_mean[1:12] + ank_gpp_mean[13:24] + ank_gpp_mean[37:48]) / 3


## Plot variables
x            <- 1:12
x_lab        <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
y_lab_gpp    <- list(bquote("GPP (g C m"^"-2"*" d"^"-1"*")"))
y_limit_gpp  <- c(4,12)

mag.cols <- magma(7)
vir.cols <- viridis(7)


### PLOT
cairo_pdf(out_name, width = 7.5, height = 3.25)

par(mfrow = c(1, 1), oma=c(2,3,2.5,0))

# Amazon GPP
op <- par(mar = c(0.5,0,0,0.5))

plot(x, fsat_k34_month, col = vir.cols[6], type = "l", axes = FALSE, lwd = 2, xaxs="i", ylim = y_limit_gpp, ylab = "")

lines(x, fcom_k34_month, col = vir.cols[4], lwd = 2)
lines(x, k34_gep, col = "black", lwd = 2)
lines(x, cf_k34_month, col = vir.cols[2], lwd = 2)

axis(1, tck = 0.04, labels = x_lab, at = x, mgp=c(3, 0.2, 0))
axis(2, tck = 0.04, mgp=c(3, 0.2, 0), las = 2)


mtext(2, text = do.call(expression, y_lab_gpp), line = 1.5)
mtext(3, text = "K34 Tower Site in the Amazon Basin", line = 0.5)

legend("topleft", legend=c("K34 ", "ChloFluo", "FluxCom", "FluxSat"),
       col=c("black", vir.cols[2], vir.cols[4], vir.cols[6]), lty = c(1,1,1,1),
       lwd = c(2,2,2,2), horiz = TRUE, y.intersp = 1, cex = 0.75)
box()

mtext(1, text = "Month", outer = TRUE, line = 1)

dev.off()


cf_reg <- lm(cf_k34_month~k34_gep)
summary(cf_reg)

fsat_reg <- lm(fsat_k34_month~k34_gep)
summary(fsat_reg)

fcom_reg <- lm(fcom_k34_month~k34_gep)
summary(fcom_reg)
