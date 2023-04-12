library(terra)
library(viridis)
library(rgdal)
library(RColorBrewer)

out_name <- "G:/ChloFluo/comps/tropical_comparisons_v2.pdf"

cf_file    <- "G:/ChloFluo/product/v02/clipfill/ChloFluo.GPP.v02.1deg.CF80.2019.clipfill.nc"
fcom_file  <- "G:/FluxCom/RS/GPP.RS_V006.FP-ALL.MLM-ALL.METEO-NONE.720_360.8daily.2019.nc"
fsat_files <- list.files("G:/FluxSat/original", full.names = TRUE, pattern = "*.nc")
sif_file   <- "G:/ChloFluo/input/SIF/1deg/SIFqc.8day.1deg.CF80.2019.nc"
sify_file  <- "G:/ChloFluo/input/yield/1deg/yield.2019.8-day.1deg.nc"
apar_file  <- "G:/ChloFluo/input/APARchl/1deg/apar.2019.8-day.1deg.nc"

roi_amazon <- vect("G:/Amazon_shp/Amazon_poly.shp") # Amazon

congo_lon <- c(10, 29, 29, 10, 10)
congo_lat <- c(3, 3, -3, -3, 3)
roi_congo <- cbind(congo_lon, congo_lat)
roi_congo <- vect(roi_congo, crs = "+proj=longlat +ellps=WGS84", type = "polygons")

seasia_lon <- c(95, 155, 155, 95, 95)
seasia_lat <- c(10, 10, -11, -11, 10)
roi_seasia <- cbind(seasia_lon, seasia_lat)
roi_seasia <- vect(roi_seasia, crs = "+proj=longlat +ellps=WGS84", type = "polygons")

cf   <- rast(cf_file, subds = "gpp")
fcom <- rast(fcom_file, subds = "GPP")
fsat <- rast(fsat_files, subds = "GPP")
sif  <- rast(sif_file, subds = "sif743_qc")
sify <- rast(sify_file, subds = "yield")
apar <- rast(apar_file, subds = "apar")


###

### Aggregate FluxSat

# Do temporal aggregation (non-leap year)
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
  gc()
}

fsat <- aggregate(fsat_8day, 20, fun = mean, na.rm = TRUE)

# Clip to regions
cf_amazon   <- mask(cf, roi_amazon)
fcom_amazon <- mask(fcom, roi_amazon)
fsat_amazon <- mask(fsat, roi_amazon)
sif_amazon  <- mask(sif, roi_amazon)
sify_amazon <- mask(sify, roi_amazon)
apar_amazon <- mask(apar, roi_amazon)
cf_amazon_ts   <- global(cf_amazon, fun = "mean", na.rm = TRUE)[,1]
fcom_amazon_ts <- global(fcom_amazon, fun = "mean", na.rm = TRUE)[,1]
fsat_amazon_ts <- global(fsat_amazon, fun = "mean", na.rm = TRUE)[,1]
sif_amazon_ts  <- global(sif_amazon, fun = "mean", na.rm = TRUE)[,1]
sify_amazon_ts <- global(sify_amazon, fun = "mean", na.rm = TRUE)[,1]*100
apar_amazon_ts <- global(apar_amazon, fun = "mean", na.rm = TRUE)[,1]

cf_congo   <- mask(cf, roi_congo)
fcom_congo <- mask(fcom, roi_congo)
fsat_congo <- mask(fsat, roi_congo)
sif_congo  <- mask(sif, roi_congo)
sify_congo <- mask(sify, roi_congo)
apar_congo <- mask(apar, roi_congo)
cf_congo_ts   <- global(cf_congo, fun = "mean", na.rm = TRUE)[,1]
fcom_congo_ts <- global(fcom_congo, fun = "mean", na.rm = TRUE)[,1]
fsat_congo_ts <- global(fsat_congo, fun = "mean", na.rm = TRUE)[,1]
sif_congo_ts  <- global(sif_congo, fun = "mean", na.rm = TRUE)[,1]
sify_congo_ts <- global(sify_congo, fun = "mean", na.rm = TRUE)[,1]*100
apar_congo_ts <- global(apar_congo, fun = "mean", na.rm = TRUE)[,1]

cf_seasia   <- mask(cf, roi_seasia)
fcom_seasia <- mask(fcom, roi_seasia)
fsat_seasia <- mask(fsat, roi_seasia)
sif_seasia  <- mask(sif, roi_seasia)
sify_seasia <- mask(sify, roi_seasia)
apar_seasia <- mask(apar, roi_seasia)
cf_seasia_ts   <- global(cf_seasia, fun = "mean", na.rm = TRUE)[,1]
fcom_seasia_ts <- global(fcom_seasia, fun = "mean", na.rm = TRUE)[,1]
fsat_seasia_ts <- global(fsat_seasia, fun = "mean", na.rm = TRUE)[,1]
sif_seasia_ts  <- global(sif_seasia, fun = "mean", na.rm = TRUE)[,1]
sify_seasia_ts <- global(sify_seasia, fun = "mean", na.rm = TRUE)[,1]*100
apar_seasia_ts <- global(apar_seasia, fun = "mean", na.rm = TRUE)[,1]

## Plot variables
x            <- 1:46
xpos         <- seq(1, 46, by = 4)
xlabs        <- seq(1, 365, by = 32)
y_lab_gpp    <- list(bquote("GPP"), bquote("(g C m"^"-2"*" d"^"-1"*")"))
y_lab_sif    <- list(bquote("SIFdaily"), bquote("(mW/m"^"2"*"/sr/nm)"))
y_lab_sify   <- list(bquote("ΦF"), bquote("(%)"))
y_lab_apar   <- list(bquote("APARchl"), bquote("(W/m"^"2"*")"))
y_lab_nirv   <- list(bquote("NIRv"), bquote("(Reflectance)"))
y_limit_gpp  <- c(4,11)
y_limit_sif  <- c(0.2,0.6)
y_limit_sify <- c(0.35,0.55)
y_limit_apar <- c(50,140)

mag.cols <- magma(7)
vir.cols <- viridis(7)


### PLOT
cairo_pdf(out_name, width = 7.5, height = 5.25)

par(mfrow = c(4, 3), oma=c(2.0,6,2.5,0.1))

# Amazon GPP
op <- par(mar = c(0.5,0,0,0.5))

plot(x, cf_amazon_ts, col = vir.cols[2], type = "l", axes = FALSE, lwd = 1.5, xaxs="i", ylim = y_limit_gpp, ylab = "")

lines(x, fcom_amazon_ts, col = vir.cols[4], lwd = 1.5)
lines(x, fsat_amazon_ts, col = vir.cols[6], lwd = 1.5)

axis(1, tck = 0.03, labels = FALSE, at = x)
axis(1, tck = 0.06, labels = FALSE, at = xpos)
axis(2, tck = 0.03, mgp=c(3, 0.2, 0), las = 2)

mtext(2, text = do.call(expression, y_lab_gpp), line = c(4.25, 2.25))
mtext(3, text = "Amazon Basin", line = 0.5)
legend("topleft", legend=c("ChloFluo", "FluxCom", "FluxSat"),
       col=c(vir.cols[2], vir.cols[4], vir.cols[6]), lty = c(1,1,1),
       horiz = TRUE, y.intersp = 1, cex = 0.75)
box()


# CONGO GPP
op <- par(mar = c(0.5,0,0,0.5))

plot(x, cf_congo_ts, col = vir.cols[2], type = "l", axes = FALSE, lwd = 1.5, xaxs="i", ylim = y_limit_gpp, ylab = "")

lines(x, fcom_congo_ts, col = vir.cols[4], lwd = 1.5)
lines(x, fsat_congo_ts, col = vir.cols[6], lwd = 1.5)

axis(1, tck = 0.03, labels = FALSE, at = x)
axis(1, tck = 0.06, labels = FALSE, at = xpos)
axis(2, tck = 0.03, labels = FALSE, las = 2)

mtext(3, text = "Congo Basin", line = 0.5)
box()


# SEASIA GPP
op <- par(mar = c(0.5,0,0,0.5))

plot(x, cf_seasia_ts, col = vir.cols[2], type = "l", axes = FALSE, lwd = 1.5, xaxs="i", ylim = y_limit_gpp, ylab = "")

lines(x, fcom_seasia_ts, col = vir.cols[4], lwd = 1.5)
lines(x, fsat_seasia_ts, col = vir.cols[6], lwd = 1.5)

axis(1, tck = 0.03, labels = FALSE, at = x)
axis(1, tck = 0.06, labels = FALSE, at = xpos)
axis(2, tck = 0.03, labels = FALSE, las = 2)

mtext(3, text = "SE Asia", line = 0.5)
box()

# Amazon apar
op <- par(mar = c(0.5,0,0,0.5))

plot(x, apar_amazon_ts, col = mag.cols[5], type = "l", axes = FALSE, lwd = 1.5, xaxs="i", ylim = y_limit_apar, ylab = "")

axis(1, tck = 0.03, labels = FALSE, at = x)
axis(1, tck = 0.06, labels = FALSE, at = xpos)
axis(2, tck = 0.03, mgp=c(3, 0.2, 0), las = 2)

mtext(2, text = do.call(expression, y_lab_apar), line = c(4.25, 2.25))
box()


# Congo apar
op <- par(mar = c(0.5,0,0,0.5))

plot(x, apar_congo_ts, col = mag.cols[5], type = "l", axes = FALSE, lwd = 1.5, xaxs="i", ylim = y_limit_apar, ylab = "")

axis(1, tck = 0.03, labels = FALSE, at = x)
axis(1, tck = 0.06, labels = FALSE, at = xpos)
axis(2, tck = 0.03, labels = FALSE, las = 2)

box()


# SEASIA apar
op <- par(mar = c(0.5,0,0,0.5))

plot(x, apar_seasia_ts, col = mag.cols[5], type = "l", axes = FALSE, lwd = 1.5, xaxs="i", ylim = y_limit_apar, ylab = "")

axis(1, tck = 0.03, labels = FALSE, at = x)
axis(1, tck = 0.06, labels = FALSE, at = xpos)
axis(2, tck = 0.03, labels = FALSE, las = 2)

box()


# Amazon SIF
op <- par(mar = c(0.5,0,0,0.5))

plot(x, sif_amazon_ts, col = mag.cols[4], type = "l", axes = FALSE, lwd = 1.5, xaxs="i", ylim = y_limit_sif, ylab = "")

axis(1, tck = 0.03, labels = FALSE, at = x)
axis(1, tck = 0.06, labels = FALSE, at = xpos)
axis(2, tck = 0.03, mgp=c(3, 0.2, 0), las = 2)

mtext(2, text = do.call(expression, y_lab_sif), line = c(4.25, 2.25))
box()


# Congo SIF
op <- par(mar = c(0.5,0,0,0.5))

plot(x, sif_congo_ts, col = mag.cols[4], type = "l", axes = FALSE, lwd = 1.5, xaxs="i", ylim = y_limit_sif, ylab = "")

axis(1, tck = 0.03, labels = FALSE, at = x)
axis(1, tck = 0.06, labels = FALSE, at = xpos)
axis(2, tck = 0.03, labels = FALSE, las = 2)

box()


# SEASIA SIF
op <- par(mar = c(0.5,0,0,0.5))

plot(x, sif_seasia_ts, col = mag.cols[4], type = "l", axes = FALSE, lwd = 1.5, xaxs="i", ylim = y_limit_sif, ylab = "")

axis(1, tck = 0.03, labels = FALSE, at = x)
axis(1, tck = 0.06, labels = FALSE, at = xpos)
axis(2, tck = 0.03, labels = FALSE, las = 2)

box()

# Amazon sify
op <- par(mar = c(0.5,0,0,0.5))

plot(x, sify_amazon_ts, col = mag.cols[6], type = "l", axes = FALSE, lwd = 1.5, xaxs="i", ylim = y_limit_sify, ylab = "")

axis(1, tck = 0.03, labels = FALSE, at = x)
axis(1, tck = 0.06, labels = xlabs, at = xpos, mgp=c(3, 0.2, 0))
axis(2, tck = 0.03, mgp=c(3, 0.2, 0), las = 2)

mtext(2, text = do.call(expression, y_lab_sify), line = c(4.25, 2.25))
box()


# Congo sify
op <- par(mar = c(0.5,0,0,0.5))

plot(x, sify_congo_ts, col = mag.cols[6], type = "l", axes = FALSE, lwd = 1.5, xaxs="i", ylim = y_limit_sify, ylab = "")

axis(1, tck = 0.03, labels = FALSE, at = x)
axis(1, tck = 0.06, labels = xlabs, at = xpos, mgp=c(3, 0.2, 0))
axis(2, tck = 0.03, labels = FALSE, las = 2)

box()


# SEASIA sify
op <- par(mar = c(0.5,0,0,0.5))

plot(x, sify_seasia_ts, col = mag.cols[6], type = "l", axes = FALSE, lwd = 1.5, xaxs="i", ylim = y_limit_sify, ylab = "")

axis(1, tck = 0.03, labels = FALSE, at = x)
axis(1, tck = 0.06, labels = xlabs, at = xpos, mgp=c(3, 0.2, 0))
axis(2, tck = 0.03, labels = FALSE, las = 2)

box()

mtext(1, text = "Day of Year", outer = TRUE, line = 1)

dev.off()

