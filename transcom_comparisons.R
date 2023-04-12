library(terra)
library(viridis)
library(rgdal)
library(RColorBrewer)

save_name <- "G:/ChloFluo/comps/transcom_comparisons_v2.pdf"

cf_file    <- "G:/ChloFluo/product/v02/clipfill/ChloFluo.GPP.v02.1deg.CF80.2019.clipfill.nc"
fcom_file  <- "G:/FluxCom/RS/GPP.RS_V006.FP-ALL.MLM-ALL.METEO-NONE.720_360.8daily.2019.nc"
fsat_files <- "G:/FluxSat/monthly/1deg/GPP_FluxSat_8day_1deg_v2_2019.nc"
regions    <- rast("G:/TransCom/regions.nc", subds = "transcom_regions")
names_csv  <- read.csv("G:/TransCom/transcom_names.csv", header = FALSE)

cf   <- rast(cf_file, subds = "gpp")
fcom <- rast(fcom_file, subds = "GPP")
fsat <- rast(fsat_files, subds = "GPP")

# Collapse and place in right order names
names_csv <- names_csv[-(25)]
names_csv <- names_csv[-23,]
for (i in 1:nrow(names_csv)) {
  for (j in 1:ncol(names_csv)) {
    letter <- names_csv[i, j]
    if (j == 1) {
      out_name <- letter
    } else {
      out_name <- c(out_name, letter)
    }
  }
  out_name <- paste(out_name, sep="", collapse="")
  out_name <- trimws(out_name)
  if (i == 1) {
    names <- data.frame(out_name)
  } else {
    names <- rbind(names, out_name)
  }
}
names <- names[-c(12:22),]

# Set oceans to NA
regions[regions > 11] <- NA

# Reproject
ext(regions) <- ext(cf)
crs(regions) <- crs(cf)

# Fcom to 1 deg
fcom_one <- aggregate(fcom, 2, fun = mean, na.rm = TRUE)


## Plot variables
x            <- 1:46
xpos         <- seq(1, 46, by = 4)
xlabs        <- seq(1, 365, by = 64)
y_lab_gpp    <- bquote("Gross Primary Production (g C m"^"-2"*" d"^"-1"*")")
y_limit_gpp  <- c(0,14)

mag.cols <- magma(11)
vir.cols <- viridis(7)

map.cols <- colorRampPalette(c("#648FFF", "#785EF0", "#DC267F", "#FE6100","#FFB000"))
map.cols <- map.cols(11)
map.cols <- c(map.cols[1], map.cols[11], map.cols[2], map.cols[10], 
              map.cols[3], map.cols[9], map.cols[4], map.cols[8],
              map.cols[5], map.cols[7], map.cols[6])

### PLOT
cairo_pdf(save_name, width = 7.5, height = 7.25)

par(mfrow = c(4, 3), oma=c(3.0,3.5,0,0.1))

for (i in 1:length(names)) {
  op <- par(mar = c(0.5,0,2.5,0.5))
  
  # Mask
  cf_trans   <- mask(cf, regions, maskvalues = i, inverse = TRUE)
  fcom_trans <- mask(fcom_one, regions, maskvalues = i, inverse = TRUE)
  fsat_trans <- mask(fsat, regions, maskvalues = i, inverse = TRUE)
  
  # Time series
  cf_trans_ts   <- global(cf_trans, fun = "mean", na.rm = TRUE)[,1]
  fcom_trans_ts <- global(fcom_trans, fun = "mean", na.rm = TRUE)[,1]
  fsat_trans_ts <- global(fsat_trans, fun = "mean", na.rm = TRUE)[,1]
  
  plot(x, cf_trans_ts, col = vir.cols[2], type = "l", axes = FALSE, lwd = 1.5, xaxs="i", ylim = y_limit_gpp, ylab = "")
  
  lines(x, fcom_trans_ts, col = vir.cols[4], lwd = 1.5)
  lines(x, fsat_trans_ts, col = vir.cols[6], lwd = 1.5)
  
  axis(1, tck = 0.03, labels = FALSE, at = x)
  axis(1, tck = 0.06, labels = FALSE, at = xpos)
  
  if (i == 10 || i == 11){
    axis(1, tck = FALSE, mgp=c(3, 0.2, 0), labels = xlabs, at = seq(1, 46, by = 8))
  }
  
  if (i == 1 || i == 4 || i == 7 || i == 10){
    axis(2, tck = 0.03, mgp=c(3, 0.2, 0), las = 2)
  } else {
    axis(2, tck = 0.03, mgp=c(3, 0.2, 0), las = 2, labels = FALSE)
  }
  
  mtext(3, text = paste0(names[i], " (", i, ")"), line = 0.5)
  
  if (i == 1) {
    legend("topleft", legend=c("ChloFluo", "FluxCom", "FluxSat"),
           col=c(vir.cols[2], vir.cols[4], vir.cols[6]), lty = c(1,1,1),
           horiz = TRUE, y.intersp = 1, cex = 0.75)
  }
  
  box()
  
}

# Map
par(xpd=TRUE)
plot(regions, axes = FALSE, type = "classes", mar = NA, col = map.cols, 
     ext = c(-180, 180, -60, 80), legend = FALSE)
legend(-165, -50, legend = seq(1,11), horiz = FALSE, fill = map.cols, bty = "n",
       ncol = 6, x.intersp = 0.5)
mtext(3, text = "TransCom Regions", line = 0.5)
# box()

mtext(2, text = as.expression(y_lab_gpp), outer = TRUE, line = 1.5)
mtext(1, text = "Day of Year", outer = TRUE, line = 1.5)

dev.off()

