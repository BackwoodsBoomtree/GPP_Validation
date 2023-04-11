library(raster)
library(viridis)
library(rgdal)
library(RColorBrewer)

out_file <- "G:/ChloFluo/comps/fluxsat/ChloFluo_vs_FluxSat_comparisons_v2.pdf"

# Round up
round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}

normalize <- function(x) {
  u <- cellStats(x, stat="mean", na.rm=T)
  o <- cellStats(x, stat="sd", na.rm=T)
  return((x - u) / (o))
}

#### Load Map ####

coastlines <- readOGR("C:/Russell/R_Scripts/TROPOMI_2/mapping/GSHHS_shp/c/GSHHS_c_L1.shp")

#### Load the data ####

cf_gpp           <- brick("G:/ChloFluo/product/v02/clipfill/ChloFluo.GPP.v02.1deg.CF80.2019.clipfill.nc")
gpp              <- brick("G:/FluxSat/monthly/1deg/GPP_FluxSat_8day_1deg_v2_2019.nc", varname = "GPP")
r2_map           <- raster("G:/ChloFluo/comps/fluxsat/raster_regressions/ChloFluo_vs_FluxSat.v02.1deg.CF80.2019.clipfill_Rsquare.tif")
pval_map         <- raster("G:/ChloFluo/comps/fluxsat/raster_regressions/ChloFluo_vs_FluxSat.v02.1deg.CF80.2019.clipfill_Pval.tif")

# Aggregate gpp
cf_annual_mean   <- mean(cf_gpp, na.rm = TRUE)
gpp_annual_mean  <- mean(gpp, na.rm = TRUE)

# Mask gpp by cf
gpp_annual_mean  <- mask(gpp_annual_mean, cf_annual_mean)

# Calculate difference map
diff_annual_mean <- cf_annual_mean - gpp_annual_mean

# Mask variables out by pvalue
pval_map[pval_map > 0.05] <- NA
r2_map                    <- mask(r2_map, pval_map)

# Row means for latitude mean difference
diff_lat  <- rev(rowMeans(as.matrix(diff_annual_mean), na.rm = TRUE))

# Trim row means to -60 and 80 latitude
diff_lat <- diff_lat[-(1:10)]
diff_lat <- diff_lat[-(61:80)]

# Means
mean_cf  <- round2(cellStats(cf_annual_mean, mean, na.rm=T), 2)
mean_gpp <- round2(cellStats(gpp_annual_mean, mean, na.rm=T), 2)

# Set thresholds for visual mapping
diff_annual_mean[diff_annual_mean < -3.0] <- -3.0
diff_annual_mean[diff_annual_mean > 3.0]  <- 3.0
cf_annual_mean[cf_annual_mean > 8 ]   <- 8
gpp_annual_mean[gpp_annual_mean > 8 ] <- 8


#### PLOT STUFF ####

diff.col <- rev(brewer.pal(n = 11, name = "RdBu"))
gpp.col  <- viridis(100)
r2.col   <- plasma(100)

labs <- c(expression(paste("ChloFluo Mean Annual GPP 2019")),
          expression(paste("FluxSat Mean Annual GPP 2019")),
          expression(paste("GPP Difference: ChloFluo - FluxSat")),
          expression(paste("Correlation for ChloFluo vs FluxSat")))

pdf(out_file, width=7.5, height=6, compress=FALSE)

par(mfrow=c(3,2),oma=c(0,0.25,1.25,0))

##### ChloFluo ####

op <- par(mar = c(0,0,0.25,0.25))
plot(cf_annual_mean, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, col = NA)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "black")
plot(coastlines, add = TRUE, border = NA, col = rgb(0.20,0.20,0.20))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4] ,col = NA, border = "black")
plot(cf_annual_mean, col=gpp.col,  ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, add=T)
mtext(3, text=labs[1], cex= 0.85)
mtext(3, text="a", cex= 0.85, adj=0, font=2)

plot(cf_annual_mean, legend.only=TRUE, col=gpp.col, horizontal=T, legend.width=2, legend.shrink=0.75,
     legend.args = list(text=expression(paste("g C m"^"-2"*"day"^"-1")), side = 1, line = -2, cex=0.70, col = "white"),
     axis.args = list(line = -1.05, cex.axis=1, tick=F, at=c(0,8), labels=c("0",">8"), col.axis = "white"),
     smallplot=c(0.40,0.90,0.2,0.25)); par(mar = par("mar"))

par(new=TRUE)
op <- par(mar = c(2.1,0.25,8,21)) # Set margins
hist(cf_annual_mean, col=rgb(0.30,0.30,0.30), breaks=10, ylim=c(0,4000), xlim=c(0,8), xaxs="i", yaxs="i", ann=FALSE, axes=FALSE)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "black", border = "white")
abline(v=mean_cf, col=rgb(206,97,75,max=255))
axis(3, tck=F, labels=mean_cf, at=mean_cf, mgp=c(3, 0.1, 0), col.axis = "white")
hist(cf_annual_mean, col=rgb(0.30,0.30,0.30), breaks=10, ylim=c(0,4000), xlim=c(0,8), xaxs="i", yaxs="i", ann=FALSE, axes=FALSE, add=T)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = NA, border = "white")


##### gpp ####

op <- par(mar = c(0,0,0.25,0.25))
plot(gpp_annual_mean, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, col = NA)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "black")
plot(coastlines, add = TRUE, border = NA, col = rgb(0.20,0.20,0.20))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = NA, border = "black")
plot(gpp_annual_mean, col=gpp.col,  ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, add=T)
mtext(3, text=labs[2], cex= 0.85)
mtext(3, text="b", cex= 0.85, adj=0, font=2)

plot(gpp_annual_mean, legend.only=TRUE, col=gpp.col, horizontal=T, legend.width=2, legend.shrink=0.75,
     legend.args = list(text=expression(paste("g C m"^"-2"*"day"^"-1")), side = 1, line = -2, cex=0.70, col = "white"),
     axis.args = list(line = -1.05, cex.axis=1, tick=F, at=c(minValue(gpp_annual_mean),8), labels=c("0",">8"), col.axis = "white"),
     smallplot=c(0.40,0.90,0.2,0.25)); par(mar = par("mar"))

par(new=TRUE)
op <- par(mar = c(2.1,0.25,8,21)) # Set margins
hist(gpp_annual_mean, col=rgb(0.30,0.30,0.30), breaks=10, ylim=c(0,4500), xlim=c(0,8), xaxs="i", yaxs="i", ann=FALSE, axes=FALSE)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "black", border = "white")
abline(v=mean_cf, col=rgb(206,97,75,max=255))
axis(3, tck=F, labels=mean_gpp, at=mean_gpp, mgp=c(3, 0.1, 0), col.axis = "white")
hist(gpp_annual_mean, col=rgb(0.30,0.30,0.30), breaks=10, ylim=c(0,4500), xlim=c(0,8), xaxs="i", yaxs="i", ann=FALSE, axes=FALSE, add=T)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = NA, border = "white")


##### Difference in Annual Means ####
op <- par(mar = c(0,0,0.25,0.25))

med <- round(median(diff_annual_mean, na.rm=T), 2)
diff_annual_mean[diff_annual_mean < -3.0] <- -3.0
diff_annual_mean[diff_annual_mean > 3.0]  <- 3.0
plot(diff_annual_mean, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, col = NA)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "black")
plot(coastlines, add = TRUE, border = NA, col = rgb(0.20,0.20,0.20))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA, border = "black")
plot(diff_annual_mean, col=diff.col, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", legend=F, horizontal=T, add=T)
mtext(3, text=labs[3], cex=0.85)
mtext(3, text="c", cex= 0.85, adj=0, font=2)


plot(diff_annual_mean, legend.only=TRUE, col=diff.col, horizontal=T, legend.width=2, legend.shrink=0.75,
     legend.args = list(text=expression(paste("g C m"^"-2"*"day"^"-1")), side = 1, line = -2, cex=0.85, col = "white"),
     axis.args = list(line = -1.05, cex.axis=1,tick=F, at=c(-3.0,0,3.0), labels=c("<-3.0","0",">3.0"), col.axis = "white"),
     smallplot=c(0.40,0.90,0.2,0.25)); par(mar = par("mar"))

par(new=TRUE)
op <- par(mar = c(3,2,8,21)) # Set margins
barplot(diff_lat, col = NA, axes=F, tck=F, xpd=F, mgp=c(3,0.3,0), ann=FALSE, xaxs = "i", yaxs = "i", horiz = TRUE, border = NA, space = 0, xlim = c(-3.0, 3.0))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black", border = "white")
abline(h = 300)
rect(-100,35,100,105, col = rgb(0.20,0.20,0.20), border = NA)
rect(-100,47,100,93, col = rgb(0.30,0.30,0.30), border = NA)
par(new=T)
barplot(diff_lat, col = ifelse(diff_lat < 0, rgb(79,147,194,max=255), rgb(206,97,75,max=255)),
        axes=F, tck=F, xpd=F, mgp=c(3,0.3,0), ann=FALSE, xaxs = "i", yaxs = "i", horiz = TRUE, border = NA, space = 0, xlim = c(-3.0, 3.0))
abline(v = 0, col = "white")

# axis
axis(side=2, tck=F, las=1, cex.axis=1, labels = c("-60°", "0°", "80°"), mgp=c(3,0.3,0), at=c(1, 70,140), col.axis = "white")
axis(side=1, tck=F, las=1, cex.axis=1, labels = c("-3.0", "0", "3.0"),
     mgp=c(3,0.3,0), at=c(-3.0, 0, 3.0), col.axis = "white")
title <- as.list(expression(paste("Mean " , Delta), "by Latitude"))
mtext(do.call(expression, title), side = 3, line = c(1,-0.25), cex = 0.85, col = "white")
box(col = "white")


##### R2 map ####
op <- par(mar = c(0,0,0.25,0.25))
plot(r2_map, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, col = NA)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "black")
plot(coastlines, add = TRUE, border = NA, col = rgb(0.20,0.20,0.20))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA, border = "black")
plot(r2_map, zlim=c(0.1,1), col=r2.col, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", legend=F, horizontal=T, add=T)
mtext(3, text=labs[4], cex=0.85)
mtext(3, text="d", cex= 0.85, adj=0, font=2)


plot(r2_map, zlim=c(0.1,1), legend.only=TRUE, col=r2.col, horizontal=T, legend.width=2, legend.shrink=0.75,
     legend.args = list(text=expression(paste("R"^"2")), side = 1, line = -2, cex=0.85, col = "white"),
     axis.args = list(line = -1.05, cex.axis=1,tick=F, at=c(0.1,1), labels=c("0.1","1"), col.axis = "white"),
     smallplot=c(0.40,0.90,0.2,0.25)); par(mar = par("mar"))

par(new=TRUE)
op <- par(mar = c(2.1,0.25,8,21)) # Set margins
hist(r2_map, col=rgb(0.30,0.30,0.30), breaks=9, ylim=c(0,4000), xlim=c(0.1,1), xaxs="i", yaxs="i", ann=FALSE, axes=FALSE)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "black", border = "white")
med <- round(median(as.vector(r2_map), na.rm = T), 2)
abline(v=med, col=rgb(206,97,75,max=255))
axis(3, tck=F, labels=med, at=med, mgp=c(3, 0.1, 0), col.axis = "white")
hist(r2_map, col=rgb(0.30,0.30,0.30), breaks=10, ylim=c(0,4000), xlim=c(0.1,1), xaxs="i", yaxs="i", ann=FALSE, axes=FALSE, add=T)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = NA, border = "white")


dev.off()
