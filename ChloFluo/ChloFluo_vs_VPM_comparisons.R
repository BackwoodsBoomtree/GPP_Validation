library(raster)
library(RColorBrewer)

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
class(coastlines)
extent(coastlines)
crs(coastlines)

#### Load the data ####

cf_annual_mean   <- raster("G:/ChloFluo/product/v01/1deg/clipfill/annual/ChloFluo.GPP.v01.1deg.CF80.2019.clipfill.annual.nc")
vpm_annual_mean  <- raster("G:/ChloFluo/comps/vpm/VPM.1deg.2019.annual.nc")
diff_annual_mean <- raster("G:/ChloFluo/comps/vpm/ChloFluo-VPM.v01.1deg.CF80.2019.annual.nc")
r2_map           <- raster("G:/ChloFluo/comps/vpm/raster_regressions/ChlFluo_vs_VPM.v01.1deg.CF80.2019.clipfill_Rsquare.tif")
pval_map         <- raster("G:/ChloFluo/comps/vpm/raster_regressions/ChlFluo_vs_VPM.v01.1deg.CF80.2019.clipfill_Pval.tif")

# Mask vpm by cf
vpm_annual_mean <- mask(vpm_annual_mean, cf_annual_mean)

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
mean_vpm <- round2(cellStats(vpm_annual_mean, mean, na.rm=T), 2)

# Set thresholds for visual mapping
diff_annual_mean[diff_annual_mean < -3.0] <- -3.0
diff_annual_mean[diff_annual_mean > 3.0]  <- 3.0
cf_annual_mean[cf_annual_mean > 8 ]   <- 8
vpm_annual_mean[vpm_annual_mean > 8 ] <- 8


#### PLOT STUFF ####

diff.col <- rev(brewer.pal(n = 11, name = "RdBu"))
gpp.col <- viridis(100)
r2.col <- plasma(100)

labs <- c(expression(paste("ChloFluo Mean Annual GPP 2019")),
          expression(paste("VPM Mean Annual GPP 2019")),
          expression(paste("GPP Difference: ChloFluo - VPM")),
          expression(paste("Correlation for ChloFluo vs VPM")))

pdf("G:/ChloFluo/comps/vpm/ChloFluo_vs_VPM_comparisons.pdf", width=7.5, height=6, compress=FALSE)

par(mfrow=c(3,2),oma=c(0,0.25,1.25,0))

##### ChloFluo ####

op <- par(mar = c(0,0,0.25,0.25))
plot(cf_annual_mean, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, col = NA)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA)
plot(coastlines, add = TRUE, border = NA, col = "gray")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA)
plot(cf_annual_mean, col=gpp.col,  ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, add=T)
mtext(3, text=labs[1], cex= 0.85)
mtext(3, text="a", cex= 0.85, adj=0, font=2)
plot(coastlines, add = TRUE, lwd=0.5)

plot(cf_annual_mean, legend.only=TRUE, col=gpp.col, horizontal=T, legend.width=2, legend.shrink=0.75,
     legend.args = list(text=expression(paste("g C m"^"-2"*"day"^"-1")), side = 1, line = -2, cex=0.70),
     axis.args = list(line = -1.05, cex.axis=1,tick=F, at=c(0,8), labels=c("0",">8")),
     smallplot=c(0.40,0.90,0.2,0.25)); par(mar = par("mar"))

par(new=TRUE)
op <- par(mar = c(2.1,0.25,8,21)) # Set margins
hist(cf_annual_mean, col="gray75", breaks=10, ylim=c(0,4000), xlim=c(0,8), xaxs="i", yaxs="i", ann=FALSE, axes=FALSE)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white",border = NA)
abline(v=mean_cf, col="red")
axis(3, tck=F, labels=mean_cf, at=mean_cf, mgp=c(3, 0.1, 0))
hist(cf_annual_mean, col="gray75", breaks=10, ylim=c(0,4000), xlim=c(0,8), xaxs="i", yaxs="i", ann=FALSE, axes=FALSE, add=T)
box()


##### VPM ####

op <- par(mar = c(0,0,0.25,0.25))
plot(vpm_annual_mean, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, col = NA)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA)
plot(coastlines, add = TRUE, border = NA, col = "gray")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA)
plot(vpm_annual_mean, col=gpp.col,  ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, add=T)
mtext(3, text=labs[2], cex= 0.85)
mtext(3, text="b", cex= 0.85, adj=0, font=2)
plot(coastlines, add = TRUE, lwd=0.5)

plot(vpm_annual_mean, legend.only=TRUE, col=gpp.col, horizontal=T, legend.width=2, legend.shrink=0.75,
     legend.args = list(text=expression(paste("g C m"^"-2"*"day"^"-1")), side = 1, line = -2, cex=0.70),
     axis.args = list(line = -1.05, cex.axis=1,tick=F, at=c(0,8), labels=c("0",">8")),
     smallplot=c(0.40,0.90,0.2,0.25)); par(mar = par("mar"))

par(new=TRUE)
op <- par(mar = c(2.1,0.25,8,21)) # Set margins
hist(vpm_annual_mean, col="gray75", breaks=10, ylim=c(0,4000), xlim=c(0,8), xaxs="i", yaxs="i", ann=FALSE, axes=FALSE)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white",border = NA)
abline(v=mean_vpm, col="red")
axis(3, tck=F, labels=mean_vpm, at=mean_vpm, mgp=c(3, 0.1, 0))
hist(vpm_annual_mean, col="gray75", breaks=10, ylim=c(0,4000), xlim=c(0,8), xaxs="i", yaxs="i", ann=FALSE, axes=FALSE, add=T)
box()


##### Difference in Annual Means ####
op <- par(mar = c(0,0,0.25,0.25))

med <- round(median(diff_annual_mean, na.rm=T), 2)
diff_annual_mean[diff_annual_mean < -3.0] <- -3.0
diff_annual_mean[diff_annual_mean > 3.0]  <- 3.0
plot(diff_annual_mean, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, col = NA)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA)
plot(coastlines, add = TRUE, border = NA, col = "gray")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA)
plot(diff_annual_mean, col=diff.col, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", legend=F, horizontal=T, add=T)
mtext(3, text=labs[3], cex=0.85)
mtext(3, text="c", cex= 0.85, adj=0, font=2)
plot(coastlines, add = TRUE, lwd=0.5)

plot(diff_annual_mean, legend.only=TRUE, col=diff.col, horizontal=T, legend.width=2, legend.shrink=0.75,
     legend.args = list(text=expression(paste("g C m"^"-2"*"day"^"-1")), side = 1, line = -2, cex=0.85),
     axis.args = list(line = -1.05, cex.axis=1,tick=F, at=c(-3.0,0,3.0), labels=c("<-3.0","0",">3.0")),
     smallplot=c(0.40,0.90,0.2,0.25)); par(mar = par("mar"))

par(new=TRUE)
op <- par(mar = c(3,2,8,21)) # Set margins
barplot(diff_lat, col = NA, axes=F, tck=F, xpd=F, mgp=c(3,0.3,0), ann=FALSE, xaxs = "i", yaxs = "i", horiz = TRUE, border = NA, space = 0, xlim = c(-3.0, 3.0))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white",border = NA)
abline(h = 300)
rect(-100,35,100,105, col = rgb(0.5,0.5,0.5,1/4), border = NA)
rect(-100,47,100,93, col = rgb(0.35,0.35,0.35,1/4), border = NA)
par(new=T)
barplot(diff_lat, col = ifelse(diff_lat < 0, rgb(79,147,194,max=255), rgb(206,97,75,max=255)),
        axes=F, tck=F, xpd=F, mgp=c(3,0.3,0), ann=FALSE, xaxs = "i", yaxs = "i", horiz = TRUE, border = NA, space = 0, xlim = c(-3.0, 3.0))
abline(v = 0)
box()
# axis
axis(side=2, tck=F, las=1, cex.axis=1, labels = c("-60°", "0°", "80°"), mgp=c(3,0.3,0), at=c(1, 70,140))
axis(side=1, tck=F, las=1, cex.axis=1, labels = c("-3.0", "0", "3.0"),
     mgp=c(3,0.3,0), at=c(-3.0, 0, 3.0))
title <- as.list(expression(paste("Mean " , Delta), "by Latitude"))
mtext(do.call(expression, title), side = 3, line = c(1,-0.25), cex = 0.85)


##### R2 map ####
op <- par(mar = c(0,0,0.25,0.25))
plot(r2_map, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, col = NA)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA)
plot(coastlines, add = TRUE, border = NA, col = "gray")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA)
plot(r2_map, zlim=c(0.1,1), col=r2.col, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", legend=F, horizontal=T, add=T)
mtext(3, text=labs[4], cex=0.85)
mtext(3, text="d", cex= 0.85, adj=0, font=2)
plot(coastlines, add = TRUE, lwd=0.5)

plot(r2_map, zlim=c(0.1,1), legend.only=TRUE, col=r2.col, horizontal=T, legend.width=2, legend.shrink=0.75,
     legend.args = list(text=expression(paste("R"^"2")), side = 1, line = -2, cex=0.85),
     axis.args = list(line = -1.05, cex.axis=1,tick=F, at=c(0.1,1), labels=c("0.1","1")),
     smallplot=c(0.40,0.90,0.2,0.25)); par(mar = par("mar"))

par(new=TRUE)
op <- par(mar = c(2.1,0.25,8,21))  
hist(r2_map, col="gray75", breaks=9, ylim=c(0,4000), xlim=c(0.1,1), xaxs="i", yaxs="i", ann=FALSE, axes=FALSE)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white",border = NA)
med <- round(median(as.vector(r2_map), na.rm = T), 2)
abline(v=med, col="red")
axis(3, tck=F, labels=med, at=med, mgp=c(3, 0.1, 0))
hist(r2_map, col="gray75", breaks=9, ylim=c(0,4000), xlim=c(0.1,1), xaxs="i", yaxs="i", ann=FALSE, axes=FALSE, add=T)
box()


dev.off()
