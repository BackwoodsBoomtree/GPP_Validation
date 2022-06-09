library(raster)
library(viridis)
library(rgdal)
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

# normalize <- function(x) {
#   u <- cellStats(x, stat="mean", na.rm=T)
#   o <- cellStats(x, stat="sd", na.rm=T)
#   return((x - u) / (o))
# }

# Min-Max normalization function
min_max_norm <- function(x) {
  (x - minValue(x)) / (maxValue(x) - minValue(x))
}

#### Load Map ####

coastlines <- readOGR("C:/Russell/R_Scripts/TROPOMI_2/mapping/GSHHS_shp/c/GSHHS_c_L1.shp")
class(coastlines)
extent(coastlines)
crs(coastlines)

#### Load the data ####

cf_annual_mean   <- raster("G:/ChloFluo/product/v01/1deg/clipfill/annual/ChloFluo.GPP.v01.1deg.CF80.2019.clipfill.annual.nc")
sif_annual_mean  <- raster("G:/ChloFluo/input/SIF/1deg/SIFqc.8day.1deg.CF80.2019.annual.nc")
r2_map           <- raster("G:/ChloFluo/comps/sif/raster_regressions/ChlFluo_vs_TROPOMI_SIF.v01.1deg.CF80.2019.clipfill_Rsquare.tif")
pval_map         <- raster("G:/ChloFluo/comps/sif/raster_regressions/ChlFluo_vs_TROPOMI_SIF.v01.1deg.CF80.2019.clipfill_Pval.tif")

# Mask sif by cf
sif_annual_mean <- mask(sif_annual_mean, cf_annual_mean)

# normalized difference
cf_annual_mean   <- setMinMax(cf_annual_mean)
sif_annual_mean  <- setMinMax(sif_annual_mean)
cf_norm          <- min_max_norm(cf_annual_mean)
sif_norm         <- min_max_norm(sif_annual_mean)
diff_annual_mean <- cf_norm - sif_norm

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
mean_sif <- round2(cellStats(sif_annual_mean, mean, na.rm=T), 2)

# Set thresholds for visual mapping
cf_annual_mean[cf_annual_mean > 8 ]   <- 8
sif_annual_mean[sif_annual_mean > 0.5 ] <- 0.5


#### PLOT STUFF ####

blue_ramp <- colorRampPalette(colors = c("#053061", "#F7F7F7"))
diff.col  <- c(blue_ramp(8), "#F4A582")
gpp.col   <- viridis(100)
r2.col    <- plasma(100)

labs <- c(expression(paste("ChloFluo Mean Annual GPP 2019")),
          expression(paste("TROPOMI Mean Annual SIF 2019")),
          expression(paste("Normalized Difference: ChloFluo - SIF")),
          expression(paste(" ChloFluo GPP vs TROPOMI SIF")))

pdf("G:/ChloFluo/comps/sif/ChloFluo_vs_TROPOSIF_comparisons_black.pdf", width=7.5, height=6, compress=FALSE)

par(mfrow=c(3,2),oma=c(0,0.25,1.25,0), bg = "black")

##### ChloFluo ####

op <- par(mar = c(0,0,0.25,0.25), bg = "black")
plot(cf_annual_mean, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, col = NA)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA, border = "white")
plot(coastlines, add = TRUE, border = NA, col = rgb(0.20,0.20,0.20))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA, border = "white")
plot(cf_annual_mean, col=gpp.col,  ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, add=T)
mtext(3, text=labs[1], cex= 0.85, col = "white")
mtext(3, text="a", cex= 0.85, adj=0, font=2, col = "white")

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


##### SIF ####

op <- par(mar = c(0,0,0.25,0.25), bg = "black")
plot(sif_annual_mean, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, col = NA)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA, border = "white")
plot(coastlines, add = TRUE, border = NA, col = rgb(0.20,0.20,0.20))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA, border = "white")
plot(sif_annual_mean, col=gpp.col,  ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, add=T)
mtext(3, text=labs[2], cex= 0.85, col = "white")
mtext(3, text="b", cex= 0.85, adj=0, font=2, col = "white")

plot(sif_annual_mean, legend.only=TRUE, col=gpp.col, horizontal=T, legend.width=2, legend.shrink=0.75,
     legend.args = list(text=expression(paste("mW m"^"2"*" sr nm")), side = 1, line = -2, cex=0.70, col = "white"),
     axis.args = list(line = -1.05, cex.axis=1, tick=F, at=c(minValue(sif_annual_mean),0.5), labels=c("0",">0.5"), col.axis = "white"),
     smallplot=c(0.40,0.90,0.2,0.25)); par(mar = par("mar"))

par(new=TRUE)
op <- par(mar = c(2.1,0.25,8,21)) # Set margins
hist(sif_annual_mean, col=rgb(0.30,0.30,0.30), breaks=10, ylim=c(0,4500), xlim=c(0,0.5), xaxs="i", yaxs="i", ann=FALSE, axes=FALSE)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "black", border = "white")
abline(v=mean_sif, col=rgb(206,97,75,max=255))
axis(3, tck=F, labels=mean_sif, at=mean_sif, mgp=c(3, 0.1, 0), col.axis = "white")
hist(sif_annual_mean, col=rgb(0.30,0.30,0.30), breaks=10, ylim=c(0,4500), xlim=c(0,0.5), xaxs="i", yaxs="i", ann=FALSE, axes=FALSE, add=T)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = NA, border = "white")


##### Difference in Annual Means ####
op <- par(mar = c(0,0,0.25,0.25))

med <- round(median(diff_annual_mean, na.rm=T), 2)
diff_annual_mean[diff_annual_mean < -0.4] <- -0.4
diff_annual_mean[diff_annual_mean > 0]  <- 0
plot(diff_annual_mean, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, col = NA)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA, border = "white")
plot(coastlines, add = TRUE, border = NA, col = rgb(0.20,0.20,0.20))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA, border = "white")
plot(diff_annual_mean, col=diff.col, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", legend=F, horizontal=T, add=T)
mtext(3, text=labs[3], cex=0.85, col = "white")
mtext(3, text="c", cex= 0.85, adj=0, font=2, col = "white")


plot(diff_annual_mean, legend.only=TRUE, col=diff.col, horizontal=T, legend.width=2, legend.shrink=0.75,
     legend.args = list(text=expression(paste("Difference (Normalized)")), side = 1, line = -2, cex=0.85, col = "white"),
     axis.args = list(line = -1.05, cex.axis=1,tick=F, at=c(-0.4,-0.2,0), labels=c("<-0.40","-0.2",">0"), col.axis = "white"),
     smallplot=c(0.40,0.90,0.2,0.25)); par(mar = par("mar"))

par(new=TRUE)
op <- par(mar = c(3,2,8,21)) # Set margins
barplot(diff_lat, col = NA, axes=F, tck=F, xpd=F, mgp=c(3,0.3,0), ann=FALSE, xaxs = "i", yaxs = "i", horiz = TRUE, border = NA, space = 0, xlim = c(-0.4, 0.4))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black", border = "white")
abline(h = 300)
rect(-100,35,100,105, col = rgb(0.20,0.20,0.20), border = NA)
rect(-100,47,100,93, col = rgb(0.30,0.30,0.30), border = NA)
par(new=T)
barplot(diff_lat, col = ifelse(diff_lat < 0, rgb(79,147,194,max=255), rgb(206,97,75,max=255)),
        axes=F, tck=F, xpd=F, mgp=c(3,0.3,0), ann=FALSE, xaxs = "i", yaxs = "i", horiz = TRUE, border = NA, space = 0, xlim = c(-0.4, 0.4))
abline(v = 0, col = "white")

# axis
axis(side=2, tck=F, las=1, cex.axis=1, labels = c("-60°", "0°", "80°"), mgp=c(3,0.3,0), at=c(1, 70,140), col.axis = "white")
axis(side=1, tck=F, las=1, cex.axis=1, labels = c("-0.4", "0", "0.4"),
     mgp=c(3,0.3,0), at=c(-0.4, 0, 0.4), col.axis = "white")
title <- as.list(expression(paste("Mean " , Delta), "by Latitude"))
mtext(do.call(expression, title), side = 3, line = c(1,-0.25), cex = 0.85, col = "white")
box(col = "white")


##### R2 map ####
op <- par(mar = c(0,0,0.25,0.25))
plot(r2_map, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, col = NA)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA, border = "white")
plot(coastlines, add = TRUE, border = NA, col = rgb(0.20,0.20,0.20))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA, border = "white")
plot(r2_map, zlim=c(0.1,1), col=r2.col, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", legend=F, horizontal=T, add=T)
mtext(3, text=labs[4], cex=0.85, col = "white")
mtext(3, text="d", cex= 0.85, adj=0, font=2, col = "white")


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
