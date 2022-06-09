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

#### Load Map ####

coastlines <- readOGR("C:/Russell/R_Scripts/TROPOMI_2/mapping/GSHHS_shp/c/GSHHS_c_L1.shp")
class(coastlines)
extent(coastlines)
crs(coastlines)

#### Load the data ####

cf_annual_mean   <- raster("G:/ChloFluo/product/v01/1deg/clipfill/annual/ChloFluo.GPP.v01.1deg.CF80.2019.clipfill.annual.nc")
r2_map           <- raster("G:/ChloFluo/comps/aparchl/raster_regressions/ChlFluo_vs_APARchl.v01.1deg.CF80.2019.clipfill_Rsquare.tif")
pval_map         <- raster("G:/ChloFluo/comps/aparchl/raster_regressions/ChlFluo_vs_APARchl.v01.1deg.CF80.2019.clipfill_Pval.tif")

# Mask variables out by pvalue
pval_map[pval_map > 0.05] <- NA
r2_map                    <- mask(r2_map, pval_map)

# Means
mean_cf  <- round2(cellStats(cf_annual_mean, mean, na.rm=T), 2)

# Set thresholds for visual mapping
cf_annual_mean[cf_annual_mean > 8 ]   <- 8


#### PLOT STUFF ####

blue_ramp <- colorRampPalette(colors = c("#053061", "#F7F7F7"))
gpp.col   <- viridis(100)
r2.col    <- plasma(100)

labs <- c(expression(paste("ChloFluo Mean Annual GPP 2019")),
          expression(paste("ChloFluo APARchl vs TROPOMI SIF")))

pdf("G:/ChloFluo/comps/aparchl/ChloFluo_vs_APARchl_comparisons_black.pdf", width=7.5, height=6, compress=FALSE)

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

##### R2 map ####
op <- par(mar = c(0,0,0.25,0.25))
plot(r2_map, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, col = NA)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA, border = "white")
plot(coastlines, add = TRUE, border = NA, col = rgb(0.20,0.20,0.20))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA, border = "white")
plot(r2_map, zlim=c(0.1,1), col=r2.col, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", legend=F, horizontal=T, add=T)
mtext(3, text=labs[2], cex=0.85, col = "white")
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
