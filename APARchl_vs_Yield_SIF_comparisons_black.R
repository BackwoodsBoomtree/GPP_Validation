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

yield_beta_map <- raster("G:/ChloFluo/comps/aparchl_vs_yield_sif/raster_regressions/APARchl_vs_SIFyield+SIF.v01.1deg.CF80.2019.clipfill_Beta_Standard_sif_yield.tif")
yield_pval_map <- raster("G:/ChloFluo/comps/aparchl_vs_yield_sif/raster_regressions/APARchl_vs_SIFyield+SIF.v01.1deg.CF80.2019.clipfill_Pval_sif_yield.tif")

sif_beta_map   <- raster("G:/ChloFluo/comps/aparchl_vs_yield_sif/raster_regressions/APARchl_vs_SIFyield+SIF.v01.1deg.CF80.2019.clipfill_Beta_Standard_sif743_qc.tif")
sif_pval_map   <- raster("G:/ChloFluo/comps/aparchl_vs_yield_sif/raster_regressions/APARchl_vs_SIFyield+SIF.v01.1deg.CF80.2019.clipfill_Pval_sif743_qc.tif")

dif_map        <- abs(sif_beta_map) - abs(yield_beta_map)

# Mask variables out by pvalue
yield_pval_map[yield_pval_map > 0.05] <- NA
sif_pval_map[sif_pval_map > 0.05]     <- NA

yield_beta_map <- mask(yield_beta_map, yield_pval_map)
sif_beta_map   <- mask(sif_beta_map, sif_pval_map)

# Thresholds for visual mapping
yield_beta_map[yield_beta_map < -1] <- -1
sif_beta_map[sif_beta_map > 1]      <- 1
dif_map[dif_map > 0.3]              <- 0.3
dif_map[dif_map < -0.3]             <- -0.3

#### PLOT STUFF ####

stand.col     <- viridis(12)
stand.col.neg <- stand.col[1:6]
stand.col.pos <- stand.col[6:12]
diff.col      <- brewer.pal(n = 7, name = "RdBu")

labs <- c(expression(paste("Standardized Beta for ΦF")),
          expression(paste("Standardized Beta for SIF")),
          expression(paste("Standard Beta SIF - Standard Beta ΦF")))

pdf("G:/ChloFluo/comps/aparchl_vs_yield_sif/APARchl_vs_Yield_SIF_comparisons_black.pdf", width=7.5, height=6, compress=FALSE)

par(mfrow=c(3,2),oma=c(0,0.25,1.25,0), bg = "black")

##### SIFyield ####

op <- par(mar = c(0,0,0.25,0.25), bg = "black")
plot(yield_beta_map, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, col = NA)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA, border = "white")
plot(coastlines, add = TRUE, border = NA, col = rgb(0.20,0.20,0.20))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA, border = "white")
plot(yield_beta_map, col=stand.col.neg,  ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, add=T)
mtext(3, text=labs[1], cex= 0.85, col = "white")
mtext(3, text="a", cex= 0.85, adj=0, font=2, col = "white")

plot(yield_beta_map, legend.only=TRUE, col=stand.col.neg, horizontal=T, legend.width=2, legend.shrink=0.75,
     legend.args = list(text=expression(paste("Standardized Beta Coefficient")), side = 1, line = -2, cex=0.70, col = "white"),
     axis.args = list(line = -1.05, cex.axis=1, tick=F, at=c(-1,-0.02), labels=c("< -1","0"), col.axis = "white"),
     smallplot=c(0.40,0.90,0.2,0.25)); par(mar = par("mar"))


##### SIF ####

op <- par(mar = c(0,0,0.25,0.25), bg = "black")
plot(sif_beta_map, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, col = NA)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA, border = "white")
plot(coastlines, add = TRUE, border = NA, col = rgb(0.20,0.20,0.20))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA, border = "white")
plot(sif_beta_map, col=stand.col.pos,  ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, add=T)
mtext(3, text=labs[2], cex= 0.85, col = "white")
mtext(3, text="b", cex= 0.85, adj=0, font=2, col = "white")

plot(sif_beta_map, legend.only=TRUE, col=stand.col.pos, horizontal=T, legend.width=2, legend.shrink=0.75,
     legend.args = list(text=expression(paste("Standardized Beta Coefficient")), side = 1, line = -2, cex=0.70, col = "white"),
     axis.args = list(line = -1.05, cex.axis=1, tick=F, at=c(0.5,1), labels=c("0.5","> 1"), col.axis = "white"),
     smallplot=c(0.40,0.90,0.2,0.25)); par(mar = par("mar"))


##### Difference ####

op <- par(mar = c(0,0,0.25,0.25), bg = "black")
plot(dif_map, ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, col = NA)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA, border = "white")
plot(coastlines, add = TRUE, border = NA, col = rgb(0.20,0.20,0.20))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA, border = "white")
plot(dif_map, col=diff.col,  ext=c(-180,180,-80,80), axes=F, xaxs="i", yaxs="i", horizontal=T, legend=F, add=T)
mtext(3, text=labs[3], cex= 0.85, col = "white")
mtext(3, text="c", cex= 0.85, adj=0, font=2, col = "white")

plot(dif_map, legend.only=TRUE, col=diff.col, horizontal=T, legend.width=2, legend.shrink=0.75,
     legend.args = list(text=expression(paste("Standardized Beta Coefficient")), side = 1, line = -2, cex=0.70, col = "white"),
     axis.args = list(line = -1.05, cex.axis=1, tick=F, at=c(-0.3, 0 , 0.3), labels=c("< -0.3", "0", "> 0.3"), col.axis = "white"),
     smallplot=c(0.40,0.90,0.2,0.25)); par(mar = par("mar"))


dev.off()
