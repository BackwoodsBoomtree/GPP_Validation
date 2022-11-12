library(terra)
library(raster)
library(viridis)
library(RColorBrewer)

#### Output PDF name ####
out_name    <- "G:/ChloFluo/comps/maps/tropical_regions_map.pdf"


coastlines     <- vect("C:/Russell/R_Scripts/TROPOMI_2/mapping/GSHHS_shp/c/GSHHS_c_L1.shp")

tropics_ext    <- ext(c(-180, 180, -23.5, 23.5))

mcd12_majority <- rast("G:/MCD12C1/2020/reprocessed/percent/MCD12C1.A2020001.006.Percent_LC_03.tif")
mcd12_majority <- crop(mcd12_majority, tropics_ext)
mcd12_majority[mcd12_majority < 50]  <- NA

roi_amazon <- vect("G:/Amazon_shp/Amazon_poly.shp") # Amazon

congo_lon <- c(10, 29, 29, 10, 10)
congo_lat <- c(3, 3, -3, -3, 3)
roi_congo <- cbind(congo_lon, congo_lat)
roi_congo <- vect(roi_congo, crs = "+proj=longlat +ellps=WGS84", type = "polygons")

seasia_lon <- c(95, 155, 155, 95, 95)
seasia_lat <- c(10, 10, -11, -11, 10)
roi_seasia <- cbind(seasia_lon, seasia_lat)
roi_seasia <- vect(roi_seasia, crs = "+proj=longlat +ellps=WGS84", type = "polygons")

map.cols  <- colorRampPalette(c("#e5f5e0", "#006d2c"))
map.cols  <- (map.cols(11))
map.col = "#006d2c"


#### Plot ####
cairo_pdf(out_name, width = 6.5, height = 3)

par(oma=c(0.5,0.5,1.25,0.5))

# Map
plot(coastlines, border = NA, col = rgb(0.30,0.30,0.30), mar = c(0,0,0,0), axes = FALSE)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black")
plot(coastlines, border = NA, col = rgb(0.30,0.30,0.30), mar = c(0,0,0,0), axes = FALSE, add = TRUE)
plot(mcd12_majority, col = map.col, add = TRUE, legend = FALSE, axes = FALSE )
plot(roi_amazon, border = "white", col = NA, add = TRUE, axes = FALSE)
plot(roi_congo, border = "white", col = NA, add = TRUE, axes = FALSE)
plot(roi_seasia, border = "white", col = NA, add = TRUE, axes = FALSE)

mtext(3, text = "Tropical Forest")

dev.off()