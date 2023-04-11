library(terra)
library(viridis)
library(rgdal)
library(RColorBrewer)

out_file <- "G:/ChloFluo/comps/c3_c4_map.pdf"

map        <- rast("G:/ChloFluo/input/C3C4/ISLSCP/c4_percent_1d.nc")
coastlines <- vect("C:/Russell/R_Scripts/TROPOMI_2/mapping/GSHHS_shp/c/GSHHS_c_L1.shp")

map_masked <- mask(map, coastlines)

col.map  <- viridis(10)

pdf(out_file, width=7.5, height=3.5, compress=FALSE)

plot(map_masked, col = col.map, ext = c(-180,180,-60,90))

mtext(3, text = "C4 Ratio")

dev.off()
