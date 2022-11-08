library(terra)
library(viridis)
library(rgdal)
library(RColorBrewer)

out_file <- "G:/ChloFluo/comps/fluxcom/tower_comparisons_black.pdf"

cf_file    <- "G:/ChloFluo/product/v01/1deg/clipfill/ChloFluo.GPP.v01.1deg.CF80.2019.clipfill.nc"
fcom_file  <- "G:/FluxCom/RS/GPP.RS_V006.FP-ALL.MLM-ALL.METEO-NONE.720_360.8daily.2019.nc"
fsat_files <- list.files("G:/FluxSat", full.names = TRUE, pattern = "*.nc")

cf   <- rast(cf_file, subds = "gpp")
fcom <- rast(fcom_file, subds = "GPP")
fsat <- rast(fsat_files, subds = "GPP")

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
}

fsat <- aggregate(fsat_8day, 20, fun = mean, na.rm = TRUE)

## Tower
k34 <- cbind(-60.2093,-2.6091)
k34 <- vect(k34, crs = "+proj=longlat +ellps=WGS84")

# Timeseries
cf_k34   <- extract(cf, k34)
fcom_k34 <- extract(fcom, k34)
fsat_k34 <- extract(fsat, k34)

cf_k34   <- unlist(cf_k34, use.names = FALSE)[-c(1)]
fcom_k34 <- unlist(fcom_k34, use.names = FALSE)[-c(1)]
fsat_k34 <- unlist(fsat_k34, use.names = FALSE)[-c(1)]


