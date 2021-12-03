library(terra)
library(viridis)

# To get the annual mean for a grid cell, we need to calculate 0 GPP for the time periods
# where there was an NA.

cf  <- "G:/ChloFluo/product/v01/1deg/clipfill/annual/ChloFluo.GPP.v01.1deg.CF80.2019.clipfill.annual.nc"
vpm <- "G:/ChloFluo/comps/vpm/VPM.1deg.2019.nc"

cf <- rast(cf)

vpm <- rast(vpm)
vpm <- mean(vpm, na.rm = TRUE)

annual_diff <- cf - vpm

writeCDF(annual_diff, "G:/ChloFluo/comps/vpm/ChloFluo-VPM.v01.1deg.CF80.2019.annual.nc",
         varname = "gpp_diff", longname = "ChloFluo GPP - VPM GPP", unit = "g C/m-2/day-1",
         missval = -9999, overwrite = TRUE)


writeCDF(vpm, "G:/ChloFluo/comps/vpm/VPM.1deg.2019.annual.nc",
         varname = "gpp", longname = "VPM GPP", unit = "g C/m-2/day-1",
         missval = -9999, overwrite = TRUE)

plot(annual_diff, col = viridis(15))
