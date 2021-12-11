library(terra)
library(viridis)

# Calculate difference between input data (ChloFluo minus other variable)
a_diff <- function(cf, variable, out_name, sname, lname, u) {
  
  cf       <- rast(cf)
  variable <- rast(variable)
  
  annual_diff <- cf - vpm
  
  writeCDF(annual_diff, out_name,
           varname = sname, longname = lname, unit = u,
           missval = -9999, overwrite = TRUE)
  
  return(annual_diff)
}

# Get annual mean of 8-day vpm data
vpm <- "G:/ChloFluo/comps/vpm/VPM.1deg.2019.nc"
vpm <- rast(vpm)
vpm <- mean(vpm, na.rm = TRUE)

writeCDF(vpm, "G:/ChloFluo/comps/vpm/VPM.1deg.2019.annual.nc",
         varname = "gpp", longname = "VPM GPP", unit = "g C/m-2/day-1",
         missval = -9999, overwrite = TRUE)

diff <- a_diff("G:/ChloFluo/product/v01/1deg/clipfill/annual/ChloFluo.GPP.v01.1deg.CF80.2019.clipfill.annual.nc",
               "G:/ChloFluo/comps/vpm/VPM.1deg.2019.annual.nc",
               "G:/ChloFluo/comps/vpm/ChloFluo-VPM.v01.1deg.CF80.2019.annual.nc",
               "gpp_diff",
               "ChloFluo GPP - VPM GPP",
               "g C/m-2/day-1")

plot(diff, col = viridis(15))