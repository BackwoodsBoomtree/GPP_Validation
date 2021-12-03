
###############################################################################
#
# Make timeseries difference maps between VPM and ChloFluo
#
# 
#
###############################################################################

using NCDatasets
using GriddingMachine
using Dates
using Plots, Colors
include("../../ChloFlo/ChloFluo_GPP/src/save/save_nc.jl");

cf  = Dataset("/mnt/g/ChloFluo/product/v01/1deg/ChloFluo.GPP.v01.1deg.CF80.2019.nc")["gpp"][:,:,:];

zoom = 1; # spatial resolution is 1/zoom degree
vpm  = load_LUT(GPPVPMv20{Float32}(), 2019, "5X", "8D", 1);
vpm  = vpm.data;

gpp_diff = cf .- vpm;
gpp_diff = permutedims(gpp_diff, [2,1,3])
vpm      = permutedims(vpm, [2,1,3]) 

heatmap(gpp_diff[:,:,1], title = "VPM", bg = :white, color = :RdBu, clim = (-7.5,7.5))

save_nc(gpp_diff, "/mnt/g/ChloFluo/comps/vpm/ChloFluo-VPM.v01.1deg.CF80.2019.nc", 2019, "gpp_diff", "ChloFluo GPP - VPM GPP", "g C/m-2/day-1");

# Save VPM to calculate annual mean map in R
save_nc(vpm, "/mnt/g/ChloFluo/comps/vpm/VPM.1deg.2019.nc", 2019, "gpp", "VPM GPP", "g C/m-2/day-1");