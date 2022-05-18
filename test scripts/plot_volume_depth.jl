include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise 
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools,
    PyPlot, JLD2, DrWatson, Statistics, JLD2
import CairoMakie as Mkie
import GeoMakie
using .OHC_helper
using GoogleDrive,NCDatasets, NetCDF

include(srcdir("config_exp.jl"))
pygui(false)
runpath,diagpath = listexperiments(exprootdir())
# abbreviations for each experiment for labels, etc.
shortnames = expnames()
marks = expsymbols()
nexp = length(shortnames) # number of experiments

theta_datadir = joinpath(datadir(), "θ_data") 
OHC_datadir = joinpath(datadir(), "OHC_data") 

# basin_name="Pacific"
basin_name="Atlantic"

exp1 = "const_density"
exp1_dir = joinpath(theta_datadir, "OHC", exp1)
basin_vol_pth = searchdir(OHC_datadir, basin_name*"Level")[1]
@load joinpath(OHC_datadir, basin_vol_pth) level_volumes
basin_volumes = level_volumes

tecco = collect(Float64, 1992+1/24:1/12:2018)

fig, ax = plt.subplots(1, 1, figsize=(10,10))

xlbl  = L"m^3ß"
ax.set_xlabel(xlbl)
ax.set_ylabel("Depth [m]")
ax.grid(true)

fcycle = 1 # units: yr^{-1}
overtones= 2; # semi-annual, quad-annual, etc.: helps capture an asymmetric seasonal cycle
Ecycle,Fcycle = seasonal_matrices(fcycle,tecco,overtones)

# returns "full" second
#full means that seasonal pattern has not been removed 

ax.errorbar(basin_volumes, z,xerr = 0, color = "black", marker = "o")
fig.subplots_adjust(top=0.9, left=0.1, right=0.9, bottom=0.12)
ax.legend(loc="lower center",bbox_to_anchor=(0.2, -0.2), ncol=3)
fig.suptitle(basin_name * " Total Volume", size = 10)
fig.tight_layout()
outputfile = plotsdir("OHC/volume_depth_plot_"*basin_name * ".pdf")
fig.savefig(outputfile)
clear("all")