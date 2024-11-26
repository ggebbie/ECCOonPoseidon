#  Remove the interannual frequency energy in surface forcing fields.
#  Keep the interannual energy in a given region.
#  Diagnostic plots have been removed from this version.
#
#  Script argument: region (must be defined in `src/ECCOonPoseidon.jl`)
#  If no arguments are passed, then interannual variability is removed everywhere.

include("../src/intro.jl")

using Revise
using ECCOtour, ECCOonPoseidon
using Statistics, Distributions, FFTW, LinearAlgebra, StatsBase
using MeshArrays, MITgcmTools
import PyPlot as plt
include(srcdir("config_exp.jl"))

# This could be put into src code for scientific project.
# read lat, lon at center of grid cell
(ϕC,λC) = latlonC(γ)
# on the vector (Staggered) grid
(ϕG,λG) = latlonG(γ)

yv = 10
xv = 70
fv = 4
tmplat  = ϕC[fv]; lat_point = tmplat[xv,yv]
tmplon  = λC[fv]; lon_point = tmplon[xv,yv]
years = 1992:2000

inputdir1 = fluxdir()
midname = "_6hourlyavg_"
# varnames = ("atmPload","oceFWflx","oceQsw","oceSflux","oceSPflx","oceTAUE","oceTAUN","oceTAUX",
#             "oceTAUY","sIceLoad","sIceLoadPatmPload","sIceLoadPatmPload_nopabar","TFLUX")
varnames = ("oceTAUX","oceTAUY")
frootsample = inputdir1*varnames[1]*midname
fluxsample_point,nseries = extract_timeseries(frootsample,years,γ,xv,yv,fv)

inputdir2 = "/vast/eccodrive/files/Version4/Release4/input_forcing/"
midname = ""
varnames = ("eccov4r4_ustr_","eccov4r4_vstr_")
frootsample = inputdir2*varnames[1]*midname
fluxsample_point2,nseries2 = extract_timeseries(frootsample,years,γ,xv,yv,fv)

inputdir3 = "/vast/eccodrive/files/Version4/Release4/other/input_forcing_unadjusted_test/"
midname = ""
varnames = ("eccov4r4_unadj_ustr_","eccov4r4_unadj_vstr_")
frootsample = inputdir3*varnames[1]*midname
fluxsample_point3,nseries3 = extract_timeseries(frootsample,years,γ,xv,yv,fv)

Δi14day = 4*14 # grid index range
nt6hr = length(fluxsample_point)
i14day = 2.5:Δi14day:nt6hr+56 # goes past end of time by 14 days to be sure

t6hr_start = 1/8 # 3Z Jan 1 1992 #in days
Δt6hr = 6/24 # 1/4 of a day
t6hr = collect(range(t6hr_start,step=Δt6hr,length=nt6hr))


fig, ax = plt.subplots(4, 1, sharey = false)
ax[1].plot(t6hr ./ 365, -fluxsample_point2, label = "Iteration 129 BF")
ax[2].plot(t6hr ./ 365, fluxsample_point, label = "Iteration 129 Fluxes")
ax[3].plot(t6hr./ 365, -fluxsample_point3, label = "Iteration 0 BF")
ax[4].plot(t6hr./ 365, -fluxsample_point2 .- (-fluxsample_point3), label = "'control adjustments'")

[a.set_ylim(-1, 1) for a in ax[1:3]]

[a.legend(frameon = false) for a in ax]
[a.grid(alpha = 0.4) for a in ax]

fig.tight_layout()
fig