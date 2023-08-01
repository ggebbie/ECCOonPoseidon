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
inputdir = fluxdir()
outputdir = datadir("forcings/iter0_mean_tau/")
# read lat, lon at center of grid cell
(ϕC,λC) = latlonC(γ)
# on the vector (Staggered) grid
(ϕG,λG) = latlonG(γ)

(!isdir(outputdir))&&(mkpath(outputdir))

midname = "_6hourlyavg_"
# varnames = ("atmPload","oceFWflx","oceQsw","oceSflux","oceSPflx","oceTAUE","oceTAUN","oceTAUX",
#             "oceTAUY","sIceLoad","sIceLoadPatmPload","sIceLoadPatmPload_nopabar","TFLUX")
varnames = ("oceTAUX","oceTAUY")
Δi129i0 = Dict()

Δi129i0["oceTAUX"] = read_bin(datadir("oceTAUX_i129_i0_diff.data"),Float32,γ)
Δi129i0["oceTAUY"] = read_bin(datadir("oceTAUY_i129_i0_diff.data"),Float32,γ)

frootsample = inputdir*varnames[end]*midname
# sample calcs at one point. Get the length of timeseries in nseries.
yv = 10
xv = 70
fv = 4
tmplat  = ϕC[fv]; lat_point = tmplat[xv,yv]
tmplon  = λC[fv]; lon_point = tmplon[xv,yv]

years = 1992:2017

fluxsample_point,nseries = extract_timeseries(frootsample,years,γ,xv,yv,fv)
                        # take biweekly mean. Use triangular filter with break points at:
# take biweekly mean. Use triangular filter with break points at:
Δi14day = 4*14 # grid index range
nt6hr = length(fluxsample_point)
i14day = 2.5:Δi14day:nt6hr+56 # goes past end of time by 14 days to be sure

t6hr_start = 1/8 # 3Z Jan 1 1992 #in days
Δt6hr = 6/24 # 1/4 of a day
t6hr = range(t6hr_start,step=Δt6hr,length=nt6hr)

# Given a full-resolution timeseries, what values at the tiepoints best represent the timeseries?
# Determine tiepoints (in time) where fluxes are adjusted by ECCO in optimization procedure.
Δt14day = 14 # units: days
t14day_start = 0 # 12Z Jan 1 1992 by inspection of figure
t14day = range(t14day_start,step=Δt14day,stop=t6hr[end]+14)# goes past the end by one
nt14day = length(t14day)

# Values are added to the t14day tiepoints, then linearly interpolated to fill gaps.
# Careful not to store high-resolution data all at same time.
# Make a function that takes values at the tiepoints and then makes a full-resolution timeseries.
@time E14to6,F6to14 = get_filtermatrix(t6hr,t14day; ival = 1)

daysperyear = 365.25 # nt6hr / (6 * len(years)
fcycle = 1/(daysperyear) # units: day^{-1}
# for removing seasonal cycle from 14-day averaged timeseries
Ecycle,Fcycle = seasonal_matrices(fcycle,t14day,3)

for vname ∈ varnames
    filein = inputdir*vname*midname
    fileout = outputdir*vname*midname
    println(filein)

    # use this F to decompose controllable/uncontrollable parts
    # i.e., 6 hourly to 14 day
    flux_14day = matrixfilter(F6to14,filein,years,γ)

    # remove seasonal cycle from 14-day averaged timeseries
    # solve for seasonal cycle parameters
    # use 14-day because it's computationally efficient
    βcycle = matmul(Fcycle,flux_14day,γ)

    # reconstruct the full seasonal cycle.
    flux_14day_seasonal = matmul(Ecycle,βcycle,γ) #also solves for time mean
    cons_offset!(flux_14day_seasonal, -Δi129i0[vname]) #remove the iteration 129 mean 
    # check for NaN's in output
    (sum(nancount(flux_14day_seasonal)) > 0) && (error("NaNs in the filtered output"));

    matrixsaveinterpolation(E14to6,flux_14day_seasonal,filein,fileout,years,γ)
end

t = LinRange(years[1], years[end], nt6hr)

vname = varnames[1]
actual,nseries = extract_timeseries(inputdir*vname*midname,years,γ,xv,yv,fv)
seasonal_saved,nseries = extract_timeseries(outputdir*vname*midname,years,γ,xv,yv,fv)

fig, ax = plt.subplots(1, figsize=(15, 5))
ax.plot(t, actual, label = "Original TFLUX [Δt=6 hours]"); ax.set_xlabel("years")
ax.plot(t, seasonal_saved, linewidth = 4, 
label = "Seasonal Cycle of TFLUX [Δt=6 hours]");

println(mean(actual .- seasonal_saved))
ax.legend()
fig