#  Remove the interannual frequency energy in surface forcing fields.
#  Keep the interannual energy in a given region.
#  Diagnostic plots have been removed from this version.
#
#  Script argument: region (must be defined in `src/ECCOonPoseidon.jl`)
#  If no arguments are passed, then interannual variability is removed everywhere.
using Pkg
Pkg.activate(".")
include("../src/intro.jl")

using Revise
using ECCOtour, ECCOonPoseidon
using Statistics, Distributions, FFTW, LinearAlgebra, StatsBase
using MeshArrays, MITgcmTools
import PyPlot as plt
include(srcdir("config_exp.jl"))

# This could be put into src code for scientific project.
inputdir = "/vast/eccodrive/files/Version4/Release4/other/flux-forced-seasonalcycle/forcing/"
outputdir = "/vast/eccodrive/files/Version4/Release4/other/flux-forced-seasonalcycle/iter0_forcing/"
# read lat, lon at center of grid cell
(ϕC,λC) = latlonC(γ)
# on the vector (Staggered) grid
(ϕG,λG) = latlonG(γ)


(!isdir(outputdir))&&(mkpath(outputdir))

midname = "_6hourlyavg_"
# varnames = ("atmPload","oceFWflx","oceQsw","oceSflux","oceSPflx","oceTAUE","oceTAUN","oceTAUX",
#             "oceTAUY","sIceLoad","sIceLoadPatmPload","sIceLoadPatmPload_nopabar","TFLUX")
varnames = ["oceTAUX","oceTAUY", "TFLUX", "oceQsw", "oceFWflx", "oceSflux"]

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
nt6hr = length(fluxsample_point)

t6hr_start = 1/8 # 3Z Jan 1 1992 #in days
Δt6hr = 6/24 # 1/4 of a day
t6hr = range(t6hr_start,step=Δt6hr,length=nt6hr)

daysperyear = 365.25 # nt6hr / (6 * len(ye?ars)
fcycle = 1/(daysperyear) # units: day^{-1}
w = Float32.(1 / length(t6hr))

i0 = Dict()
i0["oceTAUX"] = read_bin(datadir("oceTAUx_i0_mean.data"),Float32,γ)
i0["oceTAUY"] = read_bin(datadir("oceTAUy_i0_mean.data"),Float32,γ)
i0["TFLUX"] = read_bin(datadir("TFLUX_i0_mean.data"),Float32,γ)
i0["oceQsw"] = read_bin(datadir("oceQsw_i0_mean.data"),Float32,γ)
i0["oceFWflx"] = read_bin(datadir("oceFWflx_i0_mean.data"),Float32,γ)
i0["oceSflux"] = read_bin(datadir("oceSflux_i0_mean.data"),Float32,γ)

for vname ∈ varnames
    filein = inputdir*vname*midname
    fileout = outputdir*vname*midname
    println(filein)

    # check for NaN's in output
    true_mean = matrixmean(w,filein,years,γ)
    Δ = i0[vname] .- true_mean
    Δ = Float32.(Δ)
    nyr = length(years)

    nseries = []

    for tt = 1:nyr
        fnamein = filein*string(years[tt])
        println("reading file "*fnamein)
        field = read_bin(fnamein,Float32,γ);
        cons_offset!(field, Δ) #remove the iteration 129 mean 
        fnameout = fileout*string(years[tt])
        ECCOtour.write(fnameout,field)
        println("saved file "*fnameout)

    end
end

vname = varnames[1]
seasonal_saved,nseries = extract_timeseries(outputdir*vname*midname,years,γ,xv,yv,fv)

seasonal_saved = Float32.(seasonal_saved)
println(mean(Float32, seasonal_saved))
println(i0[vname][fv][xv, yv])
println(mean(Float32, seasonal_saved .- i0[vname][fv][xv, yv]))

seasonal_saved_input,nseries = extract_timeseries(inputdir*vname*midname,years,γ,xv,yv,fv)
seasonal_saved_input = Float32.(seasonal_saved_input)
println(mean(Float32, seasonal_saved_input))

vname = varnames[2]
seasonal_saved,nseries = extract_timeseries(outputdir*vname*midname,years,γ,xv,yv,fv)
seasonal_saved = Float32.(seasonal_saved)
println(mean(Float32, seasonal_saved))
println(i0[vname][fv][xv, yv])
println(mean(Float32, seasonal_saved .- i0[vname][fv][xv, yv]))


fig, ax = plt.subplots(1, figsize=(15, 5))
ax.plot(1:nt, seasonal_saved, linewidth = 4, 
label = "Seasonal Cycle of TFLUX [Δt=6 hours]");
ax.legend()
fig

