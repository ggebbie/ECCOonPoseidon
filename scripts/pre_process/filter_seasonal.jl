#  Remove the interannual frequency energy in surface forcing fields.
#  Keep the interannual energy in a given region.
#  Diagnostic plots have been removed from this version.
#
#  Script argument: region (must be defined in `src/ECCOonPoseidon.jl`)
#  If no arguments are passed, then interannual variability is removed everywhere.

include("../src/intro.jl")

using Revise
using ECCOtour, ECCOonPoseidon
using Statistics, PyPlot, Distributions, FFTW, LinearAlgebra, StatsBase
using MeshArrays, MITgcmTools

# from intro.jl, default is nointerannual
expt == "nointerannual" ? keepregion = false : keepregion = true
println("Experiment: ",expt)

include(srcdir("config_exp.jl"))

# This could be put into src code for scientific project.
inputdir = fluxdir()
outputdir = fluxdir(expt)

# read lat, lon at center of grid cell
(ϕC,λC) = latlonC(γ)
# on the vector (Staggered) grid
(ϕG,λG) = latlonG(γ)

if !isdir(outputdir)
    mkpath(outputdir)
end
midname = "_6hourlyavg_"
varnames = ["atmPload","oceFWflx","oceQsw","oceSflux","oceSPflx","oceTAUE","oceTAUN","oceTAUX",
            "oceTAUY","sIceLoad","sIceLoadPatmPload","sIceLoadPatmPload_nopabar","TFLUX"]
varnames = ["oceTAUX"]
# sample calcs at one point. Get the length of timeseries in nseries.
yv = 40
xv = 180
fv = 5
tmplat  = ϕC[fv]; lat_point = tmplat[xv,yv]
tmplon  = λC[fv]; lon_point = tmplon[xv,yv]
frootsample = inputdir*varnames[end]*midname
years = 1992:2017
fluxsample_point,nseries = extract_timeseries(frootsample,years,γ,xv,yv,fv)

nt6hr = length(fluxsample_point)

t6hr_start = 1/8 # 3Z Jan 1 1992
Δt6hr = 6/24 # 1/4 of a day
t6hr = range(t6hr_start,step=Δt6hr,length=nt6hr)

daysperyear = 365.25
fcycle = 1/(daysperyear) # units: day^{-1}
# for removing seasonal cycle from 14-day averaged timeseries
Ecycle,Fcycle = seasonal_matrices(fcycle,nt6hr)
Edag = Ecycle * Fcycle

#vname = varnames[1] # for interactive use
for vname ∈ varnames
    filein = inputdir*vname*midname
    fileout = outputdir*vname*midname
    println(filein)

    # reconstruct the full seasonal cycle.
    flux_6hr_seasonal = matrixfilter(Edag,filein,years,γ)

    # put tflux_14day_lopass on to 6hr
    # check for NaN's in output
    nancount_6hr_seasonal = sum(nancount(flux_6hr_seasonal))

    if nancount_6hr_seasonal > 0
        error("NaNs in the filtered output")
    end

    nyr = length(years)
    istart = 1
    nseries = []

    for tt = 1:nyr
        fnamein = frootin*string(years[tt])
        println("reading file "*fnamein)
        field = read_bin(fnamein,Float32,γ); fill!(field,0.0)

        push!(nseries,size(field,2))
        # may need to keep track of indices
        iend = istart+nseries[end]-1
        # iend = istart+nt-1

        field.f[:, :] .= flux_6hr_seasonal.f[:, istart:iend]
        istart = iend+1
        fnameout = fileout*string(years[tt])
        println("saved file "*fnameout)
        write(fnameout,field)
    end

end