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
inputdir = fluxdir()
outputdir = fluxdir("testing_climatology")
# read lat, lon at center of grid cell
(ϕC,λC) = latlonC(γ)
# on the vector (Staggered) grid

(!isdir(outputdir))&&(mkpath(outputdir))

midname = "_6hourlyavg_"
# varnames = ("atmPload","oceFWflx","oceQsw","oceSflux","oceSPflx","oceTAUE","oceTAUN","oceTAUX",
#             "oceTAUY","sIceLoad","sIceLoadPatmPload","sIceLoadPatmPload_nopabar","TFLUX")
varnames = ["oceTAUX"]

frootsample = inputdir*varnames[end]*midname
# sample calcs at one point. Get the length of timeseries in nseries.
yv = 10
xv = 70
fv = 4
tmplat  = ϕC[fv]; lat_point = tmplat[xv,yv]
tmplon  = λC[fv]; lon_point = tmplon[xv,yv]

years = 1992:1995

fluxsample_point,nseries = extract_timeseries(frootsample,years,γ,xv,yv,fv)
                        # take biweekly mean. Use triangular filter with break points at:
nt6hr = length(fluxsample_point)

t6hr_start = 1/8 # 3Z Jan 1 1992
Δt6hr = 6/24 # 1/4 of a day
t6hr = range(t6hr_start,step=Δt6hr,length=nt6hr)

daysperyear = 365.25
fcycle = 1/(daysperyear) # units: day^{-1}
# for removing seasonal cycle from 14-day averaged timeseries
Ecycle,Fcycle = seasonal_matrices(fcycle,t6hr)
Edag = Ecycle * Fcycle #is symmetric! 

#vname = varnames[1] # for interactive use
@time for vname ∈ varnames
    nyr = length(years)

    filein = inputdir*vname*midname
    fileout = outputdir*vname*midname
    println(filein)

    nout = size(Edag,1)
    nin  = size(Edag,2)

    tmp = MeshArray(γ,Float32,nout) # some nans her
    tmp = convert2array(tmp)
    flux_6hr_seasonal = zeros(Float32, (nout, size(tmp)...)) #initialize as array
    # reconstruct the full seasonal cycle.

    istart = 1
    nseries = []

    for tt = 1:nyr
        fnamein = filein*string(years[tt])
        println("reading file "*fnamein)
        field = read_bin(fnamein,Float32,γ); 
        push!(nseries,size(field,2))
        field_array = zeros(Float32, ( nseries[end], size(tmp)...)) #initialize as array
        println(size(field_array))
        println(size(field))
        println(size(field))

        for i = 1:nseries[end]
            field_array[i, :, :] .= convert2array(field[:, i])
        end
        println("done converting")
        # may need to keep track of indices
        iend = istart+nseries[end]-1
        E_dag_sub = @view(Edag[:, istart:iend])
        for jj in 1:nout
            weights = reshape(@view(E_dag_sub[jj, :]), (1, 1, nseries[end]))
            weighted_6hr = sum(Float32, field_array .* weights, dims = 3)
            flux_6hr_seasonal[jj, :, :] .+= weighted_6hr[:, :, 1]
        end
        istart = iend+1

    end

    # print("seasonal cycle obtained")
    # # put tflux_14day_lopass on to 6hr
    # # check for NaN's in output
    # nancount_6hr_seasonal = sum(nancount(flux_6hr_seasonal))

    # if nancount_6hr_seasonal > 0
    #     error("NaNs in the filtered output")
    # end

    # istart = 1

    # for (i, tt) = enumerate(1:nyr)
    #     fnamein = filein*string(years[tt])
    #     println("reading file "*fnamein)
    #     field = read_bin(fnamein,Float32,γ); fill!(field,0.0)

    #     # may need to keep track of indices
    #     iend = istart+nseries[i]-1
    #     fnameout = fileout*string(years[tt])
    #     println("saved file "*fnameout)
    #     write(fnameout,flux_6hr_seasonal[:, istart:iend])
    #     istart = iend+1

    # end

end

vname = varnames[1]
actual,nseries = extract_timeseries(inputdir*vname*midname,years,γ,xv,yv,fv)
seasonal_saved,nseries = extract_timeseries(outputdir*vname*midname,years,γ,xv,yv,fv)

nt = length(actual)
fig, ax = plt.subplots(1, figsize=(15, 5))
ax.plot(1:nt, actual, label = "Original TFLUX [Δt=6 hours]"); ax.set_xlabel("years")
ax.plot(1:nt, seasonal_saved, linewidth = 4, 
label = "Seasonal Cycle of TFLUX [Δt=6 hours]");

println(mean(actual .- seasonal_saved))
ax.legend()
fig


filein = inputdir*vname*midname
fnamein = filein*string(years[1])
println("reading file "*fnamein)
field = read_bin(fnamein,Float32,γ); 

field[5]

flux_6hr_seasonal = MeshArray(γ,Float32,size(Edag, 2)) # some nans here