#  Remove the interannual frequency energy in surface forcing fields.
#  Keep the interannual energy in a given region.
#  Diagnostic plots have been removed from this version.
#
#  Script argument: region (must be defined in `src/ECCOonPoseidon.jl`)
#  If no arguments are passed, then interannual variability is removed everywhere.

include("../src/intro.jl")

using Revise
using ECCOtour
using ECCOonPoseidon
using Statistics
using PyPlot
using Distributions
using FFTW
using LinearAlgebra
using StatsBase
using MeshArrays
using MITgcmTools
using MAT

# from intro.jl, default is nointerannual
expt == "nointerannual" ? keepregion = false : keepregion = true
println("Experiment: ",expt)

include(srcdir("config_exp.jl"));

include(srcdir("config_regularpoles.jl"))

basin1 = basin_mask("Pacific",γ)

# clash with `mask` name
if keepregion
    # plot_basinmask
    #msk = basin_mask(ϕC,λC,"Pacific","N")
    msk = basin_mask("Pacific",γ)

    msk_regpoles = var2regularpoles(msk,γ,nx,ny,nyarc,λarc,nyantarc,λantarc)

    figure()
    clf()
    cmap_seismic =get_cmap("seismic")
    lims = range(0.0,step=0.05,stop=1.0)
    contourf(λCregpoles,ϕCregpoles,msk_regpoles',lims,cmap=cmap_seismic)
    colorbar(label="weight",orientation="vertical",ticks=lims)
    mkpath(plotsdir())
    outfname = plotsdir("southpac_mask.eps")
    xlbl = "longitude "*L"[\degree E]"
    ylbl = "latitude "*L"[\degree N]"
    titlelbl = expt*" mask"
    title(titlelbl)
    xlabel(xlbl)
    ylabel(ylbl)
    savefig(outfname)
end


# This could be put into src code for scientific project.
inputdir = fluxdir()

# lets use a test directory so we don't overwrite something
outputdir = fluxdir(expt)*"test"

# read lat, lon at center of grid cell
(ϕC,λC) = latlonC(γ)
# on the vector (Staggered) grid
(ϕG,λG) = latlonG(γ)
 
## Define region of interest where interannual fluxes are kept.
# region = "test"
# if keepregion
#     latrect, lonrect = rectangle(region)

#     # should dlat, dlon go into src?
#     dlat = 10
#     dlon = 10
    
#     lonmid =  (lonrect[1]+lonrect[2])/2
#     centerlon!(λC,lonmid)
#     centerlon!(λG,lonmid)
# end

# if !isdir(outputdir)
#     mkpath(outputdir)
# end
midname = "_6hourlyavg_"
varnames = ("atmPload","oceFWflx","oceQsw","oceSflux","oceSPflx","oceTAUE","oceTAUN","oceTAUX",
            "oceTAUY","sIceLoad","sIceLoadPatmPload","sIceLoadPatmPload_nopabar","TFLUX")

# sample calcs at one point. Get the length of timeseries in nseries.
yv = 40
xv = 180
fv = 5
tmplat  = ϕC[fv]; lat_point = tmplat[xv,yv]
tmplon  = λC[fv]; lon_point = tmplon[xv,yv]
frootsample = inputdir*varnames[end]*midname
years = 1992:2017
fluxsample_point,nseries = extract_timeseries(frootsample,years,γ,xv,yv,fv)

# take biweekly mean. Use triangular filter with break points at:
Δi14day = 4*14 # grid index range
nt6hr = length(fluxsample_point)
i14day = 2.5:Δi14day:nt6hr+56 # goes past end of time by 14 days to be sure

t6hr_start = 1/8 # 3Z Jan 1 1992
Δt6hr = 6/24 # 1/4 of a day
t6hr = range(t6hr_start,step=Δt6hr,length=nt6hr)

# Given a full-resolution timeseries, what values at the tiepoints best represent the timeseries?
# Determine tiepoints (in time) where fluxes are adjusted by ECCO in optimization procedure.
Δt14day = 14 # units: days
t14day_start = 1/2 # 12Z Jan 1 1992 by inspection of figure
t14day = range(t14day_start,step=Δt14day,stop=t6hr[end]+14)# goes past the end by one
nt14day = length(t14day)

# Values are added to the t14day tiepoints, then linearly interpolated to fill gaps.
# Careful not to store high-resolution data all at same time.
# Make a function that takes values at the tiepoints and then makes a full-resolution timeseries.
@time E14to6,F6to14 = get_filtermatrix(t6hr,t14day)

daysperyear = 365.25
fcycle = 1/(daysperyear) # units: day^{-1}
# for removing seasonal cycle from 14-day averaged timeseries
Ecycle,Fcycle = seasonal_matrices(fcycle,t14day)

# interannual filter is a Hann(ing) filter
Thann = 100.0 # days

#vname = varnames[1] # for interactive use
for vname ∈ varnames
    filein = inputdir*vname*midname
    fileout = outputdir*vname*midname
    println(filein)

    if keepregion
        # some tracers have a staggered grid.
        # will have to re-do the mask unfortunately to be sure.
        if vname == "oceTAUX"
            msk = regional_mask(ϕC,λG,latrect,lonrect,dlat,dlon)
        elseif vname == "oceTAUY"
            msk = regional_mask(ϕG,λC,latrect,lonrect,dlat,dlon)
        else
            msk = regional_mask(ϕC,λC,latrect,lonrect,dlat,dlon)
        end
    end
    
    # use this F to decompose controllable/uncontrollable parts
    # i.e., 6 hourly to 14 day
    flux_14day = matrixfilter(F6to14,filein,years,γ)

    # remove seasonal cycle from 14-day averaged timeseries
    # solve for seasonal cycle parameters
    # use 14-day because it's computationally efficient
    βcycle = matmul(Fcycle,flux_14day,γ)

    # reconstruct the full seasonal cycle.
    flux_14day_seasonal = matmul(Ecycle,βcycle,γ)

    # remove seasonal signal from total signal
    flux_14day_noseasonal = flux_14day - flux_14day_seasonal

    # take the noseasonal timeseries and run a filter (Thann [days] cutoff)
    # hanning filter for all locations.
    # output is interannual signal.
    flux_14day_lopass = hannfilter(flux_14day_noseasonal,t14day,t14day,Thann,γ)

    #use function that takes forcing field and multiplies each entry by a number
    #from 0 to 1 (1 being removed entirely) based on a provided lat/lon box and sponge layer widths
    # pre-compute mask before looping over all variables.
    if keepregion
        apply_regional_mask!(flux_14day_lopass,1.0 .- msk)
    end

    # put tflux_14day_lopass on to 6hr
    # check for NaN's in output
    nancount_lopass = sum(nancount(flux_14day_lopass))
    nancount_14day = sum(nancount(flux_14day))
    nancount_14day_seasonal = sum(nancount(flux_14day_seasonal))
    nancount_14day_noseasonal = sum(nancount(flux_14day_noseasonal))

    if nancount_lopass + nancount_14day + nancount_14day_seasonal + nancount_14day_noseasonal > 0
        error("NaNs in the filtered output")
    end

    # write the hi-pass filtered tflux
    # need to get it on the 6hr timesteps
    # need to write it for every year
    matrixspray(E14to6,-flux_14day_lopass,filein,fileout,years,γ) #changed rmfield to _regional

end

# MV DIAGS TO DIFFERENT SCRIPT/FUNCTION
# make a figure to see a spatial map of flux, oceTAUN, oceTAUE are most interesting
diags = false
if diags # this section requires updating to paths
    filename1 = "/poseidon/ecco.jpl.nasa.gov/drive/files/Version4/Release4/other/flux-forced/forcing/oceTAUN_6hourlyavg_2000"
    field1 = read_bin(filename1,Float32,γ)
    filename2 = "/poseidon/ecco.jpl.nasa.gov/drive/files/Version4/Release4/other/flux-forced-interannual_southpac/oceTAUN_6hourlyavg_2000"
    field2 = read_bin(filename2,Float32,γ)

    # translate to regularpoles
    field_regpoles =  var2regularpoles(field2[:,100]-field1[:,100],γ,nx,ny,nyarc,farc,iarc,jarc,warc,nyantarc,fantarc,iantarc,jantarc,wantarc)

    figure()
    clf()
    lims = range(-0.1,step=0.005,stop=0.1)
    contourf(λCregpoles,ϕCregpoles,field_regpoles',lims,cmap=cmap_seismic)
    colorbar(label="wind stress",orientation="vertical",ticks=lims)
    outfname = outputdir*"delta_oceTAUN-2000-1.eps"
    xlbl = "longitude "*L"[\degree E]"
    ylbl = "latitude "*L"[\degree N]"
    titlelbl = L"\tau_y, "*"interannual_southpac - iter129, yr 2000, time 100, "*L"[N/m^2]"
    title(titlelbl)
    xlabel(xlbl)
    ylabel(ylbl)
    savefig(outfname)
end
