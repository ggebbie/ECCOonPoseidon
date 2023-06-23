include("./src/intro.jl")

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

# plot filtered oceTAUE minus original oceTAUE
expt = "interannual_northpac"
include(srcdir("config_exp.jl"));
include(srcdir("config_regularpoles.jl"))

year = "2005"
inputdir = fluxdir()
outputdir = fluxdir(expt)*"test/"
midname = "_6hourlyavg_"
vname = "oceTAUN"
filein = joinpath(inputdir,vname*midname)
fileout = joinpath(outputdir,vname*midname)
    #filename1 = "/poseidon/ecco.jpl.nasa.gov/drive/files/Version4/Release4/other/flux-forced/forcing/oceTAUN_6hourlyavg_2000"
    #filename1 = "/batou/eccodrive/files/Version4/Release4/other/flux-forced/forcing/oceTAUN_6hourlyavg_2000"
field1 = read_bin(filein*year,Float32,γ)
    #filename2 = "/poseidon/ecco.jpl.nasa.gov/drive/files/Version4/Release4/other/flux-forced-interannual_southpac/oceTAUN_6hourlyavg_2000"
    #filename2 = "/batou/eccodrive/files/Version4/Release4/other/flux-forced-nointerannual/forcing/test/oceTAUN_6hourlyavg_2000"
field2 = read_bin(fileout*year,Float32,γ)

# translate to regularpoles
field_regpoles = var2regularpoles(field2[:,100]-field1[:,100],γ,nx,ny,nyarc,λarc,nyantarc,λantarc)

figure()
clf()
lims = range(-0.1,step=0.005,stop=0.1)
#contourf(λC,ϕC,field_regpoles',lims,cmap=cmap_seismic)
contourf(field_regpoles',lims)
#colorbar(label="wind stress",orientation="vertical",ticks=lims)
outfname = outputdir*"/testdelta_oceTAUN-"*year*".eps"
xlbl = "longitude "*L"[\degree E]"
ylbl = "latitude "*L"[\degree N]"
titlelbl = L"\tau_y, "*"interannual_northpac - iter129, yr "*year*", time 100, "*L"[N/m^2]"
title(titlelbl)
xlabel(xlbl)
ylabel(ylbl)
savefig(outfname)


#look at timeseries --------------------------------------------------------------
#look for place where mask is true
# read lat, lon at center of grid cell
(ϕC,λC) = latlonC(γ)
# on the vector (Staggered) grid
(ϕG,λG) = latlonG(γ)
fv = 4 #1 #grid face number

outfname = outputdir*"testtimeseries_"*vname*".eps"

#look at map of timeseries differences
yv = 1:20:81; 
xv = 1:25:226;
forcingdiffs = zeros(length(xv),length(yv))
for yy = 1:length(yv)
    for xx = 1:length(xv)
        tmplat  = ϕC[fv]; lat_point = tmplat[xv[xx],yv[yy]]
        tmplon  = λC[fv]; lon_point = tmplon[xv[xx],yv[yy]]
        frootsample = inputdir*vname*midname
        newrootsample = outputdir*vname*midname
        years = 1992:2017
        origsample_point,nseries = extract_timeseries(frootsample,years,γ,xv[xx],yv[yy],fv)
        newsample_point,nseries = extract_timeseries(newrootsample,years,γ,xv[xx],yv[yy],fv)
        sample_diff = newsample_point - origsample_point;
        forcingdiffs[xx,yy] = mean(abs.(sample_diff))
    end
end

figure()
contourf(forcingdiffs)
colorbar()
savefig(outputdir*"testtimeseries_difference_mean_map"*vname*".eps")


#look at one timeseries
#pick a lon and lat
yv = 90 #200
xv = 40 #20
fv = 4 #1 #grid face number

tmplat  = ϕC[fv]; lat_point = tmplat[xv,yv]
tmplon  = λC[fv]; lon_point = tmplon[xv,yv]
frootsample = inputdir*vname*midname
newrootsample = outputdir*vname*midname
years = 1992:2017
origsample_point,nseries = extract_timeseries(frootsample,years,γ,xv,yv,fv)
newsample_point,nseries = extract_timeseries(newrootsample,years,γ,xv,yv,fv)

figure()
plot(origsample_point)
plot(newsample_point,linestyle =":")
savefig(outfname)

figure()
plot(newsample_point - origsample_point)
savefig(outputdir*"testtimeseries_difference_"*vname*".eps")




#clf()
field1 = read_bin(filein*year,Float32,γ)
field2 = read_bin(fileout*year,Float32,γ)

field1_regpoles = zeros(360,275,size(field1,2))
field2_regpoles = zeros(360,275,size(field1,2))
for tt = 1:size(field1,2)
    field1_regpoles[:,:,tt] = var2regularpoles(field1[:,tt],γ,nx,ny,nyarc,λarc,nyantarc,λantarc)
    field2_regpoles[:,:,tt] = var2regularpoles(field2[:,tt],γ,nx,ny,nyarc,λarc,nyantarc,λantarc)
end

figure()
plot(field1_regpoles[longitude,latitude,:])
plot(field2_regpoles[longitude,latitude,:])
