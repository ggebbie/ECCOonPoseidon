#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings,
Plots

include(srcdir("config_exp.jl"))
include("../IsopycnalHelpers.jl")

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

region = "PAC56"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = OHC_helper.get_cell_volumes(area, cell_depths);

H = OHC_helper.smush(cell_depths); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H; inv_H = Float32.(inv_H)

#computed on standard pressures, 
#may have a different value if using pressures diagnosed from ECCO
#https://ecco-v4-python-tutorial.readthedocs.io/VectorCalculus_ECCO_barotropicVorticity.html
runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

tecco= 1992+1/24:1/12:2018 # ecco years
E,F = trend_matrices(tecco)

#read in the first time step of S and θ

sig1grid = sigma1grid()
nσ = length(sig1grid)
# σlvls = findall( σtop .<= sig1grid .<= σbot)

#get root for salinity on sigma1
var = ["THETA", "theta", "theta", "theta"]
TSroot = "_on_sigma1" 
avg_temps = zeros(nσ, 4)
all_pts = sum([length(area[ff]) for ff = 1:5])

#define exp 
mask_area = PAC_msk .* area

for (i, expname) in enumerate(keys(shortnames))
    TSroot_ = var[i] * TSroot
    # Get list of files for salinity on sigma1
    filelist = searchdir(runpath[expname]*"sigma1/",TSroot_) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    nt = length(datafilelist)

    tt = 1
    println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
    Tname = datafilelist[tt]
    # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
    @time θσ1 = γ.read(runpath[expname]*"sigma1/"*Tname,MeshArray(γ,Float32,nσ))
    θavg = area_average_by_σ(θσ1, mask_area)
    avg_temps[:, i] .= θavg
end

p = plot([NaN], [NaN], label = nothing)
for (i, expname) in enumerate(keys(shortnames))
    plot!(p, avg_temps[:, i], sig1grid, xlabel = "T [°C]", ylabel = "σ",
    label = expname, yflip = true,
    title = "Pacific Ocean THETA on isopycnal surfaces \n year 1 month 1", legend=:bottomright)
end
p

