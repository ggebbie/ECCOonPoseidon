#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings,
Plots

include(srcdir("config_exp.jl"))
include("./IsopycnalHelpers.jl")

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

region = "NPAC"; 
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

tecco= 1992+1/24:1/12:2018 # ecco years
E,F = trend_matrices(tecco)

#read in the first time step of S and θ

sig1grid = sigma1grid()
nσ = length(sig1grid)
# σlvls = findall( σtop .<= sig1grid .<= σbot)

#get root for salinity on sigma1
var = ["THETA", "theta", "theta", "theta"]
TSroot = "_on_sigma1" 
perc_good = zeros(nσ, 4)
all_pts = sum([length(area[ff]) for ff = 1:5])

#define exp 
for (i, expname) in enumerate(keys(shortnames))
    TSroot_ = var[i] * TSroot
    # Get list of files for salinity on sigma1
    filelist = searchdir(runpath[expname]*"sigma1/",TSroot_) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    nt = length(datafilelist)

    #load trend trend_matrices, F is LS estimator
    tt = 1
    println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
    Tname = datafilelist[tt]
    # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
    @time θσ1 = γ.read(runpath[expname]*"sigma1/"*Tname,MeshArray(γ,Float32,nσ))
    for k = 1:nσ
        good_pts = [0]
        for ff = 1:5
            good_pts .+= sum(-10 .< θσ1.f[ff, k] .< 50)
        end
        perc_good[k, i] = 100*(good_pts[1] ./ all_pts)
    end
end

p = plot([NaN], [NaN], label = nothing)
for (i, expname) in enumerate(keys(shortnames))
    plot!(p, sig1grid, perc_good[:, i], xlabel = "σ", ylabel = "%", 
    label = expname, 
    title = "Percentage of valid THETA points on year 1 month 1", legend=:topleft)
end
p