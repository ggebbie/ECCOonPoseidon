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
import NaNMath as nm
using PyPlot   # important!

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
nanmin(X) = minimum([nm.minimum(x) for x in X])
nanmax(X) = maximum([nm.maximum(x) for x in X])
nanextr(X) = (nanmin(X), nanmax(X))
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

function vmean(x::MeshArray)
    sm = similar(x[:, 1]); fill!(sm, 0.0)
    nz = size(x, 2)
    
    for ff =1:5
        cnt = zeros(Int, size(sm[ff])...)
        for k=1:nz
            xmsk = copy(x[ff, k])
            cnt .+= (!isnan).(xmsk)
            xmsk[isnan.(xmsk)] .= 0.0 
            sm[ff] .+= xmsk
        end
        sm[ff] = sm[ff] ./ cnt
        sm[ff][cnt .== 0] .= NaN
    end
    return sm
end

for ff ∈ [3,4] 
    λ[ff][λ[ff] .<= 0 ] = λ[ff][λ[ff] .<= 0 ] .+ 360
end

proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

fig, axes = plt.subplots(3, 3, figsize=(14,10), subplot_kw=Dict("projection"=> proj0))
cf = Vector{Any}(undef ,1)
fig.suptitle("THETA at σ₁")
expname  = "iter129_bulkformula"; i =1

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

#remove bad values and outside of region
for ij in eachindex(θσ1)
    bad_mask = Int.((θσ1[ij] .< -50 ) .+ (θσ1[ij] .> 50))
    θσ1[ij][bad_mask .== 1] .= NaN
end
lvls = collect(0:9:81)

for j = 1:9
    lower, upper = (lvls[j]+1, lvls[j+1])
    ax = axes[j]
    ax.set_title("σ₁ = " * string( round(mean(sig1grid[lower:upper]); sigdigits=4 )))
    maps = vmean(θσ1[:, lower:upper]).* PAC_msk
    # vmin, vmax = nanextr(maps)
    vmin, vmax = (0, 25)

    println(vmin, " ", vmax)
    ax.set_global(); ax.coastlines()
    ax.set_extent((110, -70, -60, 60))
    for ff in 1:5
        maps[ff][PAC_msk[ff] .== 0] .= NaN
        cf[1] = ax.pcolormesh(λ[ff], ϕ[ff],  maps[ff], shading="nearest", 
        transform=projPC, rasterized = true, vmin = vmin, vmax = vmax)            
    end
    cbar = fig.colorbar(cf[1], ax = ax,
    label = "[°C]", orientation = "vertical")
end

fig.tight_layout()
fig