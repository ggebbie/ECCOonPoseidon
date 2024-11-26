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

(!@isdefined(region)) && (region = "NPAC"); 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = OHC_helper.get_cell_volumes(area, cell_depths);

# masked_area = MeshArray(γ,Float32)
# for ff = 1:5
#     masked_area.f[ff] .= area.f[ff] .* PAC_msk.f[ff]
# end

#computed on standard pressures, 
#may have a different value if using pressures diagnosed from ECCO
#https://ecco-v4-python-tutorial.readthedocs.io/VectorCalculus_ECCO_barotropicVorticity.html
runpath,diagpath = listexperiments(exprootdir());

tecco= 1992+1/24:1/12:2018 # ecco years
E,F = trend_matrices(tecco)

#read in the first time step of S and θ
(!@isdefined(expname)) && (expname = "iter0_bulkformula"); 

sig1grid = sigma1grid()
nσ = length(sig1grid)
# σlvls = findall( σtop .<= sig1grid .<= σbot)

#get root for salinity on sigma1
#THETA for iter129
TSroot = "p_on_sigma1" 

# Get list of files for salinity on sigma1
filelist = searchdir(runpath[expname]*"sigma1/",TSroot) # first filter for state_3d_set1
datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
nt = length(datafilelist)

#load trend trend_matrices, F is LS estimator
P = zeros(Float32, nσ, nt);
println("\n obtaining timeseries for experiment " * expname * " in " * region * "\n")
σ1datadir(x="") = OHC_helper.σ1datadir(x)

@time for tt = 1:nt
    println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
    Tname = datafilelist[tt]
    # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
    @time Pσ1 = γ.read(runpath[expname]*"sigma1/"*Tname,MeshArray(γ,Float32,nσ))

    Pavg = area_average_by_σ(Pσ1, PAC_msk)
    P[:, tt] .= Pavg

end

jldsave(σ1datadir(expname * region * "_AVG_P_sigma1.jld2"), P= P)

p = contourf( tecco, sig1grid, P, yflip = true, 
color = :balance, levels = 51, ylim = (32, 33))
p