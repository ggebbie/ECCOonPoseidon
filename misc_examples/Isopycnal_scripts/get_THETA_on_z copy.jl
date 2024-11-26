#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")
using .OHC_helper

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings,
Plots
include(srcdir("config_exp.jl"))
include("./IsopycnalHelpers.jl")

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

(!@isdefined(region)) && (region = "PAC56"); 

PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = OHC_helper.get_cell_volumes(area, cell_depths);
smush_depths = OHC_helper.smush(cell_depths); smush_depths[findall(smush_depths .== 0)] = Inf
inv_depths = 1 ./ smush_depths

#computed on standard pressures, 
#may have a different value if using pressures diagnosed from ECCO
#https://ecco-v4-python-tutorial.readthedocs.io/VectorCalculus_ECCO_barotropicVorticity.html
runpath,diagpath = listexperiments(exprootdir());

tecco= 1992+1/24:1/12:2018 # ecco years
E,F = trend_matrices(tecco)

#read in the first time step of S and θ
(!@isdefined(expname)) && (expname = "iter0_bulkformula"); 

# Get list of files for salinity on sigma1
filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

nt = length(datafilelist_θ)
nz = length(z)
#load trend trend_matrices, F is LS estimator
θ = zeros(Float32, nz, nt);

println("\n obtaining timeseries for experiment " * expname * " in " * region * "\n")
@time for tt = 1:nt
    println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
    # fnameS = datafilelist_S[tt]
    fnameθ = datafilelist_θ[tt]
    @time θ_gcm = γ.read(diagpath[expname]*fnameθ,MeshArray(γ,Float32,50))
    #  sθ = OHC_helper.extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_depths)

    θavg = volume_average_by_depth(θ_gcm, cell_volumes, γ)
    θ[:, tt] .= vec(θavg)
end

jldsave(σ1datadir(expname * region * "_AVG_THETA_z.jld2"), θ= θ)

θanom = θ .- mean(θ, dims = 2)
p = contourf( tecco, abs.(z), θanom, yflip = true, 
color = :balance, clims=(-1.5,1.5))
p