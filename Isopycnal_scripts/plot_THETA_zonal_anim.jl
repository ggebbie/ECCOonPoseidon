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
import NaNMath as nm
import PyPlot as plt
include(srcdir("config_exp.jl"))
include("./IsopycnalHelpers.jl")

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

region = "PAC56"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = OHC_helper.get_cell_volumes(area, cell_depths);

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
TSroot = "THETA_on_sigma1" 

# Get list of files for salinity on sigma1
filelist = searchdir(runpath[expname]*"sigma1/",TSroot) # first filter for state_3d_set1
datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
nt = length(datafilelist); 

#load trend trend_matrices, F is LS estimator
θ = zeros(Float32, nσ, nt);
println("\n obtaining timeseries for experiment " * expname * " in " * region * "\n")
σ1datadir(x="") = OHC_helper.σ1datadir(x)

for ff = 1:5
    PAC_msk[ff][PAC_msk[ff] .== 0.0] .= NaN
end

θ_zonal = zeros(nσ, 270, nt)

@time for tt = 1:nt
    println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
    Tname = datafilelist[tt]
    # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
    @time θσ1 = γ.read(runpath[expname]*"sigma1/"*Tname,MeshArray(γ,Float32,nσ))
    θ_zonal[:, :, tt] .= zonal_avg(θσ1 .* PAC_msk)
end

θ_zonal_z = θ_zonal[40:end, :, :]
Y = sig1grid[40:end]
X = OHC_helper.ma_zonal_sum(ϕ .* area) ./ OHC_helper.ma_zonal_sum(area)

θ_zonal_anom = θ_zonal_z .- mean(θ_zonal_z, dims =3)

fig, ax = plt.subplots(1,1, figsize = (20, 12.5))
CS = ax.contourf(X, Y, θ_zonal_anom[:, : ,10], extend = "both");
ax.invert_yaxis()
ax.set_xticks(-40:10:60)
ax.set_xlim(-39, 60)
fig.colorbar(CS, orientation = "horizontal", fraction = 0.05, 
label = L"[mm^2 / s]")
fig