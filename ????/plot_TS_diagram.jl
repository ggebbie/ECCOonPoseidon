#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings,
Plots, PyCall
import NaNMath as nm
import PyPlot as plt
include(srcdir("config_exp.jl"))
include("./IsopycnalHelpers.jl")
cm = pyimport("cmocean.cm");colorway = cm.balance;

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
expname = "iter129_bulkformula"

sig2grid = sigma2grid()
nσ = length(sig2grid)
# σlvls = findall( σtop .<= sig1grid .<= σbot)

#get root for salinity on sigma1
#THETA for iter129
TSroot = "THETA_on_sigma2" 

# Get list of files for salinity on sigma1
filelist = searchdir(runpath[expname]*"sigma2/",TSroot) # first filter for state_3d_set1
datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
nt = length(datafilelist); 

#load trend trend_matrices, F is LS estimator
θ = zeros(Float32, nσ, nt);
println("\n obtaining timeseries for experiment " * expname * " in " * region * "\n")

for ff = 1:5
    PAC_msk[ff][PAC_msk[ff] .== 0.0] .= NaN
end

tt = 1
println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
Tname = datafilelist[tt]
# get S on sigma1. Way to read a slice? (Didn't figure it out yet)
@time θσ2 = γ.read(runpath[expname]*"sigma2/"*Tname,MeshArray(γ,Float32,nσ))
θ_zonal = zonal_avg(θσ2 .* PAC_msk)

bounds = nm.extrema(θ_zonal)
levels = bounds[1]:0.5:bounds[2]
Y = sig2grid
X = OHC_helper.ma_zonal_sum(ϕ .* area) ./ OHC_helper.ma_zonal_sum(area)

fig, ax = plt.subplots(1,1, figsize = (20, 12.5))
CS = ax.contourf(X, Y, θ_zonal, vmin = bounds[1], vmax = bounds[2], levels = levels, 
extend = "both", cmap = cm.balance);
ax.invert_yaxis()
ax.set_xticks(-50:10:60)
ax.set_xlim(-49, 60)
ax.set_xlabel("latitude"); ax.set_ylabel("kg per m³")
fig.colorbar(CS, orientation = "horizontal", fraction = 0.05, 
label = L"[°C ]")
ax.set_title("Zonally Averaged Temperature in Isopycnal Coordinates")
fig