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

#read in the first time step of S and θ
expname = "iter129_bulkformula"

include("./get_sigma1_levels.jl") 

sig1grid = sigma1grid()
nσ = length(sig1grid)
σlvls = findall( σtop .<= sig1grid .<= σbot)

#get root for salinity on sigma1
TSroot = "SALT_on_sigma1" 

# Get list of files for salinity on sigma1
filelist = searchdir(runpath[expname]*"sigma1/",TSroot) # first filter for state_3d_set1
datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
nt = length(datafilelist)
Ssig1slice = Array{Float32,1}(undef,nt)

# Find a latitudinal slice. i.e., 12 South. Use convert2array? Or something else? Use regular poles method.

# pre-allocate
#TS = Array{Float32, 2}(undef, nt, nz*2)

global tt = 0
tecco= 1992+1/24:1/12:2018 # ecco years

tt += 1
println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
Sname = datafilelist[tt]
# get S on sigma1. Way to read a slice? (Didn't figure it out yet)
@time S = γ.read(runpath[expname]*"sigma1/"*Sname,MeshArray(γ,Float32,nσ))

# just take one sigma1 surface
Ssig1 = S[:,σlvls]