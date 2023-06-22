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
import PyPlot as plt
using PyCall
import NaNMath as nm
cm = pyimport("cmocean.cm");colorway = cm.balance;

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

#read in the first time step of S and θ
expname = "iter129_bulkformula"

sig1grid = sigma1grid()
nσ = length(sig1grid)
# σlvls = findall( σtop .<= sig1grid .<= σbot)

#get root for salinity on sigma1
#THETA for iter129
TSroot = "THETA_on_sigma1" 

# Get list of files for salinity on sigma1
filelist = searchdir(runpath[expname]*"sigma1/",TSroot) # first filter for state_3d_set1
datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
nt = length(datafilelist)

#reshape \lambda for plotting 
for ff ∈ [3,4] 
    λ[ff][λ[ff] .<= 0 ] = λ[ff][λ[ff] .<= 0 ] .+ 360
end

proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
cf = Vector{Any}(undef ,1)

tecco= 1992+1/24:1/12:2018 # ecco years
E,F = trend_matrices(tecco)
β = MeshArray(γ,Float32); fill!(β, 0.0)

level = 70
@time for tt = 1:nt
    println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
    Tname = datafilelist[tt]
    # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
    @time θσ1 = γ.read(runpath[expname]*"sigma1/"*Tname,MeshArray(γ,Float32,nσ))
    θσ1_sel = θσ1[:, level]
    for ff in 1:5
        β[ff] .+= F[2,tt] .* θσ1_sel[ff]
    end
    #can use this to form a proper regular grid or EOF
    # #       
    # Ssig1 = S[:,isig1grid]

    # # next slice at 12 South.
    # Ssig1crop =  LLCcropC(Ssig1,γ) # regular grid by cropping     
end

#plotting climatology :-)
fig, ax = plt.subplots( figsize=(15,10), subplot_kw=Dict("projection"=> proj0))

bounds = maximum([nm.maximum(β[ff].*PAC_msk[ff]) for ff = 1:5])
β_mask = deepcopy(β); β_mask[findall(β_mask.==0)] = NaN
ax.coastlines(resolution="110m")
ax.set_extent((-250, -80, -40, 70),crs=projPC)

gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
gl.xlabels_top = false;gl.ylabels_right = false

for ff in 1:5
    cf[1] = ax.pcolormesh(λ[ff], ϕ[ff],  β_mask[ff],
    vmin = -bounds, vmax = bounds, shading="nearest", transform=projPC, 
    rasterized = true, cmap = colorway)               
end
# ax.set_title(labels[i])
ax.tick_params(bottom=true, left=true)

fig.tight_layout()
# fig.suptitle("Temperature Trends between 2-3 km depth)")
fig
