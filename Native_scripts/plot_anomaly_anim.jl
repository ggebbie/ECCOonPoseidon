#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
Printf, LinearAlgebra, LaTeXStrings, PyCall, NaNStatistics
import PyPlot as plt
import NaNMath as nm
@pyimport matplotlib.animation as anim

cm = pyimport("cmocean.cm");colorway = cm.balance;

include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

region = "PAC56"; 
ocean_mask = OHC_helper.wet_pts(Γ)
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")

#computed on standard pressures, 
#may have a different value if using pressures diagnosed from ECCO
#https://ecco-v4-python-tutorial.readthedocs.io/VectorCalculus_ECCO_barotropicVorticity.html
runpath,diagpath = listexperiments(exprootdir());

#read in the first time step of S and θ
expname = "iter129_bulkformula"

# Get list of files for salinity on sigma1
filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
nt = length(datafilelist_θ)

tecco= 1992+1/24:1/12:2018 # ecco years

reg_mask = LLCcropC(PAC_msk,γ); reg_mask[reg_mask .== 0] .= NaN
reg_λ = LLCcropC(λ,γ); reg_ϕ = LLCcropC(ϕ,γ)
data = zeros(size(reg_ϕ)..., nt)

zlvls = findall( -2500 .<= z[:] .<= -2400)
# np = length(plvls)

for tt = 1:nt
    println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
    fnameθ = datafilelist_θ[tt]
    @time θz = γ.read(diagpath[expname]*fnameθ,MeshArray(γ,Float32,50))
    data[:,:, tt] .= LLCcropC(θz[:, zlvls],γ)
end

data_anom = 100 .* (data .- nm.mean(data))
bounds = nm.maximum(abs.(data_anom .* reg_mask))
proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

fig, ax = plt.subplots(figsize=(15,10), 
                       subplot_kw=Dict("projection"=> proj0))

z_approx = string(Int(floor(-z[zlvls][1])))

fig.suptitle("Temperature Anomalies in ECCO on z = " * z_approx)

cf = ax.pcolormesh(reg_λ, reg_ϕ,  data_anom[:,:, 1], 
                    vmin = -5, vmax = 5, shading="nearest", 
                    transform=projPC, rasterized = true, cmap = cm.balance) 
                    
fig.colorbar(cf, ax=ax, orientation = "horizontal", 
fraction=0.04, pad = 0.1, extend = "both", label = "[cK]")
ax.set_title("model_time")
gl = ax.gridlines(crs=projPC, draw_labels=true,
linewidth=2, color="gray", alpha=0, linestyle="--")
gl.top_labels = false; gl.right_labels = false

fig.tight_layout()

global fig, ax

function my_anim(tt)
    ax.clear()
    ax.coastlines()
    ax.set_extent((-326, -60, -60, 70),crs=projPC)
    model_time = string(Int(floor(tecco[tt+1]))) * " month " * string(((tt)%12)+1)
    ax.set_title(model_time)
    gl = ax.gridlines(crs=projPC, draw_labels=true,
    linewidth=2, color="gray", alpha=0, linestyle="--")
    gl.top_labels = false; gl.right_labels = false
    ax.pcolormesh(reg_λ, reg_ϕ,  data_anom[:,:, tt+1], 
                        vmin = -3, vmax = 3, shading="nearest", 
                        transform=projPC, rasterized = true, cmap = colorway) 
    ax.tick_params(bottom=true, left=true)
end

myanim = anim.FuncAnimation(fig, my_anim, frames=nt, interval=100)
savename = plotsdir("TempAnomalies_z" * z_approx * ".mp4")
@time myanim.save(savename)

