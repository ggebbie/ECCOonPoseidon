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
filelist = searchdir(runpath[expname]*"sigma1/",TSroot) # first filter for state_3d_set1
datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
nt = length(datafilelist)

tecco= 1992+1/24:1/12:2018 # ecco years

reg_mask = LLCcropC(PAC_msk,γ)
reg_mask[reg_mask .== 0] .= NaN
reg_λ = LLCcropC(λ,γ); reg_ϕ = LLCcropC(ϕ,γ)

data = zeros(size(reg_ϕ)..., nt)


plvls = findall( -1500 .<= z[:] .<= -1000)
nzlev = length(plvls)

for tt = 1:nt
    θsig1reg =  zeros(size(reg_ϕ)..., np)

    println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
    Tname = datafilelist[tt]
    @time θσ1 = γ.read(runpath[expname]*"sigma1/"*Tname,MeshArray(γ,Float32,nσ))

    θsig1 = θσ1[:, plvls] #get a specific level
    [θsig1reg[:, :, k] .= LLCcropC(θsig1[:, k],γ) for k in 1:np]
    data[:,:, tt] .= nanmean(θsig1reg, dims=3)
end

data_anom = 100 .* (data .- mean(data, dims = 3))
bounds = nm.maximum(abs.(data_anom .* reg_mask))
proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
cf = Vector{Any}(undef ,1)
fig, ax = plt.subplots(figsize=(15,10), 
                       subplot_kw=Dict("projection"=> proj0))
fig.suptitle("Temperature Anomalies [cK] in ECCO between 1-1.5km")

cf = ax.pcolormesh(reg_λ, reg_ϕ,  data_anom[:,:, 1], 
                    vmin = -3, vmax = 3, shading="nearest", 
                    transform=projPC, rasterized = true, cmap = colorway) 
fig.colorbar(cf, ax=ax)
model_time = string(Int(floor(tecco[1]))) * " month " * string(((0)%12)+1)
ax.set_title(model_time)
fig.tight_layout()

global fig, ax

function my_anim(tt)
    ax.clear()
    ax.coastlines()
    model_time = string(Int(floor(tecco[tt+1]))) * " month " * string(((tt)%12)+1)
    ax.set_title(model_time)
    ax.pcolormesh(reg_λ, reg_ϕ,  data_anom[:,:, tt+1], 
                        vmin = -3, vmax = 3, shading="nearest", 
                        transform=projPC, rasterized = true, cmap = colorway) 
end

myanim = anim.FuncAnimation(fig, my_anim, frames=nt, interval=100)
@time myanim.save("test5.mp4")
