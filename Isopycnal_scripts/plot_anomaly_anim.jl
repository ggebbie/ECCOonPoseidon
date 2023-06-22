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
include("./IsopycnalHelpers.jl")

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
expname = "iter0_bulkformula"

sig1grid = sigma1grid()
nσ = length(sig1grid)

#get root for salinity on sigma1
#THETA for iter129
TSroot = "THETA_on_sigma1" 

# Get list of files for salinity on sigma1
filelist = searchdir(runpath[expname]*"sigma1/",TSroot) # first filter for state_3d_set1
datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
nt = length(datafilelist)

tecco= 1992+1/24:1/12:2018 # ecco years

reg_mask = LLCcropC(PAC_msk,γ); reg_mask[reg_mask .== 0] .= NaN

reg_λ = LLCcropC(λ,γ); reg_ϕ = LLCcropC(ϕ,γ)

data = zeros(size(reg_ϕ)..., nt)

fname = expname * region * "_AVG_P_sigma1.jld2" 
!ispath(σ1datadir(fname)) && (include("get_avgP_on_sigma1.jl"))
Pσ = load(σ1datadir(fname), "P")
sigma2p = mean(Pσ, dims = 2)
p2z = gsw_z_from_p.(abs.(sigma2p), 30., 0., 0.)[:]
plvls = findall( -2600 .<= p2z .<= -2450)
# np = length(plvls)

for tt = 1:nt
    # θsig1reg =  zeros(size(reg_ϕ)..., np)

    println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
    Tname = datafilelist[tt]
    @time θσ1 = γ.read(runpath[expname]*"sigma1/"*Tname,MeshArray(γ,Float32,nσ))

    # θsig1 = θσ1[:, plvls] #get a specific level
    # [θsig1reg[:, :, k] .= LLCcropC(θsig1[:, k],γ) for k in 1:np]
    # data[:,:, tt] .= nanmean(θsig1reg, dims=3)
    data[:,:, tt] .= LLCcropC(θσ1[:, plvls[1]],γ)
end

baseline = nm.mean(data[:, :, 1] .* reg_mask)
data_anom = 100 .* (data .- baseline)

bounds = 1
proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

fig, ax = plt.subplots(figsize=(15,10), 
                       subplot_kw=Dict("projection"=> proj0))

σ_approx = string(sig1grid[plvls][1])
z_approx = string(Int(floor(-p2z[plvls[1]])))

fig.suptitle("Temperature Anomalies in ECCO on σ = "* 
σ_approx * " (z ≈ " * z_approx * ")")

cf = ax.pcolormesh(reg_λ, reg_ϕ,  data_anom[:,:, 1], 
                    vmin = -bounds, vmax = bounds, shading="nearest", 
                    transform=projPC, rasterized = true, cmap = colorway) 
fig.colorbar(cf, ax=ax, orientation = "horizontal", 
fraction=0.04, pad = 0.1, extend = "both", label = "[K]")
model_time = string(Int(floor(tecco[1]))) * " month " * string(((0)%12)+1)
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
                        vmin = -bounds, vmax = bounds, shading="nearest", 
                        transform=projPC, rasterized = true, cmap = colorway) 
    ax.tick_params(bottom=true, left=true)
end

myanim = anim.FuncAnimation(fig, my_anim, frames=nt, interval=100)

σ_approx = string(Int(floor(sig1grid[plvls][1])))

savename = plotsdir("TempGlobAnomalies_σ" * σ_approx * "_z" * z_approx * expname * ".gif")
@time myanim.save(savename)

