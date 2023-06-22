#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
Printf, LinearAlgebra, LaTeXStrings,Plots, PyCall
import PyPlot as plt
import NaNMath as nm
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

proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
cf = Vector{Any}(undef ,1)

tecco= 1992+1/24:1/12:2018 # ecco years
E,F = trend_matrices(tecco)

regular_mask = LLCcropC(PAC_msk,γ)
regular_mask[regular_mask .== 0] .= NaN
regular_λ = LLCcropC(λ,γ)
regular_ϕ = LLCcropC(ϕ,γ)

data = zeros(length(regular_ϕ), nt)
@time for tt = 1:nt
    println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
    Tname = datafilelist[tt]
    # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
    @time θσ1 = γ.read(runpath[expname]*"sigma1/"*Tname,MeshArray(γ,Float32,nσ))

    θsig1 = θσ1[:, 70]

    Ssig1crop =  LLCcropC(θsig1,γ) .* regular_mask # regular grid by cropping     
    data[:, tt] .= vec(Ssig1crop)
end

#naive fill
#better approach is found here https://stackoverflow.com/a/35611142 
#remove mean and ignore points with NaNs
data .-= mean(data, dims = 2)
wherenan = isnan.(data); data[wherenan] .= 0

U, S, V = svd(data, full = false)
eig = S.^2 ./ sum(S.^2)

fig, ax = plt.subplots(figsize=(15,10), 
                       subplot_kw=Dict("projection"=> proj0))

#plot EOF1
EOFn = reshape(U[:, 1], 360, 180)
EOFn[vec(isnan.(regular_mask))] .= NaN

ax.coastlines(resolution="110m"); 
ax.tick_params(bottom=true, left=true)
ax.set_extent((-250, -80, -40, 70),crs=projPC)

gl = ax.gridlines(crs=projPC, draw_labels=true,
                  linewidth=2, color="gray", 
                  alpha=0, linestyle="--")
gl.top_labels = false; gl.right_labels = false

vmi, vma = nm.maximum(abs.(EOFn)) .* (-1, 1)

ax.pcolormesh(regular_λ, regular_ϕ,  EOFn, 
              vmin = vmi, vmax = vma, shading="nearest", 
              transform=projPC, rasterized = true, cmap = colorway) 
fig

#comparing EOF1 
regular_A = LLCcropC(area .* PAC_msk,γ)
#temperature timeseries
weighted_avg(x, weights) = sum(vec(weights) .* x, dims = 1) ./ sum(weights)
θ_anom = weighted_avg(data, regular_A) 

#EOF1 timeseries 
EOFn = reshape(U[:, 1], 360, 180)
βs = U[:, 1]' * data
θ_anom_EOF1 = weighted_avg(U[:, 1] .* βs, regular_A) 

fig, ax = plt.subplots(figsize=(15,10))
ax.plot(tecco, θ_anom[1, :], label = "data")
ax.plot(tecco, θ_anom_EOF1[1, :], label = "EOF1")
ax.legend()
fig

