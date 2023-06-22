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
gr()
import NaNMath as nm
using .OHC_helper
import PyPlot as plt
using PyCall

@pyimport seaborn as sns;
@pyimport pandas as pd;
colors =  sns.color_palette("deep")[1:4]
sns.set_theme(context = "talk", font_scale = 1.0,
              palette = sns.color_palette("deep"));#sns.set_context("talk")
pygui(false)
cm = pyimport("cmocean.cm");colorway = cm.balance;
include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
tecco = 1992+1/24:1/12:2018


#load in latitude mask 
ocean_mask = wet_pts(Γ)
region = "NPAC30"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "full")
suffix = "2to3"
uplvl = -2e3; botlvl = -3e3
lvls = findall( botlvl .<= z[:].<= uplvl)

#create volume mask
area = readarea(γ)
cell_depths = get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 

#load in velocities 
κunadjasted_f = "/batou/eccodrive/files/Version4/Release4/input_init/total_diffkr_r009bit11.bin"
xx_κadjasted_f = "/batou/eccodrive/files/Version4/Release4/input_init/xx_diffkr.0000000129.data"
κunadjasted = read_bin(κunadjasted_f, MeshArray(γ,Float32,50))
xx_κadjasted = γ.read(xx_κadjasted_f,MeshArray(γ,Float32,50))

κunadjasted_masked = κunadjasted .* PAC_msk
κunadjasted_masked[findall(κunadjasted_masked .== 0)] = Inf

smush_depths = smush(cell_depths[:, lvls]); 
smush_depths[findall(smush_depths .==0 )] = Inf

adjustment_ratio_depth_avg = MeshArray(γ,Float32)
fill!(adjustment_ratio_depth_avg, 0.0)
for ff=1:5, k=lvls
    adjustment_ratio_depth_avg[ff] .+= (xx_κadjasted[ff, k] .* cell_depths[ff, k]) ./ smush_depths[ff]
end

for ff ∈ [3,4] 
    λ[ff][λ[ff] .<= 0 ] = λ[ff][λ[ff] .<= 0 ] .+ 360
end

pygui(true)
proj = ECCOonPoseidon.cartopy.crs.PlateCarree()
proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
cf = Vector{Any}(undef ,1)
fig, ax = plt.subplots(1, 1, figsize=(7,5), subplot_kw=Dict("projection"=> proj0))
b1 = minimum(adjustment_ratio_depth_avg)
b2 = maximum(adjustment_ratio_depth_avg)
bounds = [b1, b2]
b1, b2 = (-maximum(abs.(bounds)), maximum(abs.(bounds)))
# ax.set_title("Depths in " * region)
ax.set_global()
ax.coastlines()
ax.set_extent((110, -70, -0, 60))
ax.gridlines(crs=proj, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
for ff in 1:5
    cf[1] = ax.pcolormesh(λ[ff], ϕ[ff],  adjustment_ratio_depth_avg[ff],
    vmin = b1, vmax = b2, shading="nearest", transform=projPC, rasterized = true, cmap = colorway)            
end

cbar = fig.colorbar(cf[1], label = L"m", orientation = "horizontal")
