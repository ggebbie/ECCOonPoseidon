#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings
import NaNMath as nm

import PyPlot as plt
using PyCall
cm = pyimport("cmocean.cm");colorway = cm.balance;
@pyimport seaborn as sns;

sns.set_theme(context = "notebook", style = "ticks", font_scale = 1.2,
              palette = sns.color_palette("colorblind"));
include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)

area = readarea(γ)

              
region = "NPAC"; 
ocean_mask = OHC_helper.wet_pts(Γ)
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
bathy_mask = Γ.Depth .> 2000; bathy_mask = Float32.(bathy_mask)
runpath,diagpath = listexperiments(exprootdir());

tecco= 1992+1/24:1/12:2018 # ecco years
E,F = trend_matrices(tecco)
proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
cf = Vector{Any}(undef ,1)
#plotting climatology :-)
fig, ax = plt.subplots( figsize=(15,10), subplot_kw=Dict("projection"=> proj0))
ax.coastlines(resolution="110m")
ax.set_extent((-250, -65, -56, 60),crs=projPC)
gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
gl.top_labels = false; gl.right_labels = false
msk = PAC_msk .* bathy_mask

OHC_helper.get_min_lat(ϕ,  msk)
OHC_helper.get_max_lat(ϕ,  msk)


for ff in 1:5
    msk[ff][msk[ff] .≈ 0] .=NaN
    cf[1] = ax.pcolormesh(λ[ff], ϕ[ff],  msk[ff],
    vmin = 0, vmax = 2, shading="nearest", transform=projPC, 
    rasterized = true, cmap = "binary")               
end
ϕ[4][80:150, :]
ax.tick_params(bottom=true, left=true)
fig.tight_layout()
# fig.savefig("mask.png")
fig




abs_dist(x, r) = abs(x) < r
ϕ_mask_min = Float32(OHC_helper.get_min_lat(ϕ, PAC_msk)); ϕ_mask = ϕ .> -Inf
[ϕ_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- ϕ_mask_min, 0.1) .* PAC_msk.f[ff] for ff in 1:5]
ϕ_mask = Float32.(ϕ_mask);ϕ_mask_min

fig, ax = plt.subplots( figsize=(15,10), subplot_kw=Dict("projection"=> proj0))
ax.coastlines(resolution="110m")
ax.set_extent((-250, -70, -56, 70),crs=projPC)
gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
gl.xlabels_top = false;gl.ylabels_right = false
for ff in 1:5
    ϕ_mask[ff][ϕ_mask[ff] .≈ 0] .= NaN
    cf[1] = ax.pcolormesh(λ[ff], ϕ[ff],  ϕ_mask[ff],
    vmin = 0, vmax = 2, shading="nearest", transform=projPC, 
    rasterized = true, cmap = "binary")               
end
ax.tick_params(bottom=true, left=true)
fig.tight_layout()
fig

PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not", include_bering = false)
bathy_mask = Γ.Depth .> 2000; bathy_mask = Float32.(bathy_mask)
abs_dist(x, r) = abs(x) < r

ϕ_mask_max = 51.5; ϕ_max_mask = ϕ .> -Inf
[ϕ_max_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- ϕ_mask_max, 0.2) .* PAC_msk.f[ff] for ff in 1:5]
# ϕ_max_mask = ϕ_max_mask .* within_lon(λ, 152, 170)
ϕ_max_mask = Float32.(ϕ_max_mask)
sum(ϕ_max_mask)
# tmp=-tmp

PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not", include_bering = true)
bathy_mask = Γ.Depth .> 2000; bathy_mask = Float32.(bathy_mask)
bathy_mask = bathy_mask .* PAC_msk

fig, ax = plt.subplots( figsize=(15,10), subplot_kw=Dict("projection"=> proj0))
ax.coastlines(resolution="110m")
ax.set_extent((-220, -120, 20, 60),crs=projPC)
gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
gl.xlabels_top = false;gl.ylabels_right = false

for ff in 1:5
    bathy_mask[ff][bathy_mask[ff] .≈ 0] .= NaN
    ϕ_max_mask[ff][ϕ_max_mask[ff] .≈ 0] .= NaN
    
    cf[1] = ax.pcolormesh(λ[ff], ϕ[ff],  bathy_mask[ff],
    vmin = 0, vmax = 1, shading="nearest", transform=projPC, 
    rasterized = true, cmap = "bwr")         

    cf[1] = ax.pcolormesh(λ[ff], ϕ[ff],  ϕ_max_mask[ff],
    vmin = 0, vmax = 2, shading="nearest", transform=projPC, 
    rasterized = true, cmap = "binary", alpha = 0.5)         
end


ax.tick_params(bottom=true, left=true)
fig.tight_layout()

