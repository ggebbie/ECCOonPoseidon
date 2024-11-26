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

              
ocean_mask = OHC_helper.wet_pts(Γ)

tecco= 1992+1/24:1/12:2018 # ecco years
E,F = trend_matrices(tecco)
proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
cf = Vector{Any}(undef ,1)
#plotting climatology :-)

fig, ax = plt.subplots( figsize=(15,10), subplot_kw=Dict("projection"=> proj0))
ax.coastlines(resolution="110m")
ax.set_extent((-90, -0, -0, 60),crs=projPC)
gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
gl.top_labels = false; gl.right_labels = false



msk = region_mask(ocean_mask, λ, ϕ, (22, 37), (-73, -60))
for ff in 1:5
    msk[ff][msk[ff] .≈ 0] .=NaN

    cf[1] = ax.pcolormesh(λ[ff], ϕ[ff],  msk[ff],
    vmin = 0, vmax = 2, shading="nearest", transform=projPC, 
    rasterized = true, cmap = "binary")               
end
ax.tick_params(bottom=true, left=true)
fig.tight_layout()
fig

