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
import NaNMath as nm

import PyPlot as plt
using PyCall
cm = pyimport("cmocean.cm");colorway = cm.balance;
@pyimport seaborn as sns;

sns.set_theme(context = "poster", style = "ticks", font_scale = 1.2,
              palette = sns.color_palette("colorblind"));
include(srcdir("config_exp.jl"))
include("./IsopycnalHelpers.jl")

(ϕ,λ) = latlonC(γ)

#reshape λ for plotting 
for ff ∈ [3,4] 
    λ[ff][λ[ff] .<= 0 ] = λ[ff][λ[ff] .<= 0 ] .+ 360
end


area = readarea(γ)

region = "PAC56"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")

runpath,diagpath = listexperiments(exprootdir());

tecco= 1992+1/24:1/12:2018 # ecco years
E,F = trend_matrices(tecco)


proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
cf = Vector{Any}(undef ,1)
#plotting climatology :-)
fig, ax = plt.subplots( figsize=(15,10), subplot_kw=Dict("projection"=> proj0))

ax.coastlines(resolution="110m")
ax.set_extent((-250, -70, -56, 70),crs=projPC)

gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
gl.xlabels_top = false;gl.ylabels_right = false

for ff in 1:5
    PAC_msk[ff][PAC_msk[ff] .≈ 0] .=NaN
    cf[1] = ax.pcolormesh(λ[ff], ϕ[ff],  PAC_msk[ff],
    vmin = 0, vmax = 2, shading="nearest", transform=projPC, 
    rasterized = true, cmap = "binary")               
end
# ax.set_title(labels[i])
ax.tick_params(bottom=true, left=true)
fig.tight_layout()
# fig.suptitle("Temperature Trends between 2-3 km depth on sigma1)")

fig.savefig(plotsdir(region * "_mask.png"))
fig
