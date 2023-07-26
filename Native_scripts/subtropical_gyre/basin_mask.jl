#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings,
    PyCall
using .OHC_helper
import PyPlot as plt
@pyimport seaborn as sns;
@pyimport pandas as pd;
sns.set_theme(context = "talk", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));

include(srcdir("config_exp.jl"))

#using . is faster 
(ϕ,λ) = latlonC(γ)
area = readarea(γ)

include("SargassoMask.jl")

tecco= 1992+1/24:1/12:2018; nz = 50

proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

#plotting climatology :-)
fig, ax = plt.subplots( figsize=(15,10), subplot_kw=Dict("projection"=> proj0))
ax.coastlines(resolution="110m")
ax.set_extent((-90, -20, 0, 60),crs=projPC)
gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
gl.top_labels = false; gl.right_labels = false

# remove_seasonal_array(x) = hcat([mean(x[:, (i):(i+ 23)], dims = 2) for i in 1:24:312]...)

for ff in 1:5
    msk[ff][msk[ff] .≈ 0] .=NaN
    ax.pcolormesh(λ[ff], ϕ[ff],  msk[ff],
    vmin = 0, vmax = 2, shading="nearest", transform=projPC, 
    rasterized = true, cmap = "binary")               
end

ax.tick_params(bottom=true, left=true)
fig.tight_layout()
fig

face, index, _ = OHC_helper.findlatlon(λ, ϕ, -65, 30);
ax.scatter(-65, 30, transform=projPC)
fig

fig.savefig(plotsdir("native/subtropical_gyre/SargassoMask.png"))
fig
