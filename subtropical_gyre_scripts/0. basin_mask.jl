#this script just plots out a rectangular mask 

include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings,
    PyCall

    import PyPlot as plt
@pyimport seaborn as sns;
@pyimport pandas as pd;
sns.set_theme(context = "talk", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));

include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

include("SargassoMask.jl")

tecco= 1992+1/24:1/12:2018; nz = 50

proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

#plotting climatology :-)
fig, ax = plt.subplots( figsize=(15,10), subplot_kw=Dict("projection"=> proj0))
ax.coastlines(resolution="110m")
ax.set_extent((-85, -70, 20, 40),crs=projPC)
gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
gl.top_labels = false; gl.right_labels = false

for ff in 1:5
    msk[ff][msk[ff] .≈ 0] .=NaN
    ax.pcolormesh(λ[ff], ϕ[ff],  msk[ff],
    vmin = 0, vmax = 2, shading="nearest", transform=projPC, 
    rasterized = true, cmap = "binary")               
end

ax.tick_params(bottom=true, left=true)
fig.tight_layout()
ax.scatter(-80, 26, transform=projPC)
fig