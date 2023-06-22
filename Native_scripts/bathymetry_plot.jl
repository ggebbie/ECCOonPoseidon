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
cm = pyimport("cmocean.cm");colorway = cm.deep;
@pyimport seaborn as sns;

sns.set_theme(context = "talk", style = "ticks",
              palette = sns.color_palette("colorblind"));
include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)

area = readarea(γ)
ocean_mask = OHC_helper.wet_pts(Γ)
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")

cell_depths = Γ.Depth
# cell_depths[cell_depths .== 0] .= NaN
proj0 = ECCOonPoseidon.cartopy.crs.Robinson(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
cf = Vector{Any}(undef ,1)
λ_wrap = OHC_helper.wrap_λ(λ)

#plotting climatology :-)
fig, ax = plt.subplots( figsize=(15,10), subplot_kw=Dict("projection"=> proj0))
ax.coastlines(resolution="110m")
ax.set_extent((120, 285, -70, 70),crs=projPC)
gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
gl.xlabels_top = false;gl.ylabels_right = false
for ff in 1:5
    cell_depths[ff][cell_depths[ff] .≈ 0] .=NaN
    PAC_msk[ff][PAC_msk[ff] .≈ 0] .=NaN

    # cf[1] = ax.pcolormesh(λ_wrap[ff], ϕ[ff],  cell_depths[ff],
    # vmin = 0, vmax = 4500, shading="nearest", transform=projPC, 
    # rasterized = true, cmap = colorway)   

    cf[1] = ax.pcolormesh(λ_wrap[ff], ϕ[ff],  PAC_msk[ff], 
    shading="nearest", transform=projPC, alpha = 0.5,
    rasterized = true, cmap = "bwr_r")   

end
fig.savefig(plotsdir("native/generals/ECCO_NPAC_mask.png"), bbox_inches = "tight")

#plotting climatology :-)
fig, ax = plt.subplots( figsize=(15,10), subplot_kw=Dict("projection"=> proj0))
ax.coastlines(resolution="110m")
ax.set_extent((120, 285, -70, 70),crs=projPC)
gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
gl.xlabels_top = false;gl.ylabels_right = false
for ff in 1:5
    cell_depths[ff][cell_depths[ff] .≈ 0] .=NaN
    PAC_msk[ff][PAC_msk[ff] .≈ 0] .=NaN

    cf[1] = ax.pcolormesh(λ_wrap[ff], ϕ[ff],  cell_depths[ff],
    vmin = 0, vmax = 4500, shading="nearest", transform=projPC, 
    rasterized = true, cmap = colorway)   

    ax.pcolormesh(λ_wrap[ff], ϕ[ff],  PAC_msk[ff], 
    shading="nearest", transform=projPC, alpha = 0.5,
    rasterized = true, cmap = "bwr_r")   

end
ax.tick_params(bottom=true, left=true)
fig.colorbar(cf[1], ax = ax, orientation = "horizontal", fraction = 0.05, extend = "both")
fig.savefig(plotsdir("native/generals/ECCO_bathymetry_wNPAC.png"), bbox_inches = "tight")
fig


#plotting climatology :-)
fig, ax = plt.subplots( figsize=(15,10), subplot_kw=Dict("projection"=> proj0))
ax.coastlines(resolution="110m")
ax.set_extent((120, 285, -70, 70),crs=projPC)
gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
gl.xlabels_top = false;gl.ylabels_right = false
for ff in 1:5
    cell_depths[ff][cell_depths[ff] .≈ 0] .=NaN
    PAC_msk[ff][PAC_msk[ff] .≈ 0] .=NaN

    cf[1] = ax.pcolormesh(λ_wrap[ff], ϕ[ff],  cell_depths[ff],
    vmin = 0, vmax = 4500, shading="nearest", transform=projPC, 
    rasterized = true, cmap = colorway)   

end
ax.tick_params(bottom=true, left=true)
fig.colorbar(cf[1], ax = ax, orientation = "horizontal", fraction = 0.05, extend = "both")
fig.savefig(plotsdir("native/generals/ECCO_bathymetry.png"), bbox_inches = "tight")
fig
