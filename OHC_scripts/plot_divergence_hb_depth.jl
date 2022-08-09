using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools,
PyPlot, JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf, RollingFunctions
using PyPlot   # important!
using PyCall
using DataFrames, GLM
import Plots: Animation, frame, gif, cgrad
import Plots: contourf as jcontourf
@pyimport imageio
@pyimport matplotlib.animation as animation

function replacema(ma, vals)
    temp = ma
    for ff in 1:5
        temp[ff] .= replace(ma[ff], vals[1]=>vals[2])
    end
    return temp
end
nanmask = replacema(msk, 0=>NaN)
smush_depth = OHC_helper.smush(cell_depths) .* nanmask
for ff ∈ [3,4] 
    λ[ff][λ[ff] .<= 0 ] = λ[ff][λ[ff] .<= 0 ] .+ 360
end

proj = ECCOonPoseidon.cartopy.crs.PlateCarree()
proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
filenames = [plotsdir() * "/OHC_Divergence/temp"*string(tt)*".png"
             for tt in 1:length(tecco)]

cf = Vector{Any}(undef ,1)
fig, ax = plt.subplots(1, 1, figsize=(7,5), subplot_kw=Dict("projection"=> proj0))
ax.set_title("Depths in " * region)
# θ_flat["iter129_bulkformula"][tt][findall(θ_flat["iter129_bulkformula"][tt] .== 0)] .= NaN
b1 = nanmin(smush_depth)
b2 = nanmax(smush_depth)
ax.set_global()
ax.coastlines()
ax.set_extent((110, -70, -0, 60))
ax.gridlines(crs=proj, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
for ff in 1:5
        cf[1] = ax.pcolormesh(λ[ff], ϕ[ff],  smush_depth[ff],
        vmin = b1, vmax = b2, shading="nearest", transform=projPC, rasterized = true, cmap = colorway)            
end

cbar = fig.colorbar(cf[1], label = L"m", orientation = "horizontal")
tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/" * "θDepths_" * region * ".png",
            dpi = 1500)
close("all")