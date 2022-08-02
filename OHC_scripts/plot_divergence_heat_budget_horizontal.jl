using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools,
PyPlot, JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf, RollingFunctions
using PyPlot   # important!
using PyCall
using DataFrames, GLM

θ_depth_anom = remove_anomaly(θ_depths)
θ_depth_anom = θ_depths

fig, axs = plt.subplots(4, 1, figsize = (7, 10))
fig.suptitle(region * " θ Profile over model run")
tmat = repeat(tecco', length(lvls), 1)
zmat = repeat(z[lvls], 1, length(tecco))
extr = [Inf, -Inf]

for key in keys(shortnames)
    (minimum(θ_depth_anom[key]) < extr[1]) && (extr[1] = minimum(θ_depth_anom[key]))
    (maximum(θ_depth_anom[key]) > extr[2]) && (extr[2] = maximum(θ_depth_anom[key]))
end
clevs = collect(LinRange(extr[1], extr[2], 12))
cs = Any[0]
for (i, expname) in enumerate(keys(shortnames))
    cs[1] = axs[i].contourf(tmat, zmat, θ_depth_anom[expname],
    vmin = extr[1], vmax = extr[2], cmap = colorway, levels = clevs)
    axs[i].set_title(expname)
    axs[i].set_ylabel("depth")
end
axs[4].set_xlabel("time")

fig.subplots_adjust(right=0.8, hspace=0.4)
cbar_ax = fig.add_axes([0.83, 0.1, 0.05, 0.8])
cbar = fig.colorbar(cs[1], cax=cbar_ax,orientation="vertical")
cbar.set_label("ºC")

fig.savefig(plotsdir() * "/OHC_Divergence/" * "θDepthTS_" * region * suffix * ".png")

