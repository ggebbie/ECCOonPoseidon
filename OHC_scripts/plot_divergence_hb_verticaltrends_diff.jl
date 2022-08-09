using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools,JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf, RollingFunctions
using PyCall
using DataFrames, GLM

# θ_depth_anom = remove_anomaly(θ_depths)
θ_depth_anom = θ_depths

fig, axs = plt.subplots(2, 2, figsize = (12, 6))
fig.suptitle(region * " θ Anomaly Profile over experiment run \n Baseline: iter0")
tmat = repeat(tecco', length(lvls), 1)
zmat = repeat(z[lvls], 1, length(tecco))
extr = [Inf, -Inf]
for key in keys(shortnames)
    baseline = θ_depth_anom["iter0_bulkformula"] 
    var = θ_depth_anom[key] .- baseline
    (minimum(var)  < extr[1]) && 
    (extr[1] = minimum(var))
    (maximum(var) > extr[2]) && 
    (extr[2] = maximum(var))
end
extr .= [-maximum(abs.(extr)), maximum(abs.(extr))]
clevs = collect(LinRange(extr[1], extr[2], 15))

cs = Any[0]
for (i, expname) in enumerate(keys(shortnames))
    baseline = θ_depth_anom["iter0_bulkformula"] 
    cs[1] = axs[i].contourf(tmat, zmat, θ_depth_anom[expname] .- baseline,
    vmin = extr[1], vmax = extr[2], cmap = colorway, levels = clevs)
    axs[i].set_title(expname)
end
axs[4].set_xlabel("time")
for i=1:2
    axs[i].set_ylabel("depth [m]")
end
fig.subplots_adjust(bottom=0.2, hspace=0.4)
cbar_ax = fig.add_axes([0.1, 0.1, 0.8, 0.05])
cbar = fig.colorbar(cs[1], cax=cbar_ax,orientation="horizontal")
cbar.set_label("C")

fig.savefig(plotsdir() * "/OHC_Divergence/" * "θVertTrendAnom_" * region * suffix * ".png",
dpi = 1500)
close("all")