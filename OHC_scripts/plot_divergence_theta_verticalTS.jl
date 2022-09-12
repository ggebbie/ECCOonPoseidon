#this function creates a contour plot, with 
#depth on the y-Axis, time on the x-axis and contour values 
#as temperature values 
using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf, RollingFunctions,
LaTeXStrings

#finds the min and max of all ts
function plot_zonal_contours!(X, Y, var_zonal, axs, cs, cmap)
    extremas = map(x-> nm.extrema(x), values(var_zonal))
    extr = [minimum(minimum.(extremas)), 
            maximum(maximum.(extremas))]
    clevs = collect(LinRange(extr..., 15))
    for (i, exp) in enumerate(keys(var_zonal))
        cs[1] = axs[i].contourf(X, Y, var_zonal[exp],
        vmin = extr[1], vmax = extr[2], cmap = cmap, levels = clevs)
        axs[i].set_title(exp)
        axs[i].set_xlabel("time")
    end
    for i=1:2
        axs[i].set_ylabel("depth [m]")
    end
end

# θ_depth_anom = remove_anomaly(θ_depths)
X = repeat(tecco', length(lvls), 1); Y = repeat(z[lvls], 1, length(tecco))
cs = Vector{Any}(missing, 1)

θ_depth_anom = θ_depths
fig, axs = plt.subplots(2, 2, figsize = (18, 12))
fig.suptitle(region * " θ Profile over experiment run")
plot_zonal_contours!(X, Y, θ_depth_anom, axs, cs, cm.balance)
cbar_ax = fig.add_axes([0.1, 0.1, 0.8, 0.05]);cbar = fig.colorbar(cs[1], cax=cbar_ax,orientation="horizontal")
cbar.set_label("C");fig.subplots_adjust(bottom=0.22, hspace=0.4)
fig.savefig(plotsdir() * "/OHC_Divergence/Contours/" * "θVertTrend_" * region * "_" * suffix * ".png",
dpi = 1500)

θ_depth_anom = Dict(key => θ_depths[key] .- θ_depths["iter0_bulkformula"] 
                    for key in keys(θ_depths))
fig, axs = plt.subplots(2, 2, figsize = (18, 12))
fig.suptitle(region * " θ Profile over experiment run")
plot_zonal_contours!(X, Y, θ_depth_anom, axs, cs, cm.balance)
cbar_ax = fig.add_axes([0.1, 0.1, 0.8, 0.05]);cbar = fig.colorbar(cs[1], cax=cbar_ax,orientation="horizontal")
cbar.set_label("C");fig.subplots_adjust(bottom=0.22, hspace=0.4)
fig.savefig(plotsdir() * "/OHC_Divergence/Contours/" * "θVertTrend_iter0anom_" * region * "_" * suffix * ".png",
dpi = 1500)

close("all")