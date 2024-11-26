include("../../../../src/intro.jl")

using Revise, MAT
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools, NCDatasets
import PyPlot as plt

include(srcdir("plot_and_dir_config.jl"))
@pyimport matplotlib.patches as patches
@pyimport cmocean.cm as cmos
@pyimport matplotlib.tri as tri

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)
region = "PAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)

cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)
ΔV = lateral_sum(cell_volumes)

lvls = findall( -3000 .<= -z[:].<= -2000)
mid_depths(x) = vec(sum(Float32, x[lvls, :] .* ΔV[lvls], dims = 1) / sum(Float32, ΔV[lvls]))

tecco = 1992+1/24:1/12:2018; nz = 50
get_datafiles(expname, key) = filter(x -> occursin("data",x),searchdir(diagpath[expname],key) )

adjust_exps =  jldopen(datadir(region * "_temperature_sens_exps.jld2"))["adjust_exps"]

lw = 2.5
E,F = trend_matrices(Float32.(tecco))
compute_depth_trends(x) = (100 * 100) .* (x * F[2, :])

include(srcdir("plot_and_dir_config.jl"))

sns.set_style("darkgrid", Dict("axes.facecolor" => ".98", "grid.color" => "0.5", "grid.alpha" => 0.3))

ds_obs  = NCDataset(datadir("Pacific_Ocean_Profiles_Trends.nc"))["obs_trends"][:]
ds_obs_depth  = NCDataset(datadir("Pacific_Ocean_Profiles_Trends.nc"))["depth"][:]

fig, ax = plt.subplots(figsize = (8., 10))
ax.tick_params(which="both", bottom=true, left = true)

#Plot ECCO V4r4
plot_exps = ["iter129_bulkformula"]
label = ["ECCO V4r4"]
for (i, expt) in enumerate(plot_exps)
    println(plot_exps)
    ax.plot(compute_depth_trends(adjust_exps[expt]), z, label =  label[i], color = exp_colors[expt], linewidth = 5.0, zorder = 10)
end

# ax.plot(ds_obs, ds_obs_depth, label =  "WOCE/GO-SHIP\n(1992 - 2017)", color = exp_colors["iter129_bulkformula"], linewidth = lw, linestyle = "--", alpha = 0.5, zorder = 0)

file = matopen(datadir("Challenger_WOCE_Temperature_basinwide_avg_FINAL.mat"))
z_Challenger = read(file, "depthlist")[:]
challenger_trend = -100. * 100 .* read(file, "mTpz_LS")[:] / 115
err = 100. * 100 .* read(file, "mTpz_LSerr")[:] ./ 115
ax.errorbar(challenger_trend[:], z_Challenger[:], xerr = err[:],
fmt ="o", label = "WOCE/Challenger", color = "darkblue", linewidth = 2.5, alpha = 0.3, zorder = 0)

#Plot OPT-15
f =  jldopen(datadir("OPT-0015_GH19_PAC_FINAL.jld2"))
ΔT_GH19 = (100 * 100 * f["ΔT_GH19"]) 
ΔT_GH19 = ΔT_GH19 
depth_GH19 = f["depth_GH19"]
ax.plot(ΔT_GH19, depth_GH19, label = "GH19 Reconstruction", color = "darkblue", linewidth = 5.0, zorder = 10)

# rect = patches.Rectangle((-100*0.15, 2000), 100*0.3, 1000, linewidth=3, edgecolor="none",facecolor="black", alpha = 0.05)
# ax.add_patch(rect)
# ax.axvline(0, color = "black", zorder = 0, alpha = 0.55, linewidth = 2)
ax.set_title("Pacific Ocean Temperature Trends")
ax.set_xlabel("cK per century"); ax.set_ylabel("Depth [m]", fontweight = "bold")
ax.set_xlim(-1.5*10, 1.5*10)
ax.set_ylim(1001, 3999); ax.invert_yaxis()
ax.legend(loc = "lower left", frameon = false, bbox_to_anchor = (0.00, 0-.2), ncols = 2)
# ax.vlines(0, 0, 5000, color = "grey", alpha = 1, linestyle = "--", zorder = 0, lw = 3.5)
ax.grid(alpha = 0.3)
fig.savefig(plotsdir("wind_rewrite/0.GH19_ECCO_Trend_Comparison.png"), dpi = 400, bbox_inches = "tight")
fig

fig, axs = plt.subplots(1, 2, figsize = (10, 8), subplot_kw=Dict("projection"=> proj0))

for (i, ax) in enumerate(axs)
    ax.coastlines(resolution="110m", color = "#949494")
    ax.set_extent((110, 295, -65, 65),crs=projPC)
    gl = ax.gridlines(crs=projPC, draw_labels=true,
    linewidth=2, color="gray", alpha=0, linestyle="--")
    gl.top_labels = false
    gl.right_labels = false
    if i == 2
        gl.left_labels = false
    end
    ax.add_feature( ECCOonPoseidon.cartopy.feature.LAND, facecolor="#949494")

end

fig.tight_layout()
fig
############ Making 
# fig, axs = plt.subplots(2, 1, figsize=(17,12), subplot_kw=Dict("projection"=> proj0))
bounds = 0.14 * 100
levels = collect(-bounds:2:bounds)

fname = datadir("native/_THETA_spatial_trend_" * "2to3" * ".jld2")
β =  jldopen(fname)["theta_trends"]["FULL"]
data = β .* 100 * 100
x = vcat([λ.f[ff][:] for ff = 1:5]...);  x[x .< 0] .+= 360
y = vcat([ϕ.f[ff][:] for ff = 1:5]...)
z = vcat([data.f[ff][:] for ff = 1:5]...)
mask = 1 .* isnan.(z)
triang = tri.Triangulation(x, y)
mask = isnan.(z[triang.triangles[:, 2] .+ 1] .+ z[triang.triangles[:, 1] .+ 1] .+ z[triang.triangles[:, 3] .+ 1])
triang.set_mask(mask)
z[isnan.(z)] .= 0
cf = axs[1].tricontourf(triang, z, transform=projPC,     
cmap = cmos.balance, levels = levels, 
vmin = -bounds, vmax = bounds, extend = "both")
CS = axs[1].tricontour(triang, z, transform=projPC,     
c = "k", levels = [0])
axs[1].clabel(CS, CS.levels, inline=true, fontsize=17, inline_spacing = 10)
axs[1].set_title("ECCO V4r4")

fname = "modern_OPT-0015_θ_trends_2to3km_2.jld2"
β = jldopen(datadir(fname))["β"]
LONS = jldopen(datadir(fname))["λ"]
LATS = jldopen(datadir(fname))["ϕ"]
cbar = axs[2].contourf(LONS, LATS,  -((β .* 100 * 100)), transform=projPC, 
    cmap = cmos.balance, vmin = -bounds, vmax = bounds, levels=levels, extend = "both")   
axs[2].set_title("GH19 Reconstruction (OPT-15)")
# fig.subplots_adjust(wspace = 0.001)
fig.colorbar(cbar, ax = axs[:], label = "cK per century",
orientation = "horizontal", fraction = 0.03, pad = 0.07, extend = "both")
CS = axs[2].contour(LONS, LATS,  (β .* 100 * 100) .- 0.5, transform=projPC,     
c = "k", levels = [0])
axs[2].clabel(CS, CS.levels, inline=true, fontsize=17, inline_spacing = 10)


fig_labs = uppercase.(["a", "b", "d", "e"])
for (i, a) in enumerate(axs)
    a.annotate(fig_labs[i], (0.05, 0.85), fontsize = 25, 
    xycoords="axes fraction", fontweight = "bold")
end

fig.savefig(plotsdir("wind_rewrite/0.GH19_ECCO_Comparison.png"), dpi = 400, bbox_inches = "tight")
fig