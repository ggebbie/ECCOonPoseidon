include("../../src/intro.jl")

using Revise, MAT
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools
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

lw = 2
E,F = trend_matrices(Float32.(tecco))
compute_depth_trends(x) = (100 * 100) .* (x * F[2, :])

sns.set_theme(context = "poster", style = "ticks",
              palette = colors, rc = custom_params);

fig = plt.figure(constrained_layout=true, figsize = (20, 12))
gs = fig.add_gridspec(4, 3)
obj(i,j) = get(gs, (i,j))
slice(i,j) = pycall(pybuiltin("slice"), PyObject, i,j)
f3_ax1 = fig.add_subplot(obj(slice(0, 4), 0))
f3_ax2 = fig.add_subplot(obj(slice(0, 2) ,slice(1, 4)), projection=proj0)
f3_ax3 = fig.add_subplot(obj(slice(2, 4) ,slice(1, 3)), projection=proj0)
ax = f3_ax1
# fig, ax = plt.subplots(1, 1, figsize=(9,12))
ax.vlines(0, 0, 5, color = "grey", alpha = 0.9)
plot_exps = ["iter129_bulkformula"]
label = ["ECCO V4r4"]
for (i, expt) in enumerate(plot_exps)
    println(plot_exps)
    ax.plot(compute_depth_trends(adjust_exps[expt]), z, label =  label[i], color = exp_colors[expt], linewidth = lw)
end

f =  jldopen(datadir("OPT-0015_GH19_" * region * "4.jld2"))
ΔT_GH19 = (-100 * 100 * f["ΔT_GH19"]) .- 0.5
ΔT_GH19 = ΔT_GH19 
depth_GH19 = f["depth_GH19"]
ax.plot(ΔT_GH19, depth_GH19, label = "OPT-15", color = "k", linewidth = 2)
ax.grid()
file = matopen(datadir("Challenger_WOCE_Temperature_basinwide_avg_FINAL.mat"))
z_Challenger = read(file, "depthlist")[:]
challenger_trend = -100. * 100 .* read(file, "mTpz_LS")[:] / 140
err = 100. * 100 .* read(file, "mTpz_LSerr")[:] ./ 140
ax.errorbar(challenger_trend[:], z_Challenger[:], xerr = err[:],
fmt ="o", label = "WOCE/Challenger", color = "darkblue", linewidth = 2)

rect = patches.Rectangle((-100*0.15, 2000), 100*0.3, 1000, linewidth=3, edgecolor="none",facecolor="black", alpha = 0.05)
ax.add_patch(rect)
ax.axvline(0, color = "black", zorder = 0, alpha = 0.55, linewidth = 2)
ax.set_title("Pacific Ocean\nTemperature Trend Profile")
ax.set_xlabel("cK per century"); ax.set_ylabel("Depth [m]")
ax.set_xlim(-1.5*10, 1.5*10)
ax.set_ylim(1001, 3999); ax.invert_yaxis()
ax.legend(frameon = true, fontsize = 18, loc = "lower right")
fig
axs = [f3_ax2, f3_ax3]
for ax in axs
    ax.coastlines(resolution="110m")
    ax.set_extent((110, 295, -65, 65),crs=projPC)
    gl = ax.gridlines(crs=projPC, draw_labels=true,
    linewidth=2, color="gray", alpha=0, linestyle="--")
    gl.top_labels = false
    gl.right_labels = false
end

fig.tight_layout()

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
cbar = axs[2].contourf(LONS, LATS,  -((β .* 100 * 100) .- 0.5), transform=projPC, 
    cmap = cmos.balance, vmin = -bounds, vmax = bounds, levels=levels, extend = "both")   
axs[2].set_title("OPT-15")
fig.subplots_adjust(wspace = 0.001)
fig.colorbar(cbar, ax = axs[:], label = "cK per century",
orientation = "horizontal", fraction = 0.03, pad = 0.07, extend = "both")
CS = axs[2].contour(LONS, LATS,  (β .* 100 * 100) .- 0.5, transform=projPC,     
c = "k", levels = [0])
axs[2].clabel(CS, CS.levels, inline=true, fontsize=17, inline_spacing = 10)
f3_ax1.annotate("A", (0.05, 0.92), fontsize = 40, 
xycoords="axes fraction", fontweight = "bold")

fig_labs = uppercase.(["b", "c", "d", "e"])
for (i, a) in enumerate([f3_ax2, f3_ax3])
    a.annotate(fig_labs[i], (0.05, 0.85), fontsize = 40, 
    xycoords="axes fraction", fontweight = "bold")
end
fig
fig.savefig(plotsdir("native/paper_figures/0.GH19_ECCO_Comparison.png"), dpi = 400)