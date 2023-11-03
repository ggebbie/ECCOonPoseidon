include("../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools
import PyPlot as plt

include(srcdir("plot_and_dir_config.jl"))
@pyimport matplotlib.patches as patches

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)
region = "NPAC"; 
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
compute_depth_trends(x) = (100 * 10) .* (x * F[2, :])

fig, ax = plt.subplots(figsize = (10, 10))
plot_exps = ["iter0_bulkformula", "iter129_bulkformula", "only_wind", "only_kappa", "only_init"]
for expt in plot_exps
    println(plot_exps)
    ax.plot(compute_depth_trends(adjust_exps[expt]), z, label = plot_labels[expt], color = exp_colors[expt], linewidth = lw)
end

f =  jldopen(datadir("OPT-0015_GH19_NPAC.jld2"))
ΔT_GH19 = 100 * 10 * f["ΔT_GH19"]
depth_GH19 = f["depth_GH19"]

ax.plot(ΔT_GH19, depth_GH19, label = "OPT-15", color = "k", linewidth = lw)

rect = patches.Rectangle((-10*0.1, 2000), 10*0.2, 1000, linewidth=3, edgecolor="none",facecolor="black", alpha = 0.05)
ax.add_patch(rect)
ax.axvline(0, color = "black", zorder = 0, alpha = 0.55, linewidth = 1)
ax.legend(frameon = false)
ax.set_title("North Pacific Temperature Trend Profile", fontweight = "bold")
ax.set_xlabel("cK per decade", fontweight = "bold"); ax.set_ylabel("Depth [m]", fontweight = "bold")
ax.set_xlim(-1, 1)
ax.set_ylim(1000, 5000); ax.invert_yaxis()
fig
# fig.savefig(plotsdir("native/sensitivity_exps/North_Pacific_Temp_Trend.png"))
