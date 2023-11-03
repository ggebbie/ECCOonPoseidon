include("../../../src/intro.jl")

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
lw = 3
fig, ax = plt.subplots(figsize = (10, 7.5))
plot_exps = ["iter0_bulkformula", "only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]
for expt in plot_exps 
    println(plot_exps)
    ax.plot(tecco, mid_depths(adjust_exps[expt]), label = plot_labels[expt], color = exp_colors[expt], linewidth = lw)
end
ax.legend(frameon = false)
ax.set_xlabel("time", fontweight = "bold"); ax.set_ylabel(" θ [°C]", fontweight = "bold")
fig
fig.savefig(plotsdir("native/sensitivity_exps/0.North_Pacific_Sens_Exps_T.png"), bbox_inches = "tight")

Δadjust_exps = deepcopy(adjust_exps)
Dexps = ["only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]
[Δadjust_exps[nexp] .-= adjust_exps["iter0_bulkformula"] for nexp in Dexps]
Δadjust_exps["SUM"] = Δadjust_exps["only_init"] .+ Δadjust_exps["only_kappa"] .+ 
                      Δadjust_exps["only_wind"] .+ Δadjust_exps["only_buoyancy"]


fig, ax = plt.subplots(figsize = (10, 7.5))
plot_exps = ["only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula", "SUM"]
exp_colors["SUM"] = "grey"; plot_labels["SUM"] = "SUM"
for expt in plot_exps 
    println(plot_exps)
    ax.plot(tecco, 100 .* mid_depths(Δadjust_exps[expt]), label = plot_labels[expt], color = exp_colors[expt], linewidth = lw)
end
# ax.axhline(0, c = "k",= "--")
ax.legend(frameon = false)
ax.set_xlabel("time", fontweight = "bold"); ax.set_ylabel("Δθ [cK]", fontweight = "bold")
ax.set_ylim(-0.8, 0.8)
fig
fig.savefig(plotsdir("native/sensitivity_exps/0.North_Pacific_Sens_Exps_ΔT.png"), bbox_inches = "tight")

E,F = trend_matrices(Float32.(tecco))
compute_depth_trends(x) = (100 * 10) .* (x * F[2, :])

adjust_exps_trends = Dict()
adjust_exps_trends = Dict(key => compute_depth_trends(adjust_exps[key]) for key in keys(adjust_exps))
exps = ["Initial", "Kappa", "Forcing"]
[adjust_exps_trends[key] = adjust_exps_trends[key] .- adjust_exps_trends["CTRL"] for key in exps] #remove the effect of the baseline
adjust_exps_trends["SUM"] = adjust_exps_trends["CTRL"] .+ adjust_exps_trends["Initial"] .+ 
                            adjust_exps_trends["Kappa"] .+ adjust_exps_trends["Forcing"]

fig, ax = plt.subplots(figsize = (10, 10))
ax.plot(adjust_exps_trends["CTRL"], z, label = "CTRL (Iteration 0)")
ax.plot(adjust_exps_trends["Initial"], z, label = "INIT Effect")
ax.plot(adjust_exps_trends["Kappa"], z, label = "MIXING Effect")
ax.plot(adjust_exps_trends["Forcing"], z, label = "FORCING Effect")
ax.plot(adjust_exps_trends["SUM"], z, label = "SUM", color = "k", linestyle = "--", alpha = 0.6)
ax.plot(adjust_exps_trends["FULL"], z, label = "FULL (Iteration 129)", color = "k", linewidth = 3)
rect = patches.Rectangle((-10*0.1, 2000), 10*0.2, 1000, linewidth=3, edgecolor="none",facecolor="black", alpha = 0.05)
ax.add_patch(rect)
ax.axvline(0, color = "black", zorder = 0, alpha = 0.55, linewidth = 1)
ax.legend(frameon = false)
ax.set_title("North Pacific Temperature Trend Profile", fontweight = "bold")
ax.set_xlabel("cK per decade", fontweight = "bold"); ax.set_ylabel("Depth [m]", fontweight = "bold")
ax.set_xlim(-1, 1)
ax.set_ylim(1500, 5000); ax.invert_yaxis()
fig
fig.savefig(plotsdir("native/sensitivity_exps/North_Pacific_Temp_Trend.png"))
# fig.savefig(plotsdir("native/"), bbox_inches = "tight")

compute_depth_diff(x) = 100 * (x[:, end] .- x[:, 1])
adjust_exps_trends = Dict()
adjust_exps_trends = Dict(key => compute_depth_diff(adjust_exps[key]) for key in keys(adjust_exps))
exps = ["Initial", "Kappa", "Forcing"]
[adjust_exps_trends[key] = adjust_exps_trends[key] .- adjust_exps_trends["CTRL"] for key in exps] #remove the effect of the baseline
adjust_exps_trends["SUM"] = adjust_exps_trends["CTRL"] .+ adjust_exps_trends["Initial"] .+ 
                            adjust_exps_trends["Kappa"] .+ adjust_exps_trends["Forcing"]

fig, ax = plt.subplots(figsize = (10, 10))
ax.plot(adjust_exps_trends["CTRL"], z, label = "CTRL (Iteration 0)")
ax.plot(adjust_exps_trends["Initial"], z, label = "INIT Effect")
ax.plot(adjust_exps_trends["Kappa"], z, label = "MIXING Effect")
ax.plot(adjust_exps_trends["Forcing"], z, label = "FORCING Effect")
ax.plot(adjust_exps_trends["SUM"], z, label = "SUM", color = "k", linestyle = "--", alpha = 0.6)
ax.plot(adjust_exps_trends["FULL"], z, label = "FULL (Iteration 129)", color = "k", linewidth = 3)
rect = patches.Rectangle((-3, 2000), 6, 1000, linewidth=3, edgecolor="none",facecolor="black", alpha = 0.05)
ax.add_patch(rect)
ax.axvline(0, color = "black", zorder = 0, alpha = 0.55, linewidth = 1)
ax.legend(frameon = false)
ax.set_title("North Pacific Temperature Change (2017 minus 1992)", fontweight = "bold")
ax.set_xlabel("cK", fontweight = "bold"); ax.set_ylabel("Depth [m]", fontweight = "bold")
ax.set_xlim(-3, 3)
ax.set_ylim(1500, 5000); ax.invert_yaxis()
fig
fig.savefig(plotsdir("native/sensitivity_exps/North_Pacific_Temp_Diff.png"))




fig, ax = plt.subplots(figsize = (10, 10))
ax.plot((adjust_exps_trends["Initial"] .+ adjust_exps_trends["CTRL"])./10, z, label = "INIT")
ax.plot(adjust_exps_trends["FULL"] ./10, z, label = "INIT")

rect = patches.Rectangle((-0.1, 2000), 0.2, 1000, linewidth=3, edgecolor="none",facecolor="black", alpha = 0.05)
ax.add_patch(rect)
ax.axvline(0, color = "black", zorder = 0, alpha = 0.55, linewidth = 1)
ax.legend(frameon = false)
ax.set_title("North Pacific Temperature Trend Profile", fontweight = "bold")
ax.set_xlabel("cK per decade", fontweight = "bold"); ax.set_ylabel("Depth [m]", fontweight = "bold")
ax.set_xlim(-.1, .1)
ax.set_ylim(1500, 5000); ax.invert_yaxis()
fig
fig.savefig(plotsdir("native/sensitivity_exps/0.North_Pacific_Sens_Exps.png"))
