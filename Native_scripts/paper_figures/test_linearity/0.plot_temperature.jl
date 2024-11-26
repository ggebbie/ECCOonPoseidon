include("../../../src/intro.jl")
include("./src/intro.jl")

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

adjust_exps_NPAC = jldopen(datadir("NPAC" * "_temperature_sens_exps_seasonal_differences.jld2"))["adjust_exps"]
adjust_exps_NPAC_anomaly = Dict()
adjust_exps_NPAC_anomaly = Dict(key => (adjust_exps_NPAC[key] .- adjust_exps_NPAC[key][:, 1]) for key in keys(adjust_exps_NPAC))


fig, ax = plt.subplots(figsize = (10, 10))
ax.plot(mid_depths(adjust_exps_NPAC["iter0_bulkformula"]), label = "Iteration 0", c = exp_colors["iter0_bulkformula"])
ax.plot(mid_depths(adjust_exps_NPAC["only_wind"]), label = "Wind Stress", c = exp_colors["only_wind"])
ax.plot(mid_depths(adjust_exps_NPAC["mean_tau_adjusts"]), label = "Wind Stress (Climatological)", c = "k", alpha = 0.7)
# ax.plot(mid_depths(adjust_exps_NPAC["mean_tau_noadjusts"]), label = "CTRL (CLIMATOLOGICAL)")
ax.legend()
fig


adjust_exps_NPAC = jldopen(datadir("NPAC_temperature_sens_exps_linearity_test.jld2"))["adjust_exps"]
adjust_exps_NPAC_anomaly = Dict()
adjust_exps_NPAC_anomaly = Dict(key => (adjust_exps_NPAC[key] .- adjust_exps_NPAC[key][:, 1]) for key in keys(adjust_exps_NPAC))


fig, ax = plt.subplots(figsize = (12.5, 5))
r1 = 100 * mid_depths(adjust_exps_NPAC_anomaly["only_init"] .- adjust_exps_NPAC_anomaly["iter0_bulkformula"])
r2 = 100 * mid_depths(adjust_exps_NPAC_anomaly["iter129_bulkformula"] .- adjust_exps_NPAC_anomaly["noinitadjust"])
ax.plot(tecco, r1, label = "Iter0Lin", c = exp_colors["only_wind"])
ax.plot(tecco, r2, label = "Iter129Lin", c = exp_colors["only_wind"])
ax.legend(frameon = false)
ax.set_ylabel("[cK]",  weight = "bold")
ax.set_xlabel("time",  weight = "bold")
ax.grid()
ax.set_title("Mid-Depth North Pacific Temperature Anomaly " * L"\langle \theta' \rangle")

ax2 = ax.twinx()
ax2.plot(tecco, (r1 .- r2))
fig.savefig(plotsdir("native/paper_figures/10.clim_theta_differences.png"), dpi = 200, bbox_inches = "tight")

fig, ax = plt.subplots(figsize = (10, 10))
ax.plot(tecco, 100 .* mid_depths(adjust_exps_NPAC_anomaly["only_wind"] .- adjust_exps_NPAC_anomaly["iter0_bulkformula"]), 
        label = "Effect of Full Iteration 129 Winds")
ax.plot(tecco, 100 .* mid_depths(adjust_exps_NPAC_anomaly["mean_tau_adjusts"] .- adjust_exps_NPAC_anomaly["iter0_bulkformula"]), 
        label = "Effect of Climatological Iteration 129 winds")
ax.legend()
ax.set_ylabel("[cK]")
ax.set_xlabel("time")

x1 = 100 .* mid_depths(adjust_exps_NPAC_anomaly["only_init"] .- adjust_exps_NPAC_anomaly["iter0_bulkformula"])
x2 = 100 .* mid_depths(adjust_exps_NPAC_anomaly["iter129_bulkformula"] .- adjust_exps_NPAC_anomaly["noinitadjust"])
maximum(abs.(x1 .- x2))

ax.set_title("Effect of Initial Condition Control Adjustments \n on North Pacific Temperature Anomaly")
fig.savefig(plotsdir("native/paper_figures/100.clim_theta_differences.png"), dpi = 200, bbox_inches = "tight")
