include("../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools
import PyPlot as plt

include(srcdir("plot_and_dir_config.jl"))
sns.set_theme(context = "paper", style = "ticks",
              palette = colors, rc = custom_params);
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
adjust_exps["Difference"] = adjust_exps["iter129_bulkformula"] .- adjust_exps["iter0_bulkformula"]
lw = 2

fig, ax = plt.subplots(figsize = (6.0, 5))
plot_exps = ["iter0_bulkformula", "iter129_bulkformula", "only_init", "only_kappa", "only_buoyancy", "only_wind", "Difference"]
exp_colors["Difference"] = "k"
plot_labels["Difference"] = "Effect of All Control Adjustments"
ax.set_title("Mid-depth North Pacfic Temperature Anomaly", fontweight = "bold")
for expt in plot_exps 
    label = titlecase(plot_labels[expt])
    println(plot_exps)
    anomaly = mid_depths(adjust_exps[expt]) .- mid_depths(adjust_exps[expt])[1]
    ax.plot(tecco, 100 .* anomaly, label =label, color = exp_colors[expt], linewidth = lw)
end
ax.grid()
ax.legend(frameon = false)
ax.set_xlabel("time", fontweight = "bold"); ax.set_ylabel(" θ [cK]", fontweight = "bold")
fig
fig.savefig(plotsdir("native/paper_figures/0.North_Pacific_Sens_Exps_T_simple.png"), bbox_inches = "tight")

# Δadjust_exps = deepcopy(adjust_exps)
# Dexps = ["iter129_bulkformula"]
# [Δadjust_exps[nexp] .-= adjust_exps["iter0_bulkformula"] for nexp in Dexps]
# Δadjust_exps["SUM"] = Δadjust_exps["only_init"] .+ Δadjust_exps["only_kappa"] .+ 
#                       Δadjust_exps["only_wind"] .+ Δadjust_exps["only_buoyancy"]


# fig, ax = plt.subplots(figsize = (10, 7.5))
# plot_exps = [ "iter129_bulkformula"]
# exp_colors["SUM"] = "grey"; plot_labels["SUM"] = "SUM"
# for expt in plot_exps 
#     println(plot_exps)
#     ax.plot(tecco, 100 .* mid_depths(Δadjust_exps[expt]), label = plot_labels[expt], color = exp_colors[expt], linewidth = lw)
# end
# # ax.axhline(0, c = "k",= "--")
# ax.legend(frameon = false)
# ax.set_xlabel("time", fontweight = "bold"); ax.set_ylabel("Δθ [cK]", fontweight = "bold")
# ax.set_ylim(-0.8, 0.8)
# fig
# fig.savefig(plotsdir("native/paper_figures/0.North_Pacific_Sens_Exps_ΔT_simple.png"), bbox_inches = "tight")
