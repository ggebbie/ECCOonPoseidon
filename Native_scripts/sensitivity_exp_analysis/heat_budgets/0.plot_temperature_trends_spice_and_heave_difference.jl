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

E,F = trend_matrices(Float32.(tecco))
compute_depth_diffs(x) = 100 * (x[:, end] .- x[:, 1])[:]

eff_exps = ["only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula", "iter0_bulkformula"]

full_exps_trends = Dict()
full_exps_trends = Dict(key => compute_depth_diffs(adjust_exps[key]) for key in keys(adjust_exps))
# [full_exps_trends[key] = full_exps_trends[key] .- full_exps_trends["iter0_bulkformula"] for key in eff_exps] #remove the effect of the baseline

eff_exps = ["iter0_bulkformula", "only_init", "only_kappa", "only_wind", "only_buoyancy","iter129_bulkformula"]

spice_exps = jldopen(datadir(region * "_temp_diffs_sigma2_to_z.jld2"))["diffs"]
spice_exps_trends = Dict()
mid_depth_mean(x) = vec(sum(Float32, x[lvls] .* ΔV[lvls], dims = 1) / sum(Float32, ΔV[lvls]))
[spice_exps_trends[key] = 100 * mid_depth_mean(spice_exps[key]) for key in eff_exps] #remove the effect of the baseline


eff_exps = ["iter0_bulkformula", "only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]

fig, ax = plt.subplots(1, 3, figsize = (17, 7.5), sharey = true)
for (i, expt) in enumerate(eff_exps)
    spice_component = spice_exps_trends[expt]
    heave_component = mid_depth_mean(full_exps_trends[expt]) .- spice_component
    println(expt)
    println(mid_depth_mean(full_exps_trends[expt] .-  full_exps_trends["iter0_bulkformula"]) * 2.6[1])
    println("   ")

    ax[1].bar(i, mid_depth_mean(full_exps_trends[expt]), width=0.95,color = exp_colors[expt])
    ax[2].bar(i, heave_component, width=0.95,color = exp_colors[expt])
    ax[3].bar(i,  spice_component, width=0.95,color = exp_colors[expt])
end
ax[1].set_title("Mid-Depth North Pacific \n Temperature Trend  (z = 2-3km)", fontweight = "bold")
ax[2].set_title("Heave Contribution", fontweight = "bold")
ax[3].set_title("Spice Contribution", fontweight = "bold")

for a in ax
    a.spines["top"].set_visible(false)
    a.spines["right"].set_visible(false)
    a.spines["left"].set_visible(false)
    a.spines["bottom"].set_color("#DDDDDD")
    
    # Second, remove the ticks as well.
    a.tick_params(bottom=false, left=false)
    
    # Third, add a horizontal grid (but keep the vertical grid hidden).
    # Color the lines a light gray as well.
    a.set_axisbelow(true)
    a.yaxis.grid(true, color="#EEEEEE")
    a.xaxis.grid(false)
end
ax[1].set_ylabel("cK", fontweight = "bold")
[a.set_ylim(-1.7, 1.7) for a in ax]

exp_labs = [plot_labels[expt] for expt in eff_exps]
exp_labs[1] = "CTRL"
exp_labs[end] = "FULL"

[a.set_xticks(1:length(exp_labs), exp_labs) for a in ax]
[a.set_xticklabels(exp_labs, rotation = 45, va="center", position=(0,-0.05), fontweight = "bold") for a in ax]
fig
fig.savefig(plotsdir("native/sensitivity_exps/0.North_Pacific_Temp_Diff_HeaveSpiceDecomp.png"), bbox_inches = "tight")
fig
