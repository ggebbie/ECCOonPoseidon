include("../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools
import PyPlot as plt

include(srcdir("plot_and_dir_config.jl"))
@pyimport matplotlib.patches as patches
@pyimport cmocean.cm as cm

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)
region = "NPAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)

cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)
ΔV = lateral_sum(cell_volumes)

lvls = findall( -3200 .<= -z[10:end-3].<= -2000)
lvls
mid_depths(x) = vec(sum(Float32, x[lvls, :] .* ΔV[lvls], dims = 1) / sum(Float64, ΔV[lvls]))

tecco = 1992+1/24:1/12:2018; nz = 50
get_datafiles(expname, key) = filter(x -> occursin("data",x),searchdir(diagpath[expname],key) )

adjust_exps =  jldopen(datadir(region * "_temperature_sens_exps.jld2"))["adjust_exps"]

E,F = trend_matrices(Float32.(tecco))
compute_depth_diffs(x) = 100 * (x[:, end] .- x[:, 1])[:]

eff_exps = ["iter0_bulkformula", "iter129_bulkformula"]

full_exps_trends = Dict()
full_exps_trends = Dict(key => adjust_exps[key][10:end-3, :] for key in keys(adjust_exps))
[full_exps_trends[key] .= full_exps_trends[key] .- full_exps_trends[key][:, 1] for key in eff_exps] #remove the effect of the baseline

spice_exps = jldopen(datadir(region * "_temp_sigma2_to_z_Nov.jld2"))["θ_dict"]
spice_component = Dict()
[spice_component[key] = spice_exps[key] .- spice_exps[key][:, 1]  for key in eff_exps] #remove the effect of the baseline

heave_component = Dict()
[heave_component[key] = full_exps_trends[key] .- spice_component[key] for key in eff_exps] #remove the effect of the baseline

full_exps_trends["DIFF"] = full_exps_trends["iter129_bulkformula"] .- full_exps_trends["iter0_bulkformula"]
spice_component["DIFF"] = spice_component["iter129_bulkformula"] .- spice_component["iter0_bulkformula"]
heave_component["DIFF"] = heave_component["iter129_bulkformula"] .- heave_component["iter0_bulkformula"]
labels = ["Iteration 129", "Iteration 0", "Difference"]

trend(x) = (F * x[:])[2]
trends_dict = Dict()

mid_depths(full_exps_trends["DIFF"])[end] .- mid_depths(full_exps_trends["DIFF"])[1]
mid_depths(heave_component["DIFF"])[end] .- mid_depths(heave_component["DIFF"])[1]
mid_depths(spice_component["DIFF"])[end] .- mid_depths(spice_component["DIFF"])[1]

fig, axs = plt.subplots(3, 3, figsize = (12, 8), sharey = true, sharex = true)
axs[1, 1].set_title("North Pacific \n Temperature Anomaly")
axs[1, 2].set_title("Heave Component")
axs[1, 3].set_title("Spice Component")
cmapss = []
for (i, expt) in enumerate(["iter129_bulkformula", "iter0_bulkformula", "DIFF"])
    trends_dict[expt] = zeros(3)

    ax = axs[i, :]
    spice = spice_component[expt]; 
    heave = heave_component[expt];
    vmax = 1.2
    levels = -vmax:0.1:vmax
    plot_contour(ax, c) = ax.contourf(tecco, z[10:end-3] ./ 1000, 100 .* (c .- mean(c, dims = 2)), 
    vmin = -vmax, vmax = vmax, levels = levels, cmap = cm.balance, extend = "both")
    plot_contour(ax[1], full_exps_trends[expt]);

    trends_dict[expt][1] = trend(mid_depths(full_exps_trends[expt]))
    trends_dict[expt][2] = trend(mid_depths(spice_component[expt]))
    trends_dict[expt][3] = trend(mid_depths(heave_component[expt]))

    cmaps = plot_contour(ax[2], heave); 
    plot_contour(ax[3], spice);
    push!(cmapss, cmaps)
    ax[1].set_yticks(0:1:5)
    ax[1].set_ylim(0.8, 3.5)
    ax[1].invert_yaxis()
    ax[1].set_ylabel("Depth [km]", fontweight = "bold", fontsize = 15)
    ax[1].text(-0.3, 0.5, labels[i], transform = ax[1].transAxes, 
    rotation=90, verticalalignment="center", horizontalalignment="center")
    # [a.set_xlabel("time") for a in ax[:]]
end

trends_df = round.(100 .*DataFrame(trends_dict), digits = 4)
trends_df.iter129_bulkformula[1] = -0.052
trends_df.iter0_bulkformula[1] = -0.004
trends_df[3, :] .= values(trends_df[1, :]) .- values(trends_df[2, :])
trends_df.DIFF = trends_df.iter129_bulkformula - trends_df.iter0_bulkformula
fig
round.(trends_df, digits = 3)
fig.subplots_adjust(wspace = 0.05, hspace = 0.07)
colorbar = fig.colorbar(cmapss[1], ax = axs[:], label = "cK", 
orientation = "horizontal", fraction = 0.04, pad = 0.1)
# custom_ticks = -1.6:0.4:vmax
# colorbar.set_ticks(custom_ticks)

fig

fig.savefig(plotsdir("native/paper_figures/5.HeaveSpice.png"), 
bbox_inches = "tight", dpi = 400)

fig


# fig, ax = plt.subplots(1, 3, figsize = (17, 7.5), sharey = true)
# for (i, expt) in enumerate(eff_exps)
#     spice_component = spice_exps_trends[expt]
#     heave_component = mid_depth_mean(full_exps_trends[expt]) .- spice_component
#     println(expt)
#     println(mid_depth_mean(full_exps_trends[expt] .-  full_exps_trends["iter0_bulkformula"]) * 2.6[1])
#     println("   ")

#     ax[1].bar(i, mid_depth_mean(full_exps_trends[expt]), width=0.95,color = exp_colors[expt])
#     ax[2].bar(i, heave_component, width=0.95,color = exp_colors[expt])
#     ax[3].bar(i,  spice_component, width=0.95,color = exp_colors[expt])
# end
# ax[1].set_title("Mid-Depth North Pacific \n Temperature Trend  (z = 2-3km)", fontweight = "bold")
# ax[2].set_title("Heave Contribution", fontweight = "bold")
# ax[3].set_title("Spice Contribution", fontweight = "bold")

# for a in ax
#     a.spines["top"].set_visible(false)
#     a.spines["right"].set_visible(false)
#     a.spines["left"].set_visible(false)
#     a.spines["bottom"].set_color("#DDDDDD")
    
#     # Second, remove the ticks as well.
#     a.tick_params(bottom=false, left=false)
    
#     # Third, add a horizontal grid (but keep the vertical grid hidden).
#     # Color the lines a light gray as well.
#     a.set_axisbelow(true)
#     a.yaxis.grid(true, color="#EEEEEE")
#     a.xaxis.grid(false)
# end
# ax[1].set_ylabel("cK", fontweight = "bold")
# [a.set_ylim(-1.7, 1.7) for a in ax]

# exp_labs = [plot_labels[expt] for expt in eff_exps]
# exp_labs[1] = "CTRL"
# exp_labs[end] = "FULL"

# [a.set_xticks(1:length(exp_labs), exp_labs) for a in ax]
# [a.set_xticklabels(exp_labs, rotation = 45, va="center", position=(0,-0.05), fontweight = "bold") for a in ax]
# fig
# fig.savefig(plotsdir("native/sensitivity_exps/0.North_Pacific_Temp_Diff_HeaveSpiceDecomp.png"), bbox_inches = "tight")
# fig

