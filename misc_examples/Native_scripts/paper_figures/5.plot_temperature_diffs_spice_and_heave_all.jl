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

eff_exps = ["iter129_bulkformula", "only_wind", "only_init", "only_kappa", "only_buoyancy", "iter0_bulkformula"]

full_exps_trends = Dict()
full_exps_trends = Dict(key => adjust_exps[key] for key in keys(adjust_exps))
[full_exps_trends[key] .= full_exps_trends[key] .- full_exps_trends[key][:, 1] for key in eff_exps] #remove the effect of the baseline

spice_exps = jldopen(datadir(region * "_temp_sigma2_to_z.jld2"))["θ_dict"]
spice_component = Dict()
[spice_component[key] = spice_exps[key] .- spice_exps[key][:, 1]  for key in eff_exps] #remove the effect of the baseline

heave_component = Dict()
[heave_component[key] = full_exps_trends[key] .- spice_component[key] for key in eff_exps] #remove the effect of the baseline

# full_exps_trends["DIFF"] = full_exps_trends["iter129_bulkformula"] .- full_exps_trends["iter0_bulkformula"]
# spice_component["DIFF"] = spice_component["iter129_bulkformula"] .- spice_component["iter0_bulkformula"]
# heave_component["DIFF"] = heave_component["iter129_bulkformula"] .- heave_component["iter0_bulkformula"]
labels = ["Iteration 129", "Wind Forcing", "Initial Condition", "Mixing Parameter", "Buoyancy Forcing", "Iteration 0"]


sns.set_theme(context = "notebook", style = "ticks",
              palette = colors, rc = custom_params);

fig, ax = plt.subplots(3, 6, figsize = (17, 7.5), sharey = true, sharex = true)
for (i, expt) in enumerate(eff_exps)
    spice = spice_component[expt]; 
    heave = heave_component[expt];
    vmax = 3
    levels = -vmax:0.25:vmax
    plot_contourf(ax, c) = ax.contourf(tecco,  z, 100 .* c, vmin = -vmax * 1.0, vmax = vmax * 1.0, levels = levels, cmap = "bwr", extend = "both")
    plot_contour(ax, c) = ax.contour(tecco,  z, 100 .* c, levels = levels,colors="k")  # Negative contours default to dashed.

    plot_contourf(ax[1, i], full_exps_trends[expt]); 
    # CS = plot_contour(ax[1, i], full_exps_trends[expt]); 
    # ax[1, i].clabel(CS, fontsize=9, inline=true)

    ax[1, i].set_title(labels[i])
    # 
    cmaps = plot_contourf(ax[2, i], heave); 
    # CS = plot_contour(ax[2, i], heave); 
    # ax[2, i].clabel(CS, fontsize=9, inline=true)

    # ax[2, i].set_title("Heave Contribution")
    plot_contourf(ax[3, i], spice); 
    # CS = plot_contour(ax[3, i], spice); 
    # ax[3, i].clabel(CS, fontsize=9, inline=true)

    # ax[3, i].set_title("Spice Contribution")
end

ECCOtour.sigma2grid()
[a.set_yticks(0:1000:3500) for a in ax[:, 1]]
[a.set_ylim(1001, 3499) for a in ax[:, 1]]
[a.invert_yaxis() for a in ax[:, 1]]
fig.tight_layout()
fig

# ax[1].invert_yaxis()
# ax[1].set_ylabel("Depth [m]", fontweight = "bold")
# [a.set_xlabel("time") for a in ax[:]]
fig.colorbar(cmaps, ax = ax[:], label = "cK", orientation = "horizontal", fraction = 0.06)
# fig.savefig(plotsdir("native/paper_figures/5.HeaveSpice.png"), bbox_inches = "tight")

fig
