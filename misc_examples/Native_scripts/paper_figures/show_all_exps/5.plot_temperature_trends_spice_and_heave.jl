include("../../../src/intro.jl")

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

lvls = findall( -3000 .<= -z[10:end-3].<= -2000)
lvls
mid_depths(x) = vec(sum(Float32, x[lvls, :] .* ΔV[lvls], dims = 1) / sum(Float64, ΔV[lvls]))

tecco = 1992+1/24:1/12:2018; nz = 50
get_datafiles(expname, key) = filter(x -> occursin("data",x),searchdir(diagpath[expname],key) )

adjust_exps =  jldopen(datadir(region * "_temperature_sens_exps.jld2"))["adjust_exps"]

E,F = trend_matrices(Float32.(tecco))
compute_depth_diffs(x) = 100 * (x[:, end] .- x[:, 1])[:]

eff_exps =  ["iter0_bulkformula", "iter129_bulkformula", "only_buoyancy", 
"only_init", "only_kappa", "only_wind"]

sns.set_theme(context = "paper", style = "ticks",
              palette = colors, rc = custom_params);

full_exps_trends = Dict()
full_exps_trends = Dict(key => adjust_exps[key][10:end-3, :] for key in keys(adjust_exps))
[full_exps_trends[key] .= full_exps_trends[key] .- full_exps_trends[key][:, 1] for key in eff_exps] #remove the effect of the baseline

spice_exps = jldopen(datadir(region * "_temp_sigma2_to_z_Nov.jld2"))["θ_dict"]
spice_component = Dict()
[spice_component[key] = spice_exps[key] .- spice_exps[key][:, 1]  for key in eff_exps] #remove the effect of the baseline

heave_component = Dict()
[heave_component[key] = full_exps_trends[key] .- spice_component[key] for key in eff_exps] #remove the effect of the baseline

fig, axs = plt.subplots(1, 3, figsize = (8, 5), sharey = true, sharex = true)
axs[1].set_title("North Pacific \n Temperature Trend")
axs[2].set_title("Heave Component")
axs[3].set_title("Spice Component")
cmapss = []
plot_labels_list = ["Iteration 0", "Iteration 129", "Buoyancy Forcing", "Initial Condition", "Mixing Parameters", "Wind Stress"]
lw = 2
for (i, expt) in enumerate(eff_exps)
    # trends_dict[expt] = zeros(3)
    println(expt)
    spice = (F * spice_component[expt]')[2, :]
    
    heave =(F * heave_component[expt]')[2, :]
    total = heave .+ spice
    axs[3].plot(100 .* spice, z[10:end-3] ./ 1000, c = exp_colors[expt], linewidth = lw)
    axs[2].plot(100 .* heave, z[10:end-3] ./ 1000, label = plot_labels_list[i], c = exp_colors[expt], linewidth = lw)
    axs[1].plot(100 .* total, z[10:end-3] ./ 1000, c = exp_colors[expt], linewidth = lw)
    println("Heave: ",     round.(mid_depths(100 .* heave), digits = 3))
    println("Spice: ",      round.(mid_depths(100 .* spice), digits = 3))
    println("Total: ",     round.(mid_depths(100 .* total), digits = 3))

end

axs[1].set_yticks(0:1:5)
axs[1].set_ylim(0.5, 3.5)
axs[1].set_ylabel("Depth [km]", fontweight = "bold")
[a.set_xlabel("°C per century", fontweight = "bold") for a in axs]

axs[1].invert_yaxis()
axs[1].set_xlim(-0.175, +0.175)
fig
fig_labs = uppercase.(["a", "b", "c", "d", "e"])
for (i, a) in enumerate(axs)
    a.annotate(fig_labs[i], (0.05, 0.05), fontsize = 15, 
    xycoords="axes fraction", fontweight = "bold")
end
axs[2].legend(loc = "lower center", frameon = false, bbox_to_anchor = (0.5, -0.27), ncols = 3)
[a.grid() for a in axs]
# fig.tight_layout()
fig.subplots_adjust(wspace = 0.1)

fig

fig.savefig(plotsdir("native/paper_figures/5.HeaveSpice_Trends.png"), 
bbox_inches = "tight", dpi = 400)

fig