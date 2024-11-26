include("../../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools
using  DSP
import PyPlot as plt

pad_arr(x) = vcat(reverse(x), x, reverse(x))
unpad_arr(x, arr_len) = x[arr_len+1:2*arr_len]
function low_pass(signal)
    nt = length(signal)
    signal_mean = mean(signal)
    ff = digitalfilter(Lowpass(1/(7*12 + 1), fs = 1),Butterworth(4))
    filtered_signal = filtfilt(ff, pad_arr(signal .- signal_mean))
    return unpad_arr(filtered_signal, nt) .+ signal_mean
end

function low_pass_2d(signal)
    tmp = 1 .* signal 
    for k in 1:size(tmp)[1]
        tmp[k, :] .= low_pass(tmp[k, :])
    end
    return tmp
end


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

E,F = trend_matrices(Float32.(tecco))
compute_depth_diffs(x) = 100 * (x[:, end] .- x[:, 1])[:]

eff_exps = ["only_wind", "only_init", "only_kappa", "iter0_bulkformula"]

adjust_exps =  jldopen(datadir(region * "_temperature_sens_exps.jld2"))["adjust_exps"]
full_exps_trends = Dict()
full_exps_trends = Dict(key => adjust_exps[key] for key in keys(adjust_exps))
[full_exps_trends[key] .= full_exps_trends[key] .- full_exps_trends[key][:, 1] for key in eff_exps] #remove the effect of the baseline

spice_exps = jldopen(datadir(region * "_temp_sigma2_to_z_reference_middle.jld2"))["θ_dict"]
spice_component = Dict()
[spice_component[key] = spice_exps[key] .- spice_exps[key][:, 1]  for key in eff_exps] #remove the effect of the baseline

heave_component = Dict()
[heave_component[key] = full_exps_trends[key] .- spice_component[key] for key in eff_exps] #remove the effect of the baseline

eff_exps = ["only_wind", "only_init", "only_kappa"]
# [full_exps_trends[key] = full_exps_trends[key] .- full_exps_trends["iter0_bulkformula"]  for key in eff_exps] #remove the effect of the baseline
# [spice_component[key] = spice_exps[key] .- spice_exps["iter0_bulkformula"]  for key in eff_exps] #remove the effect of the baseline
# [heave_component[key] = heave_component[key] .- heave_component["iter0_bulkformula"] for key in eff_exps] #remove the effect of the baseline

# full_exps_trends["DIFF"] = full_exps_trends["iter129_bulkformula"] .- full_exps_trends["iter0_bulkformula"]
# spice_component["DIFF"] = spice_component["iter129_bulkformula"] .- spice_component["iter0_bulkformula"]
# heave_component["DIFF"] = heave_component["iter129_bulkformula"] .- heave_component["iter0_bulkformula"]


sns.set_theme(context = "poster", style = "ticks",
              palette = colors, rc = custom_params);
eff_exps = ["only_wind", "only_init", "only_kappa", "iter0_bulkformula"]
labels_before_befire = [L"\Delta^{\tau}",L"\Delta^{\mathbf{I}}", L"\Delta^{\kappa}", ""]

labels_before = [L"\Delta^{\tau}",L"\Delta^{\mathbf{I}}", L"\Delta^{\kappa}", ""]
labels = ["", "", "", L"_0"]

fig, ax = plt.subplots(4, 2, figsize = (17.5,24), sharey = "row", sharex = false)
fig.tight_layout()
for (i, expt) in enumerate(eff_exps)
    spice = spice_component[expt]; 
    heave = heave_component[expt];
    vmax = 2.8; levels = -vmax:0.2:vmax

    plot_contourf(ax, c) = ax.contourf(tecco,  z[20:46], 100 .* c[20:46, :], vmin = -vmax, vmax = vmax, levels = levels, cmap = cm.balance, extend = "both")
    plot_contour(ax, c) = ax.contour(tecco,  z[20:45], 100 .* c[20:45, :], levels = levels,colors="k")  # Negative contours default to dashed.
    plot_contour2(ax, c) = ax.contour(tecco,  z[45:46], 100 .* c[45:46, :], levels = levels,colors="k")  # Negative contours default to dashed.

    cmaps = plot_contourf(ax[i, 1],  low_pass_2d(heave)); 
    cs2 = plot_contour(ax[i, 1], low_pass_2d(heave)); 
    labels = ax[i,1].clabel(cs2, cs2.levels, fontsize=25, inline=true, fmt = "%.1f", 
    inline_spacing = 30, rightside_up = true, use_clabeltext = true)
    for label in labels
        text = label.get_text()
        label.set_fontweight("bold")
    end
    cs2 = plot_contour2(ax[i, 1], low_pass_2d(heave)); 

    plot_contourf(ax[i, 2], low_pass_2d(spice)); 
    cs2 = plot_contour(ax[i, 2], low_pass_2d(spice)); 
    labels = ax[i, 2].clabel(cs2, cs2.levels, fontsize=25, inline=true, fmt = "%.1f", 
    inline_spacing = 30, rightside_up = true, use_clabeltext = true)
    for label in labels
        text = label.get_text()
        label.set_fontweight("bold")
    end
    cs2 = plot_contour2(ax[i, 2], low_pass_2d(spice)); 

    ax[i, 1].set_ylabel("Depth [m]", fontweight = "bold")
    ax[i, 1].set_xlabel("Time", fontweight = "bold")
    ax[i, 2].set_xlabel("Time", fontweight = "bold")

end
fig
[a.set_yticks(0:1000:4000) for a in ax[:, 1]]
[a.set_ylim(900, 4100) for a in ax[:, 1]]
[a.invert_yaxis() for a in ax[:, 1]]
fig

ax[1, 1].annotate("Heave Anomaly [cK]", (0.5, 1.1), fontsize =  30, rotation = 0, 
xycoords="axes fraction", ha="center", fontweight = "bold", color = "k")

ax[1, 2].annotate("Spice Anomaly [cK]", (0.5, 1.1), fontsize =  30, rotation = 0, 
xycoords="axes fraction", ha="center", fontweight = "bold", color = "k")

fig
ax[1, 2].annotate("Wind Stress\nAdjustment", (1.1, 0.3), fontsize =  30, rotation = 90, 
xycoords="axes fraction", ha="center", fontweight = "normal", color = exp_colors["only_wind"])

ax[2, 2].annotate("Initial Condition\nAdjustment", (1.1, 0.2), fontsize = 30, rotation = 90, 
xycoords="axes fraction", ha="center", fontweight = "normal", color = exp_colors["only_init"])

ax[3, 2].annotate("Mixing Parameter\nAdjustment", (1.1, 0.1), fontsize =  30, rotation = 90, 
xycoords="axes fraction", ha="center", fontweight = "normal", color = exp_colors["only_kappa"])

ax[4, 2].annotate("Iteration 0", (1.1, 0.4), fontsize =  30, rotation = 90, 
xycoords="axes fraction", ha="center", fontweight = "normal", color = exp_colors["iter0_bulkformula"])
fig

fig.savefig(plotsdir("native/paper_figures/heave_spice_differences.png"), bbox_inches = "tight", dpi = 400)