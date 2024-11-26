include("../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools, Interpolations
using  DSP
import PyPlot as plt
using ColorSchemes
pad_arr(x) = vcat(reverse(x), x, reverse(x))
unpad_arr(x, arr_len) = x[arr_len+1:2*arr_len]
function low_pass(signal)
    nt = length(signal)
    signal_mean = mean(signal)
    ff = digitalfilter(Lowpass(1/25, fs = 1),Butterworth(5))
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
@pyimport matplotlib.colors as mcolors
@pyimport matplotlib.cm as mcm

mScalarMappable  = mcm.ScalarMappable
mLinearSegmentedColormap = mcolors.LinearSegmentedColormap
mNormalize = mcolors.Normalize
mhex2color = mcolors.hex2color


# Step 1: Convert hex codes to RGB tuples
hex_colors = ["#67001f", "#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7",
"#d1e5f0","#92c5de","#4393c3","#2166ac", "#053061"]  # Blue to yellow to red
hex_colors = reverse(hex_colors)
rgb_colors = [mhex2color(color) for color in hex_colors]
# Step 2: Create your colormap
cmap_name = "custom_blue_white_red"
custom_cm = mLinearSegmentedColormap.from_list(cmap_name, rgb_colors)

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
tecco = 1992+1/24:1/12:2018; nz = 50
E,F = trend_matrices(Float32.(tecco))

mid_depths(x) = vec(sum(Float32, x[lvls, :] .* ΔV[lvls], dims = 1) / sum(Float64, ΔV[lvls]))
get_datafiles(expname, key) = filter(x -> occursin("data",x),searchdir(diagpath[expname],key) )
compute_depth_diffs(x) = 100 * (x[:, end] .- x[:, 1])[:]

eff_exps =  ["iter0_bulkformula", "only_wind",
         "iter129_bulkformula", "only_buoyancy", "only_init",
         "only_kappa"]

ctrl_exp = ["only_wind","only_buoyancy", "only_init","only_kappa"]

spice_exps = jldopen(datadir(region * "_temp_sigma2_to_z_wind_rewrite.jld2"))["θ_dict"]
z_dict =  jldopen(datadir(region * "_temp_sigma2_to_z_wind_rewrite.jld2"))["z_dict"]
spice_component = Dict()
[spice_component[key] = spice_exps[key] .- spice_exps[key][:, 1]  for key in eff_exps] #remove the effect of the baseline
spice_component["lin_recon"] = 1. .* spice_component["iter0_bulkformula"]
[spice_component["lin_recon"] .+= (spice_component[key] .- spice_component["iter0_bulkformula"])  for key in ctrl_exp] #remove the effect of the baseline
spice_component["difference"] = spice_component["iter129_bulkformula"] .- spice_component["iter0_bulkformula"]

adjust_exps =  jldopen(datadir(region * "_temperature_sens_exps.jld2"))["adjust_exps"]
adjust_exps2 =  jldopen(datadir(region * "_temperature_sens_seasonal.jld2"))["adjust_exps"]
# adjust_exps["mean_tau_noadjust_redo_bf"] = adjust_exps2["mean_tau_noadjust_redo_bf"]
# adjust_exps["mean_tau_yesadjust_redo_bf"] = adjust_exps2["mean_tau_yesadjust_redo_bf"]

full_exps_trends = Dict()
full_exps_trends = Dict(key => adjust_exps[key] for key in eff_exps)
[full_exps_trends[key] = adjust_exps[key] .- adjust_exps[key][:, 1]  for key in eff_exps] #remove the effect of the baseline
full_exps_trends["lin_recon"] = 1. .* full_exps_trends["iter0_bulkformula"]
[full_exps_trends["lin_recon"] .+= (full_exps_trends[key] .- full_exps_trends["iter0_bulkformula"])  for key in ctrl_exp] #remove the effect of the baseline
full_exps_trends["difference"] = full_exps_trends["iter129_bulkformula"] .- full_exps_trends["iter0_bulkformula"]

heave_component = Dict()
[heave_component[key] = full_exps_trends[key] .- spice_component[key] for key in keys(spice_component)] 

sns.set_style("darkgrid", Dict("axes.facecolor" => ".95", "grid.color" => "0.5"))

eff_exps =  ["iter0_bulkformula","iter129_bulkformula", "difference"]

ncols = 3
fig, ax = plt.subplots(length(eff_exps), ncols, figsize = (15,12), 
sharey = "row", sharex = true)
fig.tight_layout()
include(srcdir("plot_and_dir_config.jl"))
zs = z[:]
for (i, expt) in enumerate(eff_exps)
    spice = spice_component[expt]; 
    heave = heave_component[expt];
    vmax = 3; levels = -vmax:0.25:vmax
    # Step 3: Normalize your data
    norm = mNormalize(vmin=-vmax, vmax=vmax)
    # Step 4: Create a ScalarMappable object
    sm = mScalarMappable(cmap=custom_cm, norm=norm)
    ax[i, 1].contourf(tecco,  zs, 100 .* low_pass_2d(spice .+ heave), 
    vmin = -vmax, vmax = vmax, norm = norm,
    levels = levels, cmap =custom_cm, extend = "both")
    ax[i, 2].contourf(tecco,  zs, 100 .* low_pass_2d(heave), 
    vmin = -vmax, vmax = vmax, norm = norm,
    levels = levels, cmap =custom_cm, extend = "both")
    ax[i, 3].contourf(tecco,  zs, 100 .* low_pass_2d(spice), 
    vmin = -vmax, vmax = vmax, norm = norm,
    levels = levels, cmap =custom_cm, extend = "both")

    for (j, var) in enumerate([heave.+ spice, heave, spice])
        cs2 = ax[i, j].contour(tecco,  zs[1:38], 100 .* low_pass_2d(var)[1:38, :], 
        levels = levels,colors="k", alpha = 0.4, zorder = 10) 
        cs2 = ax[i, j].contour(tecco,  zs[end-5:end], 100 .* low_pass_2d(var)[end-5:end, :], 
        levels = levels,colors="k", alpha = 0.4, zorder = 10) 
        cs2 = ax[i, j].contour(tecco,  zs[38:end-5], 100 .* low_pass_2d(var)[38:end-5, :], 
        levels = levels,colors="k", alpha = 0.4, zorder = 20)  # Negative contours default to dashed.
        
        labels = ax[i,j].clabel(cs2, cs2.levels, fontsize=15.0, inline=true, fmt = "%.2f", 
        inline_spacing = 10, rightside_up = true, use_clabeltext = true)
        # Set the text color and alpha for the labels
        for label in labels
            label.set_color("k")
            label.set_alpha(0.7)
            # label.set_fontweight("bold")
        end
    end
end
[a.set_yticks(0:1000:4000) for a in ax[:, 1]]
[a.set_ylim(1500, 4100) for a in ax[:, 1]]
[a.invert_yaxis() for a in ax[:, 1]]
[a.set_ylabel("Depth [m]") for a in ax[:, 1]]
[a.set_xlabel("time") for a in ax[end, :]]

[a.tick_params(which="both", left = true) for a in ax[:, 1][:]]
[a.tick_params(which="both", bottom=true) for a in ax[end, :][:]]
ax[end, 1].tick_params(which="both", bottom=true, left = true)

ax[1, 1].annotate("North Pacific \n Temperature Anomaly [cK]", (0.5, 1.1), fontsize =  20, rotation = 0, 
xycoords="axes fraction", ha="center", color = "k")
ax[1, 2].annotate("Heave \n Contribution [cK]", (0.5, 1.1), fontsize =  20, rotation = 0, 
xycoords="axes fraction", ha="center", color = "k")
ax[1, 3].annotate("Spice \n Contribution [cK]", (0.5, 1.1), fontsize =  20, rotation = 0, 
xycoords="axes fraction", ha="center", color = "k")
fig.subplots_adjust(hspace = 0.1)

xticks = [1993, 2000, 2008, 2016]
xticks_lab = string.(xticks)
for a in ax[end, :][:]
    a.set_xticks(xticks)
    a.set_xticklabels(xticks_lab, rotation=40)
end

plot_labels_list = ["Iteration 0", "Iteration 129", "Difference"]
exp_colors["mean_tau_noadjust_redo_bf"] = exp_colors["iter0_bulkformula"] 
exp_colors["mean_tau_yesadjust_redo_bf"] = exp_colors["only_wind"] 
exp_colors["lin_recon"] = "k"
exp_colors["difference"] = "k"

alphas = [1, 1, 0.7, 0.7, 1, 1, 1, 1]
for (i, expt) in enumerate(eff_exps)
    ax[i, 3].annotate(plot_labels_list[i], (1.05, 0.5), fontsize =  19, rotation = 270, 
    xycoords="axes fraction", ha="center",va = "center", 
    fontweight = "bold", color = exp_colors[expt], alpha = alphas[i])
end
fig.tight_layout()
fig.subplots_adjust(wspace = 0.05)
# Define colors for the rectangles
# colors = [ exp_colors["only_wind"],  exp_colors["iter0_bulkformula"], "k"]  # Different colors for each row



# Wrap a rectangle around each row using transforms
for (i, color) in enumerate(repeat(["k"], ncols))
    # Get the position and size of the first and last axes in the row

    pos1 = ax[1, i].get_position()
    rect = patches.Rectangle(
        (pos1.x0, pos1.y1 + 0.005),  # Bottom-left corner
        (pos1.x1 - pos1.x0) ,  # Width (span from first to last subplot)
         0.07,  # Height (of the subplot row)
        linewidth=2.5,
        edgecolor="none",
        facecolor="lightgrey",
        transform=fig.transFigure,  # Use figure's coordinate systemm 
        zorder = 0
    )
    # Add the rectangle to the figure
    fig.add_artist(rect)
end

yticks = [2000, 3000, 4000]
yticks_lab = string.(yticks)
for a in ax[:, 1][:]
    a.set_yticks(yticks)
    a.set_yticklabels(yticks_lab, rotation=0)
end
fig

fig.savefig(plotsdir("wind_rewrite/6.heave_spice_differences_iter0_iter129.png"), 
bbox_inches = "tight", dpi = 200)

fig

fig, ax = plt.subplots()
ax.plot(spice_component["iter129_bulkformula"][40, :], c = "b")
ax.plot(spice_component["lin_recon"][40, :], c = "k")
ax.plot(spice_component["iter0_bulkformula"][40, :], c = "r")
fig

fig, ax = plt.subplots()
ax.plot(heave_component["iter129_bulkformula"][40, :], c = "b")
ax.plot(heave_component["lin_recon"][40, :], c = "k")
ax.plot(heave_component["iter0_bulkformula"][40, :], c = "r")
fig