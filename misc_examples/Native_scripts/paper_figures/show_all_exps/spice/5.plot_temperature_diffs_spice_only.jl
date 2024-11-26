include("../../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools
using  DSP
import PyPlot as plt
using GibbsSeaWater

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
sns.set_theme(context = "poster", style = "ticks",
              palette = colors, rc = custom_params);
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

tecco = 1992+1/24:1/12:2018; nz = 50
get_datafiles(expname, key) = filter(x -> occursin("data",x),searchdir(diagpath[expname],key) )

E,F = trend_matrices(Float32.(tecco))
compute_depth_diffs(x) = 100 * (x[:, end] .- x[:, 1])[:]

sig1grid = ECCOtour.sigma2grid()
θ_dict = Dict(); z_dict = Dict(); θσ_dict = Dict()
eff_exps =  ["iter0_bulkformula", "only_wind", "only_buoyancy", "iter129_bulkformula", "only_kappa",  "only_init"]
for expname in eff_exps
    Pσ = jldopen(datadir(expname * region * "_AVG_P_sigma2.jld2"))["P"]
    θσ = jldopen(datadir(expname * region * "_AVG_THETA_sigma2.jld2"))["θ"]

    # avg_Pσ = mean(Pσ, dims = 2); wherenan = isfinite.(avg_Pσ .* 1)
    avg_Pσ = Pσ[:, 1]; wherenan = isfinite.(avg_Pσ .* 1)

    avg_zσ = gsw_z_from_p.(avg_Pσ, 45, 0, 0)

    θ_dict[expname] = θσ
    z_dict[expname] = avg_zσ[:]
end

spice_component = Dict()
[spice_component[key] = θ_dict[key] .- θ_dict[key][:, 1]  for key in eff_exps] #remove the effect of the baseline

plot_labels_list = ["Iteration 0", "Wind Stress", "Buoyancy Forcing", "Iteration 129",
"Mixing Parameters", "Initial Condition"]

vmax = 2.4; levels = -vmax:0.4:vmax

fig,axs=plt.subplots(3, 2, figsize = (17, 20))
for (i, expt) in enumerate(eff_exps)
    axs[i].set_title(plot_labels_list[i])
    spice = 100 .* spice_component[expt][30:61, :]; 
    spice = low_pass_2d(spice)

    z_coord = abs.(z_dict[expt][30:61])

    cmaps = axs[i].contourf(tecco,  z_coord,  spice, vmin = -vmax, vmax = vmax, 
                            levels = levels, cmap = cm.balance, extend = "both")

    cs2 = axs[i].contour(tecco,  z_coord[1:15], spice[1:15, :], levels = levels,colors="k")
    cs2 = axs[i].contour(tecco,  z_coord[15:25], spice[15:25, :], levels = levels,colors="k")
    labels = axs[i].clabel(cs2, cs2.levels, fontsize=20, inline=true, fmt = "%.1f", 
    inline_spacing = 20, rightside_up = true, use_clabeltext = true)

    for label in labels
        text = label.get_text()
        label.set_fontweight("bold")
    end

    cs2 = axs[i].contour(tecco,  z_coord[25:end], spice[25:end, :], levels = levels,colors="k")
    labels = axs[i].clabel(cs2, cs2.levels, fontsize=20, inline=true, fmt = "%.1f", 
    inline_spacing = 20, rightside_up = true, use_clabeltext = true)

    for label in labels
        text = label.get_text()
        label.set_fontweight("bold")
    end
    axs[i].set_ylabel("Depth [m]", fontweight = "bold")
    axs[i].set_xlabel("Time", fontweight = "bold")
end
[a.set_yticks(0:500:3500) for a in axs[:]]
[a.set_ylim(900, 3600) for a in axs[:]]
[a.invert_yaxis() for a in axs[:]]
fig.tight_layout()
fig

fig.colorbar(cmaps, ax = ax[:], label = "cK", orientation = "horizontal", fraction = 0.06)
# fig.savefig(plotsdir("native/paper_figures/5.HeaveSpice.png"), bbox_inches = "tight")

fig

fig, ax = plt.subplots()
ax.plot( spice_component["iter0_bulkformula"][40, :])
ax.plot( low_pass(spice_component["iter0_bulkformula"][40, :]))

fig