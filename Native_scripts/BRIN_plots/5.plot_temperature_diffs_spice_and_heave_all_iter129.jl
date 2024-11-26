include("../../src/intro.jl")

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
mid_depths(x) = vec(sum(Float32, x[lvls, :] .* ΔV[lvls], dims = 1) / sum(Float64, ΔV[lvls]))

tecco = 1992+1/24:1/12:2018; nz = 50
get_datafiles(expname, key) = filter(x -> occursin("data",x),searchdir(diagpath[expname],key) )

E,F = trend_matrices(Float32.(tecco))
compute_depth_diffs(x) = 100 * (x[:, end] .- x[:, 1])[:]

eff_exps = ["iter129_bulkformula", "only_wind","iter0_bulkformula"]

adjust_exps =  jldopen(datadir(region * "_temperature_sens_exps_BRIN2.jld2"))["adjust_exps"]
full_exps_trends = Dict()
full_exps_trends = Dict(key => adjust_exps[key] for key in keys(adjust_exps))

spice_exps = jldopen(datadir(region * "_temp_sigma2_to_z_reference_middle_BRIN.jld2"))["θ_dict"]

sns.set_theme(context = "talk", style = "ticks");

fig, ax = plt.subplots(1, 3, figsize = (12,5), sharey = "row", sharex = false)
fig.tight_layout()
titles = ["Iteration 0", "Wind Stress", "Iteration 129"]
for (i, expt) in enumerate(["iter0_bulkformula", "only_wind", "iter129_bulkformula"])
    full = 100 .* mid_depths(1 .* full_exps_trends[expt])
    full = full .- full[1]
    spice = 100 .* mid_depths(1 .* spice_exps[expt]); 
    spice = spice .- spice[1]
    heave = full .- spice

    ax[i].plot(tecco, spice, alpha = 0.2, color = "red")
    ax[i].plot(tecco, heave, alpha = 0.2, color = "blue")
    ax[i].plot(tecco, full, alpha = 0.2, color = "black")

    spice = low_pass(spice)
    heave = low_pass(heave)
    full = low_pass(full)
    ax[i].plot(tecco, spice, alpha = 0.95, color = "red")
    ax[i].plot(tecco, heave, alpha = 0.95, color = "blue")
    ax[i].plot(tecco, full, alpha = 0.95, color = "black")
    ax[i].set_xlabel("time")
    ax[i].set_title(titles[i])
    ax[i].grid()
    ax[i].annotate("Spice", (0.5, 0.87), fontsize =  17, rotation = 0, 
xycoords="axes fraction", ha="center", fontweight = "bold", color = "red")


end
fig
fig.subplots_adjust(wspace = 0.07)

ax[1].annotate("Heave", (0.5, 0.56), fontsize =  17, rotation = 0, 
xycoords="axes fraction", ha="center", fontweight = "bold", color = "blue")
ax[2].annotate("Heave", (0.7, 0.25), fontsize =  17, rotation = 0, 
xycoords="axes fraction", ha="center", fontweight = "bold", color = "blue")
ax[3].annotate("Heave", (0.7, 0.15), fontsize =  17, rotation = 0, 
xycoords="axes fraction", ha="center", fontweight = "bold", color = "blue")

ax[1].annotate("θ'", (0.8, 0.76), fontsize =  19, rotation = 0, 
xycoords="axes fraction", ha="center", fontweight = "bold", color = "black")
fig
ax[2].annotate("θ'", (0.8, 0.5), fontsize =  19, rotation = 0, 
xycoords="axes fraction", ha="center", fontweight = "bold", color = "black")
ax[3].annotate("θ'", (0.8, 0.43), fontsize =  19, rotation = 0, 
xycoords="axes fraction", ha="center", fontweight = "bold", color = "black")
ax[1].set_ylabel("[cK]")

fig

fig.savefig(plotsdir("BRIN/heave_spice_differences_all.png"), bbox_inches = "tight", dpi = 400)