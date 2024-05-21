include("../../../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, DSP, PyCall, 
    RollingFunctions, LinearAlgebra, Statistics, FFTW, LaTeXStrings
import PyPlot as plt 
import NaNMath as nm

include(srcdir("config_exp.jl"))
@pyimport matplotlib.patches as patches

include(srcdir("MeshArraysPlots.jl"))

tecco = 1992+1/24:1/12:2018
nz = 50

pad_arr(x) = vcat(reverse(x), x, reverse(x))
unpad_arr(x, arr_len) = x[arr_len+1:2*arr_len]
function low_pass(signal)
    nt = length(signal)
    signal_mean = mean(signal)
    ff = digitalfilter(Lowpass(1/(7*12 + 1), fs = 1),Butterworth(4))
    filtered_signal = filtfilt(ff, pad_arr(signal .- signal_mean))
    return unpad_arr(filtered_signal, nt) .+ signal_mean
end

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
λ_wrap = wrap_λ(λ)
ocean_mask = wet_pts(Γ)
cmo = pyimport("cmocean.cm")

region = "PAC"
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)



region = "NPAC"
NPAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)

cell_depths = get_cell_thickness(NPAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)
area_mask = area .* NPAC_msk

maskW = 1 .* NPAC_msk
[maskW.f[i][iszero.(maskW.f[i])] .= NaN for i in 1:5]
[maskW.f[i][isfinite.(maskW.f[i])] .= 1 for i in 1:5]

include(srcdir("plot_and_dir_config.jl"))

W_mean = jldopen(datadir("W_NPAC_all_times.jld2"))["W"]
lvls = findall( -3000 .<= -z[:].<= -2000)
nlevel = length(lvls)

cm²pday² = (86400^2) * (100^2)
W_diff = Dict()
years = [12*5, 12*15, 12*25]; nyears = length(years)
for expname in ["only_wind", "only_init", "only_kappa"]
    W_diff[expname] = []
    for tt in 1:312
        filt_W = 1 .* (W_mean[expname][:, (nlevel*(tt-1) + 1):(nlevel*(tt))] .- W_mean["iter0_bulkformula"][:, (nlevel*(tt-1) + 1):(nlevel*(tt))])
        filt_W = filt_W .* filt_W

        diff_W = sum(filt_W .* cell_volumes[:, lvls] ) / sum(cell_volumes[:, lvls])
        push!(W_diff[expname], diff_W)
    end
end

# fig, axs = plt.subplots(3, 1, figsize=(14,10))
# axs[1, 1].plot(W_diff["only_wind"] .* cm²pday², alpha = 0.2)
# axs[1, 1].plot(low_pass(W_diff["only_wind"]) .* cm²pday²)

# axs[2, 1].plot(W_diff["only_init"] .* cm²pday², alpha = 0.2)
# axs[2, 1].plot(low_pass(W_diff["only_init"]) .* cm²pday²)

# axs[3, 1].plot(W_diff["only_kappa"] .* cm²pday², alpha = 0.2)
# axs[3, 1].plot(low_pass(W_diff["only_kappa"]) .* cm²pday²)


vars =  ["only_wind", "only_init", "only_kappa"]
fig, ax = plt.subplots(3, 1, figsize = (10, 15), sharex = true)
lw = 2
for (i, expname) in enumerate(vars)
    Wres = W_diff[expname] .* cm²pday²
    ax[i, 1].plot(tecco, Wres, label = L"W_{res}' |_{z = 2000}", 
    linewidth = lw, color = exp_colors[expname], alpha = 0.5)

    Wres = low_pass(Wres)

    ax[i, 1].plot(tecco, Wres, label = L"W_{res}' |_{z = 2000}", 
    linewidth = lw, color = exp_colors[expname])

    ax[i, 1].spines["top"].set_visible(false)
    # ax[i, 1].spines["bottom"].set_visible(false)
    # ax[i, 1].set_xticks([])
end
fig
for i in [[1]]
 ax[1].set_ylim(-100, 600)
    ax[1].set_yticks(-0:100.0:500)
    ax[2].set_ylabel(L"\Delta E_w [Sv]"); ax[2].set_ylim(-10, 60)    
    ax[2].set_yticks(-0:10:50)
    ax[3].set_ylabel(L"\Delta E_w [Sv]"); ax[3].set_ylim(-5, 20)
    ax[3].set_yticks(-0:5.0:15)
    ax[2].annotate("Initial Condition\nAdjustment Response", (0.7, 0.73), fontsize = 17.5, 
    xycoords="axes fraction", ha="center", fontweight = "normal", color = exp_colors["only_init"])
    
    ax[3].annotate("Mixing Parameter\nAdjustment Response", (0.7, 0.85), fontsize =  17.5, 
    xycoords="axes fraction", ha="center", fontweight = "normal", color = exp_colors["only_kappa"])
    
    ax[1].annotate("Wind Stress\nAdjustment Response", (0.7, 0.72), fontsize =  17.5, 
    xycoords="axes fraction", ha="center", fontweight = "normal", color = exp_colors["only_wind"])

    # ax[1, i].set_xticks(collect(1993:4:2017)); ax[1, i].grid(); ax[1, i].set_xticks([])
    # ax[2, i].set_xticks(collect(1993:4:2017)); ax[2, i].grid(); ax[2, i].set_xticks([])
    # ax[3, i].set_xticks(collect(1993:4:2017)); ax[3, i].grid(); ax[3, i].set_xticks([])


end
[a.grid() for a in ax]
fig.suptitle("Mid-Depth North Pacific\nVertical Kinetic Energy Responses", y = 0.93)
[a.set_ylabel(" [cm²day⁻²]") for a in ax]

fig
ax[1].spines["top"].set_visible(true)
ax[3].spines["bottom"].set_visible(true)
ax[3].set_xticks(collect(1993:4:2017))

fig.subplots_adjust(hspace = 0.0, wspace = 0.27)
fig

fig.savefig(plotsdir("native/paper_figures/ΔW_maps.png"), bbox_inches = "tight", dpi = 400)
