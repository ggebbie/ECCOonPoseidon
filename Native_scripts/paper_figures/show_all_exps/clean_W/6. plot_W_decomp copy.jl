include("../../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour, PyCall, 
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings, DataFrames, Statistics
import NaNMath as nm
import PyPlot as plt
using DSP
include(srcdir("config_exp.jl"))
@pyimport seaborn as sns;

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

include(srcdir("plot_and_dir_config.jl"))

@pyimport cmocean.cm as cmos

pad_arr(x) = vcat(reverse(x), x, reverse(x))
unpad_arr(x, arr_len) = x[arr_len+1:2*arr_len]
function low_pass(signal)
    nt = length(signal)
    signal_mean = mean(signal)
    ff = digitalfilter(Lowpass(1/(7*12), fs = 1),Butterworth(4))
    filtered_signal = filtfilt(ff, pad_arr(signal .- signal_mean))
    return unpad_arr(filtered_signal, nt) .+ signal_mean
end


tecco = Float32.(collect(1992+1/24:1/12:2018))
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= -z[:].<= uplvl)

region = "NPAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)
ctrl_vol = sum(cell_volumes[:, lvls])
area_sum = sum(area .* PAC_msk)

vars =  ["iter0_bulkformula",
        "only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]
WEul = Dict(); WBol = Dict()
ΨEul = Dict(); dθdz = Dict(); 
wθ_resid = Dict(); wθ_eul = Dict(); wθ_bol = Dict()
z_ref = findall( -3400 .<= -z[:].<= -2000)
z_ref = [z_ref[1], z_ref[end] + 1]
for (i, expname) in enumerate(vars)
    read_file = jldopen(datadir("Ψ_EulBol_timeseries_PAC_" * expname *".jld2"))
    Ψ_EulBol = read_file["Ψ_exp_timeseries"]
    read_file = jldopen(datadir("Ψ_Eul_timeseries_PAC_" * expname *".jld2"))
    Ψ_Eul = read_file["Ψ_exp_timeseries"]
    ΨEul[expname] = 1 .* Ψ_Eul
    ϕ_avg = read_file["ϕ_avg"]
    ϕ_ref = findall( 23 .<= ϕ_avg .<= 75)

    Ψ_Bol = Ψ_EulBol .- Ψ_Eul

    WEul[expname] = -Float32.(Ψ_Eul[z_ref, ϕ_ref[1], :])
    WBol[expname] = -Float32.(Ψ_Bol[z_ref, ϕ_ref[1], :])
end

effect_exp =  ["only_wind", "only_kappa", "only_init"]

lw = 3

vars =  ["only_wind", "only_init", "only_kappa", "iter0_bulkformula"]
fig, ax = plt.subplots(1, 2, figsize = (10, 7))
for (i, expname) in enumerate(vars)
    We = 1e-6 .* WEul[expname][1, :][:]
    Wb = 1e-6 .* WBol[expname][1, :][:]
    Wres = (We .+ Wb)
    ax[1, 1].plot(tecco, Wres, label = L"W_{res}' |_{z = 2000}", 
    linewidth = lw, color = exp_colors[expname], alpha = 0.2)

    We = low_pass(We)
    Wb = low_pass(Wb)
    Wres = low_pass(Wres)

    ax[1, 1].plot(tecco, Wres, label = L"W_{res}' |_{z = 2000}", 
    linewidth = lw, color = exp_colors[expname])

    ax[1, 1].spines["top"].set_visible(false)
    ax[1, 1].spines["bottom"].set_visible(false)
    ax[1, 1].set_xticks([])

    We = 1e-6 .* WEul[expname][2, :][:]
    Wb = 1e-6 .* WBol[expname][2, :][:]
    Wres = We .+ Wb
    ax[2].plot(tecco, Wres, label = L"W_{res}' |_{z = 3000}", 
    linewidth = lw, color = exp_colors[expname], alpha = 0.2)
    We = low_pass(We)
    Wb = low_pass(Wb)
    Wres = low_pass(Wres)
    println(cov(Wres, We) / sqrt(var(Wres) * ( var(We))))

    ax[2].plot(tecco, Wres, label = L"W_{res}' |_{z = 3000}", 
    linewidth = lw, color = exp_colors[expname])
    ax[2].spines["top"].set_visible(false)
    ax[2].spines["bottom"].set_visible(false)

    ax[2].set_xticks([])
end

extra_char = [L"^{\mathcal{U}}", L"^{\mathcal{B}}"]
# for i in 1:2
#     ax[2, i].set_ylabel(L"\Delta \mathcal{W}" * extra_char[i] *" [Sv]"); ax[2, i].set_ylim(-3.0, 3.0)
#     ax[2, i].set_yticks(-2:1.0:2)
#     ax[3, i].set_ylabel(L"\Delta \mathcal{W}" * extra_char[i] * " [Sv]"); ax[3, i].set_ylim(-1.5, 1.5)
#     ax[3, i].set_yticks(-1:0.5:1)
#     ax[1, i].set_ylabel(L"\Delta \mathcal{W}" * extra_char[i] * " [Sv]"); ax[1, i].set_ylim(-3., 3.)
#     ax[1, i].set_yticks(-2:1.0:2)
#     ax[2, i].annotate("Initial Condition\nAdjustment Response", (0.7, 0.3), fontsize = 15, 
#     xycoords="axes fraction", ha="center", fontweight = "normal", color = exp_colors["only_init"])
    
#     ax[3, i].annotate("Mixing Parameter\nAdjustment Response", (0.25, 0.1), fontsize = 15, 
#     xycoords="axes fraction", ha="center", fontweight = "normal", color = exp_colors["only_kappa"])
    
#     ax[1, i].annotate("Wind Stress\nAdjustment Response", (0.25, 0.1), fontsize = 15, 
#     xycoords="axes fraction", ha="center", fontweight = "normal", color = exp_colors["only_wind"])

#     # ax[1, i].set_xticks(collect(1993:4:2017)); ax[1, i].grid(); ax[1, i].set_xticks([])
#     # ax[2, i].set_xticks(collect(1993:4:2017)); ax[2, i].grid(); ax[2, i].set_xticks([])
#     # ax[3, i].set_xticks(collect(1993:4:2017)); ax[3, i].grid(); ax[3, i].set_xticks([])


# end

fig
ax[3, 1].spines["bottom"].set_visible(true)
ax[3, 2].spines["bottom"].set_visible(true)
ax[3, 1].set_xticks(collect(1993:4:2017))
ax[3, 2].set_xticks(collect(1993:4:2017))


fig.subplots_adjust(hspace = 0.0, wspace = 0.27)
fig
fig.savefig(plotsdir("native/paper_figures/ΔW_full_faces.png"), bbox_inches = "tight", dpi = 400)


[WEul[expt] .-= WEul["iter0_bulkformula"] for expt in effect_exp]
[WBol[expt] .-= WBol["iter0_bulkformula"] for expt in effect_exp]

vars =  ["only_wind", "only_init", "only_kappa"]
fig, ax = plt.subplots(3, 2, figsize = (15, 15))
for (i, expname) in enumerate(vars)
    We = 1e-6 .* WEul[expname][1, :][:]
    Wb = 1e-6 .* WBol[expname][1, :][:]
    Wres = (We .+ Wb)
    ax[i, 1].plot(tecco, Wres, label = L"W_{res}' |_{z = 2000}", 
    linewidth = lw, color = exp_colors[expname], alpha = 0.5)

    We = low_pass(We)
    Wb = low_pass(Wb)
    Wres = low_pass(Wres)

    ax[i, 1].plot(tecco, Wres, label = L"W_{res}' |_{z = 2000}", 
    linewidth = lw, color = exp_colors[expname])

    ax[i, 1].spines["top"].set_visible(false)
    ax[i, 1].spines["bottom"].set_visible(false)
    ax[i, 1].set_xticks([])

    We = 1e-6 .* WEul[expname][2, :][:]
    Wb = 1e-6 .* WBol[expname][2, :][:]
    Wres = We .+ Wb
    ax[i, 2].plot(tecco, Wres, label = L"W_{res}' |_{z = 3000}", 
    linewidth = lw, color = exp_colors[expname], alpha = 0.5)
    We = low_pass(We)
    Wb = low_pass(Wb)
    Wres = low_pass(Wres)
    println(cov(Wres, We) / sqrt(var(Wres) * ( var(We))))

    ax[i, 2].plot(tecco, Wres, label = L"W_{res}' |_{z = 3000}", 
    linewidth = lw, color = exp_colors[expname])
    ax[i, 2].spines["top"].set_visible(false)
    ax[i, 2].spines["bottom"].set_visible(false)

    ax[i, 2].set_xticks([])
end

extra_char = [L"^{\mathcal{U}}", L"^{\mathcal{B}}"]
for i in 1:2
    ax[2, i].set_ylabel(L"\Delta \mathcal{W}" * extra_char[i] *" [Sv]"); ax[2, i].set_ylim(-3.0, 3.0)
    ax[2, i].set_yticks(-2:1.0:2)
    ax[3, i].set_ylabel(L"\Delta \mathcal{W}" * extra_char[i] * " [Sv]"); ax[3, i].set_ylim(-1.5, 1.5)
    ax[3, i].set_yticks(-1:0.5:1)
    ax[1, i].set_ylabel(L"\Delta \mathcal{W}" * extra_char[i] * " [Sv]"); ax[1, i].set_ylim(-3., 3.)
    ax[1, i].set_yticks(-2:1.0:2)
    ax[2, i].annotate("Initial Condition\nAdjustment Response", (0.7, 0.3), fontsize = 15, 
    xycoords="axes fraction", ha="center", fontweight = "normal", color = exp_colors["only_init"])
    
    ax[3, i].annotate("Mixing Parameter\nAdjustment Response", (0.25, 0.1), fontsize = 15, 
    xycoords="axes fraction", ha="center", fontweight = "normal", color = exp_colors["only_kappa"])
    
    ax[1, i].annotate("Wind Stress\nAdjustment Response", (0.25, 0.1), fontsize = 15, 
    xycoords="axes fraction", ha="center", fontweight = "normal", color = exp_colors["only_wind"])

    # ax[1, i].set_xticks(collect(1993:4:2017)); ax[1, i].grid(); ax[1, i].set_xticks([])
    # ax[2, i].set_xticks(collect(1993:4:2017)); ax[2, i].grid(); ax[2, i].set_xticks([])
    # ax[3, i].set_xticks(collect(1993:4:2017)); ax[3, i].grid(); ax[3, i].set_xticks([])


end

fig
ax[3, 1].spines["bottom"].set_visible(true)
ax[3, 2].spines["bottom"].set_visible(true)
ax[3, 1].set_xticks(collect(1993:4:2017))
ax[3, 2].set_xticks(collect(1993:4:2017))


fig.subplots_adjust(hspace = 0.0, wspace = 0.27)
fig
fig.savefig(plotsdir("native/paper_figures/ΔW_faces.png"), bbox_inches = "tight", dpi = 400)
