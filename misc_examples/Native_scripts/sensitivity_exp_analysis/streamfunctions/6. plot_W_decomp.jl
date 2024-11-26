include("../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour, PyCall, 
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings, DataFrames, Statistics
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))
@pyimport seaborn as sns;

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

include(srcdir("plot_and_dir_config.jl"))

@pyimport cmocean.cm as cmos


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
    ϕ_ref = findall( 23 .<= ϕ_avg .<= 40)

    Ψ_Bol = Ψ_EulBol .- Ψ_Eul

    WEul[expname] = -Float32.(Ψ_Eul[z_ref, ϕ_ref[1], :])
    WBol[expname] = -Float32.(Ψ_Bol[z_ref, ϕ_ref[1], :])
end

effect_exp =  ["only_wind", "only_kappa", "only_init"]

[WEul[expt] .-= WEul["iter0_bulkformula"] for expt in effect_exp]
[WBol[expt] .-= WBol["iter0_bulkformula"] for expt in effect_exp]

lw = 3

vars =  ["only_init", "only_kappa", "only_wind"]
fig, ax = plt.subplots(3, 2, figsize = (15, 15), sharex = true, sharey = "row")
for (i, expname) in enumerate(vars)
    We = 1e-6 .* WEul[expname][1, :][:]
    Wb = 1e-6 .* WBol[expname][1, :][:]
    Wres = We .+ Wb

    ax[i, 1].plot(tecco, Wres, label = L"W_{res}' |_{z = 2000}", 
    linewidth = lw, color = exp_colors[expname])
    ax[i, 1].plot(tecco, We, label = L"W' |_{z = 2000}", 
    linewidth = lw, color = 0.7 .* exp_colors[expname])
    ax[i, 1].plot(tecco, Wb, label = L"W'^* |_{z = 2000}", 
    linewidth = lw, color = exp_colors[expname], alpha = 0.5)
    ax[i, 1].legend(frameon = false, loc="upper left", bbox_to_anchor=(1, 1))
    ax[i, 1].axhline(0, linestyle = "--", color = "k")

    We = 1e-6 .* WEul[expname][2, :][:]
    Wb = 1e-6 .* WBol[expname][2, :][:]
    Wres = We .+ Wb

    println(cov(Wres, We) / var(Wres))

    ax[i, 2].plot(tecco, Wres, label = L"W_{res}' |_{z = 3000}", 
    linewidth = lw, color = exp_colors[expname])
    ax[i, 2].plot(tecco, We, label = L"W' |_{z = 3000}", 
    linewidth = lw, color = 0.7 .* exp_colors[expname])
    ax[i, 2].plot(tecco, Wb, label = L"W'^* |_{z = 3000}", 
    linewidth = lw, color = exp_colors[expname], alpha = 0.5)
    lgnd = ax[i, 2].legend(frameon = false, loc="upper left", bbox_to_anchor=(1, 1))
    ax[i, 2].axhline(0, linestyle = "--", color = "k")

end
ax[1, 1].set_ylabel("[Sv]"); ax[1, 1].set_ylim(-2.5, 2.5)
ax[2, 1].set_ylabel("[Sv]"); ax[2, 1].set_ylim(-1.0, 1.0)
ax[3, 1].set_ylabel("[Sv]"); ax[3, 1].set_ylim(-3., 3.)
fig.subplots_adjust(wspace = 0.5)
fig

include(srcdir("Signals.jl"))

vars =  ["only_init", "only_kappa", "only_wind"]
fig, ax = plt.subplots(2, 3, figsize = (15, 8), sharex = true, sharey = true)
for (i, expname) in enumerate(vars)
    We = 100 * WEul[expname][1, :][:]
    Wb = 100 * WBol[expname][1, :][:]
    Wres = We .+ Wb

    ax[1, i].plot(roll_freqs, rolling_spectral_density(Wres, N, T), 
    label = L"w_{res}' |_{z = 2000}", 
    linewidth = lw, color = exp_colors[expname])
    ax[1, i].plot(roll_freqs, rolling_spectral_density(We, N, T), 
    label = L"w' |_{z = 2000}", 
    linewidth = lw, color = 0.7 * exp_colors[expname])
    ax[1, i].plot(roll_freqs, rolling_spectral_density(Wb, N, T), 
    label = L"w^*' |_{z = 2000}", 
    linewidth = lw, color = exp_colors[expname], alpha = 0.5)
    # ax[1, i].legend(frameon = false)

    We = 100 * WEul[expname][2, :][:]
    Wb = 100 * WBol[expname][2, :][:]
    Wres = We .+ Wb

    println(cov(Wres, We) / var(Wres))

    ax[2, i].plot(roll_freqs, rolling_spectral_density(Wres, N, T), 
    label = L"w_{res}' |_{z = 3000}", 
    linewidth = lw, color = exp_colors[expname])
    ax[2, i].plot(roll_freqs, rolling_spectral_density(We, N, T), 
    label = L"w' |_{z = 3000}", 
    linewidth = lw, color = 0.7 * exp_colors[expname])
    ax[2, i].plot(roll_freqs, rolling_spectral_density(Wb, N, T), 
    label = L"w^*' |_{z = 3000}", 
    linewidth = lw, color = exp_colors[expname], alpha = 0.5)
    ax[2, 2].legend(frameon = false)

end
ax[1].set_yscale("log")
ax[1].set_xscale("log")
ax[1].invert_xaxis()


fig