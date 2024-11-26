include("../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour, PyCall, 
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings, 
    FFTW, RollingFunctions
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))
@pyimport seaborn as sns;

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

include(srcdir("plot_and_dir_config.jl"))

tecco = collect(1992+1/24:1/12:2018)
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= -z[:].<= uplvl)

region = "NPAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)
ctrl_vol = sum(cell_volumes[:, lvls])

region = "PAC"; 

ϕ_avg = jldopen(datadir("Ψ_Eul_timeseries_"*region*"_" * "only_init" *".jld2"))["ϕ_avg"]
ϕ_ref = findall( 23 .<= ϕ_avg .<= 40)
z_ref = findall( -3400 .<= -z[:].<= -2000)

vars =  ["iter0_bulkformula", "iter129_bulkformula","only_init", "only_kappa", "only_wind", "only_buoyancy"]
ΨEul = Dict(); ΨBol = Dict()

for (i, expname) in enumerate(vars)
    read_file = jldopen(datadir("Ψ_EulBol_timeseries_"*region*"_" * expname *".jld2"))
    Ψ_EulBol = read_file["Ψ_exp_timeseries"]

    read_file = jldopen(datadir("Ψ_Eul_timeseries_"*region*"_" * expname *".jld2"))
    Ψ_Eul = read_file["Ψ_exp_timeseries"]

    # Ψs[expname] = Ψ_exp[z_ref, ϕ_ref, :][:]
    Ψ_Bol = Ψ_EulBol .- Ψ_Eul

    WEul = -Float32.(Ψ_Eul[z_ref, ϕ_ref[1], :])
    WBol = -Float32.(Ψ_Bol[z_ref, ϕ_ref[1], :])

    ΨEul[expname] = (WEul[1:end-1, :]  .+ WEul[2:end, :]) ./ 2
    ΨEul[expname] = sum(ΨEul[expname] .* ΔzF[z_ref[1:end-1]], dims = 1)[:] ./ sum(ΔzF[z_ref[1:end-1]])

    ΨBol[expname] = (WBol[1:end-1, :]  .+ WBol[2:end, :]) ./ 2
    ΨBol[expname] = sum(ΨBol[expname] .* ΔzF[z_ref[1:end-1]], dims = 1)[:] ./ sum(ΔzF[z_ref[1:end-1]])
end


vars =  ["iter0_bulkformula", "iter129_bulkformula","only_init", "only_kappa", "only_wind", "only_buoyancy"]
eff_expt =  ["only_wind", "only_kappa", "only_init"]

[ΨEul[expt] .-= ΨEul["iter0_bulkformula"] for expt in eff_expt]
[ΨBol[expt] .-= ΨBol["iter0_bulkformula"] for expt in eff_expt]

fig, ax = plt.subplots(2, 3, figsize = (20, 12), sharex = true, sharey = "row")

for (i, expname) in enumerate(eff_expt)
    ax[1, i].plot(tecco, 1e-6.* ΨEul[expname], c = exp_colors[expname], label = plot_labels_effects[expname], linewidth = 2.5)
    ax[1, i].axhline(0, c = "k", linestyle = "--")
    ax[1, i].set_title(L"\Delta \left<w\right>")

    ax[2, i].plot(tecco, 1e-6.* ΨBol[expname], c = exp_colors[expname], label = plot_labels_effects[expname], linewidth = 2.5)
    ax[2, i].axhline(0, c = "k", linestyle = "--")
    ax[2, i].set_title(L"\Delta \left<w^*\right>")

    # ax.set_title(plot_labels[expname] * " Effect", fontweight = "bold")
end
ax[1, 1].set_ylim(-0.9, 0.9)
ax[1].set_ylabel("Ψ [Sv]")
[a.legend(frameon = false, labelcolor="linecolor", loc = "lower left") for a in ax]
fig
fig.savefig(plotsdir("native/sensitivity_exps/5.ΨEulBol_Effect_timeseries.png"))
fig



nt  = length(ΨEul["only_wind"])
freqs = 1 ./ reverse(FFTW.rfftfreq(nt)[2:end]) 
freqs ./= 12
Δt = 1; N = nt; T = N*Δt

spectral_density(x, N, T) = reverse((2*T * inv(N^2)) .* abs2.(rfft(x .- mean(x))[2:end] ))
rolling_spectral_density(x, N, T) = rollmean(spectral_density(x, N, T), 3)
perc_var_spectral_dens(x, N, T) = rolling_spectral_density(x, N, T) ./ sum(rolling_spectral_density(x, N, T))
roll_freqs = rollmean(freqs, 3)


fig, ax = plt.subplots(2, 3, figsize = (20, 12), sharex = true, sharey = "row")

for (i, expname) in enumerate(eff_expt)
    ax[1, i].plot(roll_freqs, perc_var_spectral_dens(ΨEul[expname], N, T), 
    c = exp_colors[expname], label = plot_labels_effects[expname], linewidth = 2.5)
    ax[1, i].set_title(L"\Delta \left<w\right>" * " Spectra")

    ax[2, i].plot(roll_freqs, perc_var_spectral_dens(ΨBol[expname], N, T), 
    c = exp_colors[expname], label = plot_labels_effects[expname], linewidth = 2.5)
    ax[2, i].set_title(L"\Delta \left<w^*\right>" * " Spectra")

    # ax.set_title(plot_labels[expname] * " Effect", fontweight = "bold")
end
ax[1, 1].invert_xaxis()
[a.set_xscale("log") for a in ax]
[a.set_yscale("log") for a in ax]
[a.set_xlabel("Period (years)") for a in ax[2, :]]
[a.set_ylabel("Spectral Density " * L"[\frac{m^2}{s^2}" * " per year]") for a in ax[:, 1]]
[a.legend(frameon = false, labelcolor="linecolor", loc = "lower left") for a in ax]
fig.savefig(plotsdir("native/sensitivity_exps/5.ΨEulBol_Effect_timeseries_spectra.png"))
fig