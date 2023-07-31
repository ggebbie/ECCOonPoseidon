using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, PyCall, 
    Distributions, FFTW

import PyPlot as plt

@pyimport cmocean.cm as cmo;

θP = load(fnames(section, expname))["θ"]
SP = load(fnames(section, expname))["S"]
P = zero(SP); [P[:, i] = abs.(pstdz) for i = 1:312]
SA = gsw_sa_from_sp.(SP, P, zero(SP), zero(SP))
θ = gsw_ct_from_pt.(SA, θP)


nt = size(θ, 2)

freqs = 1 ./ reverse(rfftfreq(nt)[2:end]) 
freqs ./= 12
Δt = 1; N = nt; T = N*Δt

θ_centered = θ .- mean(θ, dims = 2)

spectral_density(x, N, T) = reverse((2*T * inv(N^2)) .* abs2.(rfft(x)[2:end] ))

perc_frac(x) = x ./ sum(x)

X_spectra  = zeros(length(z), length(freqs))
for k in 1:length(z)
    X_spectra[k, :] .= spectral_density(θ_centered[k, :], N, T)
    X_spectra[k, :] .= X_spectra[k, :] ./ sum(X_spectra[k, :])
    X_spectra[k, :] .*= 100
end

fig, ax = plt.subplots(figsize = (10, 5))
vm = 1.2
levels = -vm:0.2:vm 
ax.set_title("Temperature Anomaly Along Line " * section)

CB = ax.contourf(tecco, abs.(z[20:end]), 100 .* θ_centered[20:end, :], levels = levels, 
vmin = -vm, vmax = vm, cmap = cmo.balance, extend = "both")
fig.colorbar(CB, label = "°C", extend = "both")
ax.invert_yaxis()

ax.set_xlabel("Time"); ax.set_ylabel("Depth [m]")
fig.savefig(plotsdir("native/sections/" * section  * "_Temperature.png"), bbox_inches = "tight")
fig

fig, ax = plt.subplots(figsize = (7.5, 5))
levels = 0:5:80
ax.set_title("Spectra of Temperature Anomalies Along Line " * section)
cm = ax.contourf(freqs, abs.(z[20:end]), X_spectra[20:end, :], vmin = 0, vmax = 80, levels = levels, 
cmap = "Spectral_r", shading = "gouraud")
fig.colorbar(cm, label = "% Spectral Density")
ax.set_xscale("log"); ax.set_xlabel("Period (Years)"); ax.set_ylabel("Depth [m]")
ax.invert_yaxis()
fig.savefig(plotsdir("native/sections/" * section * "_Spectra.png"), bbox_inches = "tight")
fig
