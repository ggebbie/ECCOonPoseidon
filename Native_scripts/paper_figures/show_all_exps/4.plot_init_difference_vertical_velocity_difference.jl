include("../../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, DSP, PyCall, 
    RollingFunctions, LinearAlgebra, Statistics, FFTW, LaTeXStrings
import PyPlot as plt 
import NaNMath as nm
include(srcdir("config_exp.jl"))

include(srcdir("MeshArraysPlots.jl"))

tecco = 1992+1/24:1/12:2018
nz = 50

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
λ_wrap = wrap_λ(λ)
ocean_mask = wet_pts(Γ)
@pyimport cmocean.cm as cmo

region = "PAC"
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)

region = "NPAC"
NPAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)

cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 

H = vertical_sum(cell_depths)
for ij in eachindex(H)
    H.f[ij][H.f[ij] .< 3000] .= 0.0
    H.f[ij][H.f[ij] .> 3000] .= 1.0

end 

PAC_H_mask = PAC_msk .* H

include(srcdir("plot_and_dir_config.jl"))

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); wherenan = (!isnan).(ϕ_avg .* 1)

normalize(x, dims) = (x .- mean(x, dims = dims)) ./ std(x, dims = dims)
normalize(x) = (x[:] .- mean(x[:])) ./ std(x[:])

corr(x, y) = cov(normalize(x), normalize(y)) / (var(normalize(x)) * var(normalize(y)))
corr(x, y, dims) = cov(normalize(x, dims), normalize(y, dims); dims = dims) / (var(normalize(x); dims = dims) * var(normalize(y); dims = dims))

vars =  ["only_init",  "iter0_bulkformula"]
ϕ_avg = jldopen(datadir("Ψ_Eul_timeseries_PAC_only_init" *".jld2"))["ϕ_avg"]
Ws = Dict(); ΨEul = Dict(); ΨBol = Dict(); ΨEulBol = Dict()
z_ref = findall( -3300 .<= -z[:].<= -2000)
wθ_resid = Dict(); wθ_eul = Dict(); wθ_bol = Dict()

for (i, expname) in enumerate(vars)
    read_file = jldopen(datadir("Ψ_Eul_timeseries_PAC_" * expname *".jld2"))
    Ψ_exp = read_file["Ψ_exp_timeseries"]
    ΨEul[expname] = 1 .* Ψ_exp

    read_file = jldopen(datadir("Ψ_EulBol_timeseries_PAC_" * expname *".jld2"))
    Ψ_exp = read_file["Ψ_exp_timeseries"]
    ΨBol[expname] = Ψ_exp .- ΨEul[expname]
    ΨEulBol[expname] = 1 .* Ψ_exp
end

ΨEulBol["only_init"] .-= ΨEulBol["iter0_bulkformula"]

ΨEul["only_init"] .-= ΨEul["iter0_bulkformula"]
ΨBol["only_init"] .-= ΨBol["iter0_bulkformula"]

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); 

Wres = -(ΨEulBol["only_init"][:, 1:end-2, :] .- ΨEulBol["only_init"][:, 3:end, :])
WEul = -(ΨEul["only_init"][:, 1:end-2, :] .- ΨEul["only_init"][:, 3:end, :])
WBol = -(ΨBol["only_init"][:, 1:end-2, :] .- ΨBol["only_init"][:, 3:end, :])

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); 
ϕ_area = zonal_sum(PAC_msk .* area); 
ϕ_area=ϕ_area[(!isnan).(ϕ_avg), :][:][2:end-1]
ϕ_avg = ϕ_avg[(!isnan).(ϕ_avg), :][:][2:end-1]


cmpday = 100 * 86400
mpyear = 1 * 86400 * 365

Ψ_bounds = 3
levels = -Ψ_bounds:0.5:Ψ_bounds


z1 = 150; z2 = 2500; ϕtgt = 38
labels = [L"\Delta^{{\kappa}} \overline{W}", 
L"\Delta^{{\kappa}} \overline{W}", L"\Delta^{{\kappa}} \overline{W^*}"]

fig, ax = plt.subplots(1, 3, figsize = (17, 7), sharey = true)
cms = []
for (i, Wfill) in enumerate([Wres, WEul, WBol])
    ax2 = ax[i]
    W_ = zeros(size(Wfill)...)
    for it in 1:312
        W_[:, :, it] .= mpyear .* Wfill[:, :, it] ./ ϕ_area'
    end
    ax2.set_facecolor("black")
    mean_W =  mean(W_, dims = 3)[:,:,1]
    CM = ax2.contourf(ϕ_avg, z, mean_W,levels = levels, cmap=cmo.balance, 
    vmin = -Ψ_bounds, vmax = Ψ_bounds, extend = "both")
    push!(cms, CM)
    ax2.invert_yaxis()
    ax2.scatter(ϕtgt, z2, c = "black", marker = "*", s = 300)
    ax2.set_title(labels[i])
    ax2.set_xticks(-40:20:60)
    ax2.set_xlim(-34, 60)
    lab = string.(abs.(collect(-40:20:60)))
    lab = lab .* ["°S", "°S", "", "°N", "°N", "°N"]
    ax2.set_xticklabels(lab)
    rect = patches.Rectangle((23, 2000), 59, 1000, linewidth=3, edgecolor="black",facecolor="none", alpha = 0.7)
    ax2.add_patch(rect)
end
fig.subplots_adjust(wspace = 0.1)
ax[1].set_ylabel("Depth [m]", fontweight = "bold")

fig_labs = uppercase.(["a", "b", "c", "d", "e", "f"])
for (i, a) in enumerate(ax)
    a.annotate(fig_labs[i], (0.93, 0.02), fontsize = 25, color = "white", 
    xycoords="axes fraction", fontweight = "bold")
end

fig.colorbar(cms[1], ax = ax[:], orientation = "horizontal",
fraction = 0.05, label = " [m year" * L" ^{-1}" * "]", pad = 0.2)
# fig.savefig(plotsdir("native/paper_figures/ΔW_timemean_comp_kappa.png"), bbox_inches = "tight", dpi = 400)
fig

fig, ax = plt.subplots(1, 2, figsize = (10, 6), sharex = true, sharey = true)
iϕ = Base.findmin(abs.(ϕ_avg .- ϕtgt))[2]
iz2 = Base.findmin(abs.(z .- z2))[2]

W = mpyear .* Wres ./ ϕ_area[iϕ]
ax[1].plot(tecco, W[iz2, iϕ, :][:], c = "k", label = L"\Delta^{{\kappa}} W^{res}")

fig
W = mpyear .* WEul ./ ϕ_area[iϕ]
ax[2].plot(tecco, W[iz2, iϕ, :][:], c = "k", alpha = 0.8, label = L"\Delta^{{\kappa}} W", lw = 2)
W = mpyear .* WBol ./ ϕ_area[iϕ]
ax[2].plot(tecco, W[iz2, iϕ, :][:], c = "k", alpha = 0.5, label =  L"\Delta^{{\kappa}} W^*", lw = 2)

ax[2].set_ylim(-13, 13)
[a.set_ylabel(" [m year" * L"\mathbf{ ^{-1}}" * "]", fontweight = "bold") for a in [ax[1]]]
[a.set_xlabel("time", fontweight = "bold") for a in ax]
for a in ax[:]
    a.legend(frameon = false, ncols = 2, fontsize = 15, loc = "lower left")
end

fig_labs = uppercase.(["a", "b", "c", "d", "e", "f"])
for (i, a) in enumerate(ax)
    a.annotate(fig_labs[i], (0.05, 0.90), fontsize = 25, color = "k", 
    xycoords="axes fraction", fontweight = "bold")
end
fig
fig.subplots_adjust(wspace = 0.1)
# fig.savefig(plotsdir("native/paper_figures/ΔW_timeseries_comp_Mix.png"), bbox_inches = "tight", dpi = 400)
fig



fraction_of_variance(x, i) = sum((x[1:i] .- mean(x)).^2) / (var(x) * length(x))

fig,ax=plt.subplots(1, 3, figsize = (17, 6), sharey = true)
cms = []
for (i, W_) in enumerate([Wres, WEul, WBol])
    axs = ax[i]
    W = 1e-6 .* W_

    ΨCorr = zeros(size(W)[1:2])
    for i in 1:size(W, 1), j in 1:size(W, 2)
        x = W[i, j, :][:]
        ΨCorr[i, j] = 100 * fraction_of_variance(x[:], 12*7)
    end
    
    axs.set_facecolor("black")
    levels = 0:10:100
    CM = axs.contourf(ϕ_avg, z, ΨCorr,levels = levels, cmap=cmo.thermal, 
    vmin = 0, vmax = 100)
    push!(cms, CM)
    axs.invert_yaxis()
    axs.set_xticks(-40:20:60)
    axs.set_xlim(-34, 60)
    lab = string.(abs.(collect(-40:20:60)))
    lab = lab .* ["°S", "°S", "", "°N", "°N", "°N"]
    axs.set_xticklabels(lab)
    axs.set_title("corr(" * L"\Delta^{\tau} \mathcal{T}, ~" *  labels[i] * ")")
    rect = patches.Rectangle((23, 2000), 59, 1000, linewidth=3, edgecolor="black",facecolor="none", alpha = 0.7)
    axs.add_patch(rect)

end
fig
ax[1].set_ylabel("Depth [m]", fontweight = "bold")
fig.subplots_adjust(wspace = 0.1)
fig.colorbar(cms[1], ax = ax[:], orientation = "horizontal", fraction  =0.04, label = "correlation")
fig