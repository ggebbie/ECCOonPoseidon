include("../../src/intro.jl")

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

vars =  ["only_kappa",  "iter0_bulkformula"]
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

ΨEulBol["only_kappa"] .-= ΨEulBol["iter0_bulkformula"]

ΨEul["only_kappa"] .-= ΨEul["iter0_bulkformula"]
ΨBol["only_kappa"] .-= ΨBol["iter0_bulkformula"]

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); 

Wres = -(ΨEulBol["only_kappa"][:, 1:end-2, :] .- ΨEulBol["only_kappa"][:, 3:end, :])
WEul = -(ΨEul["only_kappa"][:, 1:end-2, :] .- ΨEul["only_kappa"][:, 3:end, :])
WBol = -(ΨBol["only_kappa"][:, 1:end-2, :] .- ΨBol["only_kappa"][:, 3:end, :])

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); 
ϕ_area = zonal_sum(PAC_msk .* area); 
ϕ_area=ϕ_area[(!isnan).(ϕ_avg), :][:][2:end-1]
ϕ_avg = ϕ_avg[(!isnan).(ϕ_avg), :][:][2:end-1]

Ψ_bounds = 3
levels = -Ψ_bounds:0.5:Ψ_bounds
cmpday = 100 * 86400

z1 = 150; z2 = 2500; ϕtgt = 38
labels = [L"\Delta^{\mathbf{M}} \overline{W^{res}}", 
L"\Delta^{\mathbf{M}} \overline{W}", L"\Delta^{\mathbf{M}} \overline{W^*}"]


iϕ = Base.findmin(abs.(ϕ_avg .- ϕtgt))[2]

fig, ax = plt.subplots(1, 3, figsize = (17, 5), sharey = true)
cms = []
W_test = zeros(50, 12)

for (i, Wfill) in enumerate([Wres, WEul, WBol])
    ax2 = ax[i]
    W_ = zeros(50, 12)
    for it in 1:12
        mnth= 1 .* Wfill[:, iϕ:end, 1:12:end]
        mnth[isnan.(mnth)] .= 0.0
        W_[:, it] .= mean(cmpday .* mnth ./ ϕ_area[iϕ:end]', dims = [2, 3])
    end
    W_test .= W_
    println(W_)
    # ax2.set_facecolor("black")
    CM = ax2.contourf(1:12, z, W_,levels = levels, cmap=cmo.balance, 
    vmin = -Ψ_bounds, vmax = Ψ_bounds, extend = "both")
    # CM = ax2.contourf(1:12, z, W_)
    push!(cms, CM)
    ax2.invert_yaxis()
    # ax2.scatter(ϕtgt, z2, c = "black", marker = "*", s = 300)
    # ax2.set_title(labels[i])
    # ax2.set_xticks(-40:20:60)
    # ax2.set_xlim(-34, 60)
    # lab = string.(abs.(collect(-40:20:60)))
    # lab = lab .* ["°S", "°S", "", "°N", "°N", "°N"]
    # ax2.set_xticklabels(lab)
    # rect = patches.Rectangle((23, 2000), 59, 1000, linewidth=3, edgecolor="black",facecolor="none", alpha = 0.7)
    # ax2.add_patch(rect)
end
fig
fig.subplots_adjust(wspace = 0.1)
ax[1].set_ylabel("Depth [m]", fontweight = "bold")

fig.colorbar(cms[1], ax = ax[:], orientation = "horizontal",
fraction = 0.05, label = " [cm day" * L" ^{-1}" * "]", pad = 0.2)
fig.savefig(plotsdir("native/paper_figures/ΔW_timemean_comp_kappa.png"), bbox_inches = "tight", dpi = 400)

fig, ax = plt.subplots(1, 2, figsize = (10, 6), sharex = true, sharey = true)
iϕ = Base.findmin(abs.(ϕ_avg .- ϕtgt))[2]
iz2 = Base.findmin(abs.(z .- z2))[2]

W = cmpday .* Wres ./ ϕ_area[iϕ]
ax[1].plot(tecco, W[iz2, iϕ, :][:], c = "k", label = L"\Delta^{\mathbf{M}} W^{res}")
fig
W = cmpday .* WEul ./ ϕ_area[iϕ]
ax[2].plot(tecco, W[iz2, iϕ, :][:], c = "k", alpha = 0.8, label = L"\Delta^{\mathbf{M}} W", lw = 2)
W = cmpday .* WBol ./ ϕ_areaϕ_area[iϕ]
ax[2].plot(tecco, W[iz2, iϕ, :][:], c = "k", alpha = 0.5, label =  L"\Delta^{\mathbf{M}} W^*", lw = 2)

ax[2].set_ylim(-4, 4)
[a.set_ylabel(" [cm day" * L"\mathbf{ ^{-1}}" * "]", fontweight = "bold") for a in [ax[1]]]
[a.set_xlabel("time", fontweight = "bold") for a in ax]
for a in ax[:]
    a.legend(frameon = false, ncols = 2, fontsize = 15, loc = "lower left")
end
fig
fig.subplots_adjust(wspace = 0.1)
# fig.savefig(plotsdir("native/paper_figures/ΔW_timeseries_comp_Mix.png"), bbox_inches = "tight", dpi = 400)
fig