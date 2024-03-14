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
# @pyimport cmocean.cm as cmo

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

filelist = searchdir(diagpath["iter0_bulkformula"],"state_2d_set1") # first filter for state_3d_set1
datafilelist_τ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
τx_tmp, τy_tmp = extract_ocnTAU(diagpath, "iter0_bulkformula" , datafilelist_τ[1], γ)

curl_mask = true_curl(τx_tmp, τy_tmp, Γ)
[curl_mask.f[i][(!isnan).(curl_mask.f[i])] .= 1.0 for i = 1:5]
[curl_mask.f[i][(isnan).(curl_mask.f[i])] .= 0.0 for i = 1:5]
curl_mask = curl_mask .* PAC_msk

function get_transports(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid)

    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    datafilelist_τ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    nt = length(datafilelist_τ);
    nϕ = length(ϕ_avg)
    curlτ = zeros(Float32, nϕ, nt);

    @time for tt = 1:nt
        println(tt)

        Tname = datafilelist_τ[tt]
        τx, τy = extract_ocnTAU(diagpath, expname , Tname, γ)
        τcurl = true_curl(τx, τy, Γ)
        [τcurl.f[i][(isnan).(τcurl.f[i])] .= 0.0 for i = 1:5]
        curlτ[:, tt] .= zonal_average(τcurl, curl_mask .* area)[:]

    end

    return curlτ
end

τ_mean = Dict()
τ_mean["only_wind"]   = get_transports(diagpath, "only_wind", γ)
τ_mean["iter0_bulkformula"]   = get_transports(diagpath, "iter0_bulkformula", γ)
τ_mean["only_wind"] .-= τ_mean["iter0_bulkformula"]
using StatsBase

vars =  ["only_wind",  "iter0_bulkformula"]
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

ΨEulBol["only_wind"] .-= ΨEulBol["iter0_bulkformula"]

ΨEul["only_wind"] .-= ΨEul["iter0_bulkformula"]
ΨBol["only_wind"] .-= ΨBol["iter0_bulkformula"]

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); 

τ = 1 .* τ_mean["only_wind"][(!isnan).(ϕ_avg), :]
τ = τ[2:end-1, :]
αs = []
for a in 1:length(τ[:, 1])
    aut = autocor(τ[a, :], [1]); push!(αs, aut[1])
end

Wres = -(ΨEulBol["only_wind"][:, 1:end-2, :] .- ΨEulBol["only_wind"][:, 3:end, :])
WEul = -(ΨEul["only_wind"][:, 1:end-2, :] .- ΨEul["only_wind"][:, 3:end, :])
WBol = -(ΨBol["only_wind"][:, 1:end-2, :] .- ΨBol["only_wind"][:, 3:end, :])


ϕ_avg = zonal_average(ϕ, PAC_msk .* area); 
ϕ_area = zonal_sum(PAC_msk .* area); 
ϕ_area = ϕ_area[(!isnan).(ϕ_avg), :][:][2:end-1]
ϕ_avg = ϕ_avg[(!isnan).(ϕ_avg), :][:][2:end-1]

cmpday = 100 * 86400
mpyear = 1 * 86400 * 365

Ψ_bounds = Int(floor(3 * mpyear / cmpday)) - 1
levels = -Ψ_bounds:1.5:Ψ_bounds

z1 = 150; z2 = 2000; ϕtgt = 39
labels = ["ΔWres", "ΔWEul", "ΔWBol"]

fig = plt.figure(figsize = (8, 8*3))
ax1 = fig.add_subplot(16, 1, (1, 2))
ax2 = fig.add_subplot(16, 1, (4, 6), sharex=ax1)
ax3 = fig.add_subplot(16, 1, (8, 10), sharex=ax2)
ax4 = fig.add_subplot(16, 1, (12, 14), sharex=ax3)
axs = [ax2, ax3, ax4]
titles = [L"\Delta^{{\tau}} \overline{W^{res}}}", 
L"\Delta^{\tau} \overline{W}}", L"\Delta^{\tau} \overline{W^{*}}}"]
ax1.set_title(L"\Delta^{\tau} \mathcal{T}}}")
ax1.set_ylabel("N m" * L" ^{-3} \times 10^{-7}")
fig
mult = 1e7
ax1.plot(ϕ_avg,mean(τ .* mult, dims = 2)[:], c = "k");
fig
ax1.set_ylim(-1, 1)
ax1.grid()
cms = []
for (i, Wfill) in enumerate([Wres, WEul, WBol])
    W_ = zeros(size(Wfill)...)
    for it in 1:312
        W_[:, :, it] .= mpyear .* Wfill[:, :, it] ./ ϕ_area'
    end
    ax = axs[i]
    ax.set_title(titles[i])
    ax.set_facecolor("black")

    mean_W =  mean(W_, dims = 3)[:,:,1]

    CM = ax.contourf(ϕ_avg, z, mean_W,levels = levels, cmap=cmo.balance, 
    vmin = -Ψ_bounds, vmax = Ψ_bounds, extend = "both")
    ax.invert_yaxis()

    ax.scatter(ϕtgt, z2, c = "black", marker = "*", s = 300)
    ax.scatter(ϕtgt, z1, c = "green", marker = "*", s = 300)

    ax.set_xticks(-40:20:60)
    ax.set_xlim(-34, 60)
    ax.set_ylabel("Depth [m]", fontweight = "bold")
    lab = string.(abs.(collect(-40:20:60)))
    lab = lab .* ["°S", "°S", "", "°N", "°N", "°N"]
    ax.set_xticklabels(lab)
    push!(cms, CM)
end
fig.subplots_adjust(hspace = 0.01)
fig
cbar = fig.colorbar(cms[1], ax = [ax1,ax2, ax3, ax4], orientation = "horizontal",
fraction = 0.04, label = " [m year" * L" ^{-1}" * "]", pad = 0.035)
# cbar.set_ticks(-3:1.5:3)

ax1.annotate("A", (0.93, 0.05), fontsize = 25, color = "black", 
xycoords="axes fraction", fontweight = "bold")
fig
fig_labs = uppercase.(["b", "c", "d", "e", "f"])
for (i, a) in enumerate([ax2, ax3, ax4])
    a.annotate(fig_labs[i], (0.93, 0.05), fontsize = 28, color = "white", 
    xycoords="axes fraction", fontweight = "bold")
end
fig.savefig(plotsdir("native/paper_figures/ΔW_timemean_comp_wind.png"), bbox_inches = "tight", dpi = 400)

fig, ax = plt.subplots(2, 2, figsize = (10, 6), sharex = true, sharey = true)
iϕ = Base.findmin(abs.(ϕ_avg .- ϕtgt))[2]
iz1 = Base.findmin(abs.(z .- z1))[2]
iz2 = Base.findmin(abs.(z .- z2))[2]

W = mpyear .* Wres ./ ϕ_area[iϕ]
ax[1].plot(tecco, W[iz1, iϕ, :][:], c = "k", label = L"\Delta^{{\tau}} W^{res}")
ax[2].plot(tecco, W[iz2, iϕ, :][:], c = "green", label = L"\Delta^{{\tau}} W^{res}")

W = mpyear .* WEul ./ ϕ_area[iϕ]
ax[3].plot(tecco, W[iz1, iϕ, :][:], c = "k", alpha = 0.8, label = L"\Delta^{{\tau}} W", lw = 2)
ax[4].plot(tecco, W[iz1, iϕ, :][:], c = "green", alpha = 0.8, label = L"\Delta^{{\tau}} W", lw = 2)

W = mpyear .* WBol ./ ϕ_area[iϕ]
ax[3].plot(tecco, W[iz1, iϕ, :][:], c = "k", alpha = 0.5, label =  L"\Delta^{{\tau}} W^*", lw = 2)
ax[4].plot(tecco, W[iz2, iϕ, :][:], c = "green", alpha = 0.8, label = L"\Delta^{{\tau}} W", lw = 2)

[a.set_ylabel(" [m year" * L"\mathbf{ ^{-1}}" * "]", fontweight = "bold") for a in [ax[1], ax[2]]]
[a.set_xlabel("time", fontweight = "bold") for a in ax]
for a in ax[:]
    a.legend(frameon = false, ncols = 2, fontsize = 15, loc = "lower left")
end

fig_labs = uppercase.(["a", "b", "c", "d", "e", "f"])
for (i, a) in enumerate([ax[1], ax[3], ax[2], ax[4]])
    a.annotate(fig_labs[i], (0.05, 0.83), fontsize = 25, color = "k", 
    xycoords="axes fraction", fontweight = "bold")
end
fig
fig.subplots_adjust(wspace = 0.1)
fig.savefig(plotsdir("native/paper_figures/ΔW_timeseries_comp_wind.png"), bbox_inches = "tight", dpi = 400)
fig

