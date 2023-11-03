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
        τcurl = MeshArrays.curl(τx, τy, Γ)
        curlτ[:, tt] .= zonal_average(τcurl, PAC_msk .* area)[:]

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
fig, ax = plt.subplots()
ax.plot(ϕ_avg, αs)
fig

αs[94]
autocor(WEul[27, 94, :], [1])

autocor(WEul[27, 94, :], [1])
Ψ_bounds = 3
levels = -Ψ_bounds:0.5:Ψ_bounds
cmpday = 100 * 86400
z1 = 150; z2 = 2000; ϕtgt = 39
labels = ["ΔWres", "ΔWEul", "ΔWBol"]

fig = plt.figure(figsize = (8, 8*3))
ax1 = fig.add_subplot(16, 1, (1, 1))
ax2 = fig.add_subplot(16, 1, (3, 6), sharex=ax1)
ax3 = fig.add_subplot(16, 1, (8, 11), sharex=ax2)
ax4 = fig.add_subplot(16, 1, (13, 16), sharex=ax3)
axs = [ax2, ax3, ax4]
titles = [L"\Delta^{\mathbf{W}} \overline{W^{res}}}", 
L"\Delta^{\mathbf{W}} \overline{W}}", L"\Delta^{\mathbf{W}} \overline{W^{*}}}"]
ax1.set_title(L"\Delta^{\mathbf{W}} \mathcal{T}}}")
ax1.set_ylabel("kg m" * L" ^{-3}")

fig
ax1.plot(ϕ_avg, mean(τ, dims = 2), c = "k");
ax1.set_ylim(-1e-12, 1e-12)
cms = []
for (i, Wfill) in enumerate([Wres, WEul, WBol])
    W_ = zeros(size(Wfill)...)
    for it in 1:312
        W_[:, :, it] .= cmpday .* Wfill[:, :, it] ./ ϕ_area'
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
cbar = fig.colorbar(cms[1], ax = [ax1,ax2, ax3, ax4], orientation = "horizontal",
fraction = 0.04, label = " [cm day" * L" ^{-1}" * "]", pad = 0.03)
cbar.set_ticks(-3:1.5:3)

fig
fig.savefig(plotsdir("native/paper_figures/ΔW_timemean_comp_wind.png"), bbox_inches = "tight", dpi = 400)

 
fig.subplots_adjust(wspace = 0.1, hspace = 0.1)
fig.savefig(plotsdir("native/paper_figures/ΔW_timeseries_comp_wind.png"), bbox_inches = "tight", dpi = 400)



nt = 312
corrs = Float64[]
alpha = 0.0
for j = 1:10000
    x = zeros(312)
    [x[i] = alpha .* x[i-1] + rand(Normal()) for i = 2:312]
    y = zeros(312)
    [y[i] = alpha .* y[i-1] + rand(Normal()) for i = 2:312]
    x = (x .- mean(x)) ./ std(x)
    y = (y .- mean(y)) ./ std(y)
    push!(corrs, corr(x, y))
end

fig, ax = plt.subplots()
ax.hist(corrs, bins = 20)
fig