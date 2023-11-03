include("../../src/intro.jl")

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


include(srcdir("plot_and_dir_config.jl"))

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); wherenan = (!isnan).(ϕ_avg .* 1)

normalize(x, dims) = (x .- mean(x, dims = dims)) ./ std(x, dims = dims)
normalize(x) = (x[:] .- mean(x[:])) ./ std(x[:])
square_error(x, y) = sum((x .- mean(x) .* (y .- mean(y))))
corr(x, y) = cov(x, y) / (std(x) * std(y))
# corr(x, y, dims) = cov(x, , normalize(y, dims); dims = dims) / (var(normalize(x); dims = dims) * var(normalize(y); dims = dims))

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

Wres = -(ΨEulBol["only_wind"][:, 1:end-2, :] .- ΨEulBol["only_wind"][:, 3:end, :])
WEul = -(ΨEul["only_wind"][:, 1:end-2, :] .- ΨEul["only_wind"][:, 3:end, :])
WBol = -(ΨBol["only_wind"][:, 1:end-2, :] .- ΨBol["only_wind"][:, 3:end, :])

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); 
area_avg = zonal_sum(PAC_msk .*area)
area_avg = area_avg[(!isnan).(ϕ_avg), :][:]
ϕ_avg = ϕ_avg[(!isnan).(ϕ_avg), :][:]
labels = [L"\Delta^{\mathbf{W}} W^{res}", 
L"\Delta^{\mathbf{W}} W", L"\Delta^{\mathbf{W}} W^*"]
save_lab = ["Wres", "WEul", "WBol"]

fig,ax=plt.subplots(1, 3, figsize = (17, 6), sharey = true)
cms = []
for (i, W_) in enumerate([Wres, WEul, WBol])
    axs = ax[i]
    W = 1e-6 .* W_

    ΨCorr = zeros(size(W)[1:2])
    for i in 1:size(W, 1), j in 1:size(W, 2)
        x = W[i, j, :][:]
        y = τ[j, :]
        ΨCorr[i, j] = corr(x, y) 
    end
    
    axs.set_facecolor("black")
    levels = -1:0.2:1
    CM = axs.contourf(ϕ_avg[2:end-1], z, ΨCorr,levels = levels, cmap=cmo.balance, 
    vmin = -1, vmax = 1, extend = "both")
    push!(cms, CM)
    axs.invert_yaxis()
    axs.set_xticks(-40:20:60)
    axs.set_xlim(-34, 60)
    lab = string.(abs.(collect(-40:20:60)))
    lab = lab .* ["°S", "°S", "", "°N", "°N", "°N"]
    axs.set_xticklabels(lab)
    axs.set_title("corr(" * L"\Delta^{\mathbf{W}} \mathcal{T}, ~" *  labels[i] * ")")
    rect = patches.Rectangle((23, 2000), 59, 1000, linewidth=3, edgecolor="black",facecolor="none", alpha = 0.7)
    axs.add_patch(rect)

end
ax[1].set_ylabel("Depth [m]", fontweight = "bold")
fig.subplots_adjust(wspace = 0.1)
fig.colorbar(cms[1], ax = ax[:], orientation = "horizontal", fraction  =0.04, label = "correlation")
fig
fig.savefig(plotsdir("native/paper_figures/ΔW_timemean_corr_wind.png"), bbox_inches = "tight", dpi = 400)


fig
fig, ax = plt.subplots()
ΨCorr = zeros(size(W)[1:2])
for i in 1:size(W, 1), j in 1:size(W, 2)
    x = W[i, j, :][:]
    y = τ[j, :]
    ΨCorr[i, j] = corr(x, y) 
end


