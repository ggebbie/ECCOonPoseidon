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
cmo = pyimport("cmocean.cm")

region = "PAC"
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)

cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)

region = "NPAC"
NPAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)

include(srcdir("plot_and_dir_config.jl"))

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); wherenan = (!isnan).(ϕ_avg .* 1)

normalize(x, dims) = (x .- mean(x, dims = dims)) ./ std(x, dims = dims)
normalize(x) = (x[:] .- mean(x[:])) ./ std(x[:])
square_error(x, y) = sum((x .- mean(x) .* (y .- mean(y))))
corr(x, y) = cov(x, y) / (std(x) * std(y))
# corr(x, y, dims) = cov(x, , normalize(y, dims); dims = dims) / (var(normalize(x); dims = dims) * var(normalize(y); dims = dims))

function get_winds(diagpath::Dict{String, String}, 
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

function get_transports(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid)

    filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
    datafilelist_τ  = filter(x -> occursin("data",x),filelist) # second filter for "data"


    nt = length(datafilelist_τ);
    nϕ = length(ϕ_avg)
    W_zonal_avg = zeros(Float32, nz, nϕ, nt);
    @time for tt = 1:nt
        println(tt)
        fnameuvw = datafilelist_τ[tt]
        u, v, w = extract_eulerian_velocities(diagpath, expname, fnameuvw, γ)
        W_zonal_avg[:, :, tt] .= zonal_average(w, cell_volumes)

    end

    return W_zonal_avg
end

τ_mean = Dict()
τ_mean["only_wind"]   = get_winds(diagpath, "only_wind", γ)
τ_mean["iter0_bulkformula"]   = get_winds(diagpath, "iter0_bulkformula", γ)
τ_mean["diff"] = τ_mean["only_wind"] .- τ_mean["iter0_bulkformula"]

W_mean = Dict()
W_mean["only_wind"]   = get_transports(diagpath, "only_wind", γ)
W_mean["iter0_bulkformula"]   = get_transports(diagpath, "iter0_bulkformula", γ)
W_mean["diff"] = W_mean["only_wind"] .- W_mean["iter0_bulkformula"]

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); 

labels = [L"\Delta^{\mathbf{W}} W^{res}", 
L"\Delta^{\mathbf{W}} W", L"\Delta^{\mathbf{W}} W^*"]
save_lab = ["Wres", "WEul", "WBol"]

fig,axs=plt.subplots(figsize = (13, 6), sharey = true)
cms = []

W = 1.0 .* W_mean["diff"]
τ = 1.0 .* τ_mean["diff"] 

ΨCorr = zeros(size(W)[1:2])
for i in 1:size(W, 1), j in 1:size(W, 2)
    x = W[i, j, :][:]
    y = τ[j, :]
    ΨCorr[i, j] = corr(x, y) 
end
    
axs.set_facecolor("black")
levels = -1:0.1:1

CM = axs.pcolormesh(ϕ_avg[(!isnan).(ϕ_avg)], z, ΨCorr[:, (!isnan).(ϕ_avg)].^2, cmap=cmo.balance, 
vmin = -1, vmax = 1)
fig
push!(cms, CM)
axs.invert_yaxis()
axs.set_xticks(-40:20:60)
axs.set_xlim(-34, 60)
lab = string.(abs.(collect(-40:20:60)))
lab = lab .* ["°S", "°S", "", "°N", "°N", "°N"]
axs.set_xticklabels(lab)
axs.set_title("corr(" * L"\Delta^{\mathbf{W}} \mathcal{T}, ~" *  labels[2] * ")")
rect = patches.Rectangle((23, 2000), 59, 1000, linewidth=3, edgecolor="black",facecolor="none", alpha = 0.7)
axs.add_patch(rect)


axs.set_ylabel("Depth [m]", fontweight = "bold")
fig.subplots_adjust(wspace = 0.1)
fig.colorbar(cms[1], ax = axs, orientation = "horizontal", fraction  =0.04, label = "correlation")
fig

fig,axs=plt.subplots(figsize = (13, 6), sharey = true)
cms = []

W = 1.0 .* W_mean["diff"]
τ = 1.0 .* τ_mean["diff"] 

ΨCoef = zeros(size(W)[1:2])
for i in 1:size(W, 1), j in 1:size(W, 2)
    X = ones(312, 2)
    X[:, 1] .= τ[j, :]
    wherenan_τ = isnan.(τ[j, :])
    y = W[i, j, :][:]
    if isnan(sum(τ[j, :]))
        ΨCoef[i, j] = NaN
    else
        ΨCoef[i, j] = (pinv(X) * y)[1]
    end
end
axs.set_facecolor("black")
levels = -1:0.1:1
nm_max = nm.maximum(ΨCoef)
CM = axs.pcolormesh(ϕ_avg[(!isnan).(ϕ_avg)], z, ΨCoef[:, (!isnan).(ϕ_avg)], cmap=cmo.balance, 
vmin = -1e6, vmax = 1e6)
fig
push!(cms, CM)
axs.invert_yaxis()
axs.set_xticks(-40:20:60)
axs.set_xlim(-34, 60)
lab = string.(abs.(collect(-40:20:60)))
lab = lab .* ["°S", "°S", "", "°N", "°N", "°N"]
axs.set_xticklabels(lab)
axs.set_title("corr(" * L"\Delta^{\mathbf{W}} \mathcal{T}, ~" *  labels[2] * ")")
rect = patches.Rectangle((23, 2000), 59, 1000, linewidth=3, edgecolor="black",facecolor="none", alpha = 0.7)
axs.add_patch(rect)


axs.set_ylabel("Depth [m]", fontweight = "bold")
fig.subplots_adjust(wspace = 0.1)
fig.colorbar(cms[1], ax = axs, orientation = "horizontal", fraction  =0.04, label = "correlation")
fig

# fig.savefig(plotsdir("native/paper_figures/ΔW_timemean_corr_wind.png"), bbox_inches = "tight", dpi = 400)
