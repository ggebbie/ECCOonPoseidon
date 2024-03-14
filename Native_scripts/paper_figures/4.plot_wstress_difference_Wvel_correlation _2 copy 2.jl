include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, DSP, PyCall, 
    RollingFunctions, LinearAlgebra, Statistics, FFTW, LaTeXStrings
using DataFrames, GLM, StatsBase
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
# corr(x, y, dims) = cov(x, , normalize(y, dims); dims = d]ad ims) / (var(normalize(x); dims = dims) * var(normalize(y); dims = dims))
var_explained(x, y) = 100 * (1 - (var(y - x) / var(y)) )

MeshA
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
τ_mean["difference"] = τ_mean["only_wind"] .- τ_mean["iter0_bulkformula"]

W_mean = Dict()
W_mean["only_wind"]   = get_transports(diagpath, "only_wind", γ)
W_mean["iter0_bulkformula"]   = get_transports(diagpath, "iter0_bulkformula", γ)
W_mean["difference"] = W_mean["only_wind"] .- W_mean["iter0_bulkformula"]

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); 
τ = 1 .* τ_mean["difference"][(!isnan).(ϕ_avg), :]

labels = [L"\Delta^{\mathbf{W}} W^{res}", 
L"\Delta^{\mathbf{W}} W", L"\Delta^{\mathbf{W}} W^*"]
save_lab = ["Wres", "WEul", "WBol"]

fig,axs =plt.subplots(figsize = (17, 6), sharey = true)

W = W_mean["difference"][:, (!isnan).(ϕ_avg), :]
ϕ_avg = ϕ_avg[(!isnan).(ϕ_avg), :][:]

area_avg = area_avg[(!isnan).(ϕ_avg), :][:]

ΨCorr = zeros(size(W)[1:2])
for i in 1:size(W, 1), j in 1:size(W, 2)
    x = W[i, j, :][:]
    y = τ[j, :]
    ΨCorr[i, j] = corr(x, y) 
end

fig,axs =plt.subplots(figsize = (17, 6), sharey = true)
axs.set_facecolor("black")
levels = -1:0.2:1
CM = axs.contourf(ϕ_avg, z, ΨCorr,levels = levels, cmap=cmo.balance, 
vmin = -1, vmax = 1, extend = "both")
axs.invert_yaxis()
axs.set_xticks(0:20:60)
axs.set_xlim(0, 60)
fig

fig,axs =plt.subplots(figsize = (17, 6), sharey = true)
i = 37; j = 100

Y = mean(W[:, j, :], dims = 2)

fig,axs =plt.subplots(figsize = (17, 6), sharey = true)
axs.plot(Y[:], -z)
fig

Y = W[i, j, :][:]
X = τ[j, :][:]

fig,axs =plt.subplots(figsize = (17, 6), sharey = true)
axs.scatter(X, Y)
fig

Ni = 5
cp1 = zeros(Ni + 1)
cm1 = zeros(Ni + 1)

normalize(x) = (x .- mean(x)) ./ std(x)
corr_normalize(x, y) = corr(normalize(x), normalize(y))
for i in 0:Ni

    cp1[i+1] = corr_normalize(X[i+1:end],Y[1:(end-i)]);
    cm1[i+1] = corr_normalize(X[1:(end-i)],Y[i+1:end]);
end
corrs = vcat(reverse(cm1[2:end]), cp1)

fig,axs =plt.subplots(figsize = (17, 6), sharey = true)
lags = collect(-Ni:Ni)
axs.scatter(lags, corrs)
fig