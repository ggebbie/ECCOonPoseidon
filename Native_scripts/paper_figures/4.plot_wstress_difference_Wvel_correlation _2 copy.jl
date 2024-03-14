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

function regress_matrix(x)

    nt = length(x)

    # would be neat to return whatever type goes in
    E = ones(Float32,nt,2)
    E[:, 1] .= x
    
    # F = (E'*E)\E' # least squares estimator

    return pinv(E)
end

function get_estimate(x, y) 
    F = regress_matrix(x)
    a, b = F * y
   return  a.* x .+ b
end

var_explained(function_get_estimate(tecco, tecco), tecco)

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
τ_mean["only_wind"] .-= τ_mean["iter0_bulkformula"]

W_mean = Dict()
W_mean["only_wind"]   = get_transports(diagpath, "only_wind", γ)
W_mean["iter0_bulkformula"]   = get_transports(diagpath, "iter0_bulkformula", γ)
W_mean["only_wind"] .-= W_mean["iter0_bulkformula"]

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); 
τ = 1 .* τ_mean["only_wind"][(!isnan).(ϕ_avg), :]

labels = [L"\Delta^{\mathbf{W}} W^{res}", 
L"\Delta^{\mathbf{W}} W", L"\Delta^{\mathbf{W}} W^*"]
save_lab = ["Wres", "WEul", "WBol"]

fig,axs =plt.subplots(figsize = (17, 6), sharey = true)

W = W_mean["only_wind"][:, (!isnan).(ϕ_avg), :]
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
CM = axs.contourf(ϕ_avg, z, ΨCorr.^2,levels = levels, cmap=cmo.balance, 
vmin = -1, vmax = 1, extend = "both")
axs.invert_yaxis()
axs.set_xticks(0:20:60)
axs.set_xlim(0, 60)
fig

fig,axs =plt.subplots(figsize = (17, 6), sharey = true)

ΨCorr = zeros(size(W)[1:2])

i = 37; j = 100
data = DataFrame(X=τ[j, :], Y=Float64.(W[i, j, :][:]))
data = data .* 1e13
ols = lm(@formula(Y ~ X), data)
ols
coefs = coef(ols)
coefs[1] /= 1e13
estimate = coefs[1] .+ (coefs[2] .* τ[j, :])
# var_explained(get_estimate(y, x), x)
corr(estimate,  W[i, j, :][:])
corr(τ[j, :],  W[i, j, :][:])

round(r2(ols); digits=5)
var_explained(estimate, W[i, j, :][:])



corr(x, y) 
fig, ax = plt.subplots(3)
ax[1].plot(x)
ax[2].plot(y)
ax[3].plot(round.(predict(ols), digits=5))

fig
ΨCorr = zeros(size(W)[1:2])
for i in 1:size(W, 1), j in 1:size(W, 2)
    y = W[i, j, :][:]
    x = τ[j, :]
    if (!isnan)(sum(y))
        data = DataFrame(X=τ[j, :], Y=Float64.(W[i, j, :][:]))
        data = data .* 1e13
        ols = lm(@formula(Y ~ X), data)
        ΨCorr[i, j] = round(r2(ols); digits=5)
    else
        ΨCorr[i, j] = NaN
    end
end

fig,axs =plt.subplots(figsize = (17, 6), sharey = true)
axs.set_facecolor("black")
CM = axs.pcolormesh(ϕ_avg, z, ΨCorr, vmin = 0, vmax = 1, cmap = cmo.thermal)
axs.invert_yaxis()
fig.colorbar(CM)
fig

push!(cms, CM)

lab = string.(abs.(collect(-40:20:60)))
lab = lab .* ["°S", "°S", "", "°N", "°N", "°N"]
axs.set_xticklabels(lab)
axs.set_title("corr(" * L"\Delta^{\mathbf{W}} \mathcal{T}, ~" *  labels[i] * ")")
rect = patches.Rectangle((23, 2000), 59, 1000, linewidth=3, edgecolor="black",facecolor="none", alpha = 0.7)
axs.add_patch(rect)

ax[1].set_ylabel("Depth [m]", fontweight = "bold")
fig.subplots_adjust(wspace = 0.1)
fig.colorbar(cms[1], ax = ax[:], orientation = "horizontal", fraction  =0.04, label = "correlation")
fig

# fig.savefig(plotsdir("native/paper_figures/ΔW_timemean_corr_wind.png"), bbox_inches = "tight", dpi = 400)


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
    levels = 0.5:0.1:1
    CM = axs.contourf(ϕ_avg[2:end-1], z, ΨCorr.^2 ,levels = levels, cmap=cmo.thermal, 
    vmin = 0.5, vmax = 1, extend = "both")
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


