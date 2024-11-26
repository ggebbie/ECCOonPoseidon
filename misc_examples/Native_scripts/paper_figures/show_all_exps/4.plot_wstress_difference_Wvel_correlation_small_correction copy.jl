include("../../../src/intro.jl")

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

filelist = searchdir(diagpath["iter0_bulkformula"],"state_2d_set1") # first filter for state_3d_set1
datafilelist_τ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
τx_tmp, τy_tmp = extract_ocnTAU(diagpath, "iter0_bulkformula" , datafilelist_τ[1], γ)

curl_mask = ma_curl(τx_tmp, τy_tmp, Γ)
[curl_mask.f[i][(!isnan).(curl_mask.f[i])] .= 1.0 for i = 1:5]
[curl_mask.f[i][(isnan).(curl_mask.f[i])] .= 0.0 for i = 1:5]
curl_mask = curl_mask .* PAC_msk

face, index, _ = findlatlon(λ, ϕ, -150, 55);
curl_mask[face][index]
function get_winds(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid)

    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    datafilelist_τ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    nt = length(datafilelist_τ);
    nϕ = length(ϕ_avg)
    curlτ = zeros(Float32, nt);
    @time for tt = 1:nt
        println(tt)

        Tname = datafilelist_τ[tt]
        τx, τy = extract_ocnTAU(diagpath, expname , Tname, γ)
        τcurl = ma_curl(τx, τy, Γ)
        [τcurl.f[i][(isnan).(τcurl.f[i])] .= 0.0 for i = 1:5]
        
        curlτ[tt] = τcurl[face][index]

    end

    return curlτ
end

extract_loc(ma, loc) = ma[loc]

function get_transports(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid)

    filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
    datafilelist_τ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    nt = length(datafilelist_τ);
    W_zonal_avg = zeros(Float32, nz, nt);
    @time for tt = 1:nt
        println(tt)
        fnameuvw = datafilelist_τ[tt]
        u, v, w = extract_eulerian_velocities(diagpath, expname, fnameuvw, γ)
        for k = 1:50
            W_zonal_avg[k, tt] = w.f[face, k][index]
        end
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

fig,ax=plt.subplots(1, 2, figsize = (17, 6))
ax[1].plot(W_mean["diff"][39, :])
ax[2].plot(τ_mean["diff"])
fig

α = zeros(50)
for i = 1:50
    α[i] = corr(W_mean["diff"][i, :], τ_mean["diff"])
end

fig, ax = plt.subplots()
ax.scatter(α, -z)
fig

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
    axs.set_title("corr(" * L"\Delta^{\tau} \mathcal{T}, ~" *  labels[i] * ")")
    rect = patches.Rectangle((23, 2000), 59, 1000, linewidth=3, edgecolor="black",facecolor="none", alpha = 0.7)
    axs.add_patch(rect)

end
ax[1].set_ylabel("Depth [m]", fontweight = "bold")
fig.subplots_adjust(wspace = 0.1)
fig.colorbar(cms[1], ax = ax[:], orientation = "horizontal", fraction  =0.04, label = "correlation")
fig

fig_labs = uppercase.(["a", "b", "c", "d", "e", "f"])
for (i, a) in enumerate(ax)
    a.annotate(fig_labs[i], (0.93, 0.03), fontsize = 25, color = "white", 
    xycoords="axes fraction", fontweight = "bold")
end
fig
fig.savefig(plotsdir("native/paper_figures/ΔW_timemean_corr_wind.png"), bbox_inches = "tight", dpi = 400)


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
    axs.set_title("corr(" * L"\Delta^{{\tau}} \mathcal{T}, ~" *  labels[i] * ")")
    rect = patches.Rectangle((23, 2000), 59, 1000, linewidth=3, edgecolor="black",facecolor="none", alpha = 0.7)
    axs.add_patch(rect)

end
ax[1].set_ylabel("Depth [m]", fontweight = "bold")
fig.subplots_adjust(wspace = 0.1)
fig.colorbar(cms[1], ax = ax[:], orientation = "horizontal", fraction  =0.04, label = "correlation")
fig

fig.savefig(plotsdir("native/paper_figures/ΔW_timemean_Rsq_wind.png"), bbox_inches = "tight", dpi = 400)


fig
fig, ax = plt.subplots()
ΨCorr = zeros(size(W)[1:2])
for i in 1:size(W, 1), j in 1:size(W, 2)
    x = W[i, j, :][:]
    y = τ[j, :]
    ΨCorr[i, j] = corr(x, y) 
end


