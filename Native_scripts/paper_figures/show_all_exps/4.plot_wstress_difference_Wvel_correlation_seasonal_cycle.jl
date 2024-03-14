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

filelist = searchdir(diagpath["climatological_tau"],"state_2d_set1") # first filter for state_3d_set1
datafilelist_τ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
τx_tmp, τy_tmp = extract_ocnTAU(diagpath, "iter0_bulkformula" , datafilelist_τ[1], γ)

curl_mask = ma_curl(τx_tmp, τy_tmp, Γ)
[curl_mask.f[i][(!isnan).(curl_mask.f[i])] .= 1.0 for i = 1:5]
[curl_mask.f[i][(isnan).(curl_mask.f[i])] .= 0.0 for i = 1:5]
curl_mask = curl_mask .* PAC_msk

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
        τcurl = ma_curl(τx, τy, Γ)
        [τcurl.f[i][(isnan).(τcurl.f[i])] .= 0.0 for i = 1:5]

        curlτ[:, tt] .= zonal_average(τcurl, curl_mask .* area)[:]

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
τ_mean["climatological_tau"]   = get_winds(diagpath, "climatological_tau", γ)
τ_mean["only_wind"]   = get_winds(diagpath, "only_wind", γ)
τ_mean["iter0_bulkformula"]   = get_winds(diagpath, "iter0_bulkformula", γ)
τ_mean["diff_clim"] = τ_mean["climatological_tau"] .- τ_mean["iter0_bulkformula"]
τ_mean["diff_full"] = τ_mean["only_wind"] .- τ_mean["iter0_bulkformula"]

W_mean = Dict()
W_mean["climatological_tau"]   = get_transports(diagpath, "climatological_tau", γ)
W_mean["only_wind"]   = get_transports(diagpath, "only_wind", γ)
W_mean["iter0_bulkformula"]   = get_transports(diagpath, "iter0_bulkformula", γ)
W_mean["diff_clim"] = W_mean["climatological_tau"] .- W_mean["iter0_bulkformula"]
W_mean["diff_full"] = W_mean["only_wind"] .- W_mean["iter0_bulkformula"]

fig, ax = plt.subplots(1, 2)
ax[1].plot(τ_mean["diff_clim"][210, :])
ax[1].plot(τ_mean["diff_full"][210, :])
ax[2].plot(W_mean["diff_clim"][43, 210, :])
ax[2].plot(W_mean["diff_full"][43, 210, :])
fig

fig, ax = plt.subplots(2, 2, sharey = false)
# ax[1].plot(τ_mean["iter0_bulkformula"][210, :])
ax[2, 1].plot(τ_mean["only_wind"][210, :])
ax[2, 1].plot(τ_mean["iter0_bulkformula"][210, :])
ax[2, 1].plot(τ_mean["climatological_tau"][210, :])
fig
# ax[1].plot(τ_mean["iter0_bulkformula"][210, :])
ax[2, 2].plot(W_mean["only_wind"][43, 210, :])
ax[2, 2].plot(W_mean["iter0_bulkformula"][43, 210, :])
ax[2, 2].plot(W_mean["climatological_tau"][43, 210, :])
mean(W_mean["climatological_tau"][43, 210, :])
mean(W_mean["only_wind"][43, 210, :])

fig


fig, ax = plt.subplots(1, 2, sharey = true)
ax[1].plot(W_mean["climatological_tau"][15, 230, :])
ax[1].plot(W_mean["iter0_bulkformula"][15, 230, :])

ax[2].plot(W_mean["only_wind"][15, 230, :])
ax[2].plot(W_mean["iter0_bulkformula"][15, 230, :])
fig


ϕ_avg = zonal_average(ϕ, PAC_msk .* area); wherenan = (!isnan).(ϕ_avg .* 1)
new_W = W_mean["only_wind"][:, (!isnan).(ϕ_avg), :]
ϕ_avg = ϕ_avg[(!isnan).(ϕ_avg)]
ϕ_avg[100]
new_W[:, 100,1]

Ws = Dict(); ΨEul = Dict(); ΨBol = Dict(); ΨEulBol = Dict()
z_ref = findall( -3300 .<= -z[:].<= -2000)
wθ_resid = Dict(); wθ_eul = Dict(); wθ_bol = Dict()
ϕ_avg[210]
ϕ_avg[230]
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
τ_mean["only_wind"] = τ_mean["only_wind"]
τ = 1 .* τ_mean["only_wind"][(!isnan).(ϕ_avg), :]
τ = τ[2:end-1, :]

Wres = -(ΨEulBol["only_wind"][:, 1:end-2, :] .- ΨEulBol["only_wind"][:, 3:end, :])
WEul = -(ΨEul["only_wind"][:, 1:end-2, :] .- ΨEul["only_wind"][:, 3:end, :])
WBol = -(ΨBol["only_wind"][:, 1:end-2, :] .- ΨBol["only_wind"][:, 3:end, :])

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); 
area_avg = zonal_sum(PAC_msk .*area)
area_avg = area_avg[(!isnan).(ϕ_avg), :][:]
ϕ_avg = ϕ_avg[(!isnan).(ϕ_avg), :][:]
labels = [L"\Delta^{\tau} W^{res}", 
L"\Delta^{\tau} W", L"\Delta^{\tau} W^*"]
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


