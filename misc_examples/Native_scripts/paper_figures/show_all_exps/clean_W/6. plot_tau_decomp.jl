include("../../../../src/intro.jl")

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

τ_mean = Dict()
τ_mean["climatological_tau"]   = get_winds(diagpath, "climatological_tau", γ)
τ_mean["only_wind"]   = get_winds(diagpath, "only_wind", γ)
τ_mean["iter0_bulkformula"]   = get_winds(diagpath, "iter0_bulkformula", γ)
τ_mean["diff_clim"] = τ_mean["climatological_tau"] .- τ_mean["iter0_bulkformula"]
τ_mean["diff_full"] = τ_mean["only_wind"] .- τ_mean["iter0_bulkformula"]

pad_arr(x) = vcat(reverse(x), x, reverse(x))
unpad_arr(x, arr_len) = x[arr_len+1:2*arr_len]
function low_pass(signal)
    nt = length(signal)
    signal_mean = mean(signal)
    ff = digitalfilter(Lowpass(1/(7*12 + 1), fs = 1),Butterworth(4))
    filtered_signal = filtfilt(ff, pad_arr(signal .- signal_mean))
    return unpad_arr(filtered_signal, nt) .+ signal_mean
end


area_longitude = zonal_sum(NPAC_msk .* area)
τ_mean["diff_full"][isnan.(τ_mean["diff_full"])] .= 0.0
τ = sum(τ_mean["diff_full"] .* area_longitude, dims = 1) ./ sum(area_longitude)
τ = τ[1, :]

lw = 3

fig, ax = plt.subplots()
ax.plot(tecco, τ, linewidth = lw, color = exp_colors["only_wind"], alpha = 0.5)
ax.plot(tecco, low_pass(τ), linewidth = lw, color = exp_colors["only_wind"], alpha = 1)
ax.spines["top"].set_visible(false)
ax.spines["right"].set_visible(false)
ax.set_xticks(collect(1993:4:2017))
ax.set_ylabel( L"\Delta^{\tau} \mathcal{T}" * " [kg/m³]"); ax.set_ylim(-2e-8, 2e-8)
ax.annotate("Wind Stress\nAdjustment Response", (0.25, 0.8), fontsize = 15, 
xycoords="axes fraction", ha="center", fontweight = "normal", color = exp_colors["only_wind"])
ax.set_yticks(collect(-1.5:0.5:1.5) .* 1e-8)
fig
fig.savefig(plotsdir("native/paper_figures/Δτ.png"), bbox_inches = "tight", dpi = 400)

We = 1e-6 .* WEul["only_wind"][2, :][:]
Wb = 1e-6 .* WBol["only_wind"][2, :][:]

Wres = We .+ Wb

fig, ax = plt.subplots()
ax.plot(low_pass(Wres))
ax.twinx().plot(low_pass(τ))
fig