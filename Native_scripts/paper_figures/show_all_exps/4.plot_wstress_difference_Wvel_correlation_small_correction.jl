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
    curlτ = Any[]
    @time for tt = 1:312
        println(tt)

        Tname = datafilelist_τ[tt]
        τx, τy = extract_ocnTAU(diagpath, expname , Tname, γ)
        τcurl = ma_curl(τx, τy, Γ)
        [τcurl.f[i][(isnan).(τcurl.f[i])] .= 0.0 for i = 1:5]
        
        push!(curlτ, τcurl)


    end

    return curlτ
end

extract_loc(ma, loc) = ma[loc]

function get_transports(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid)

    filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
    datafilelist_τ  = filter(x -> occursin("data",x),filelist) # second filter for "data"


    nt = length(datafilelist_τ);
    W_zonal_avg = Any[]
    @time for tt = 1:312
        println(tt)
        fnameuvw = datafilelist_τ[tt]
        u, v, w = extract_eulerian_velocities(diagpath, expname, fnameuvw, γ)
        push!(W_zonal_avg, w)

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


mean_W_diff = 0.0 .* W_mean["diff"][1]; fill!(mean_W_diff, 0.0)
mean_τ_diff = 0.0 .* τ_mean["diff"][1]; fill!(mean_τ_diff, 0.0)

σ_W_diff = 0.0 .* W_mean["diff"][1]; fill!(σ_W_diff, 0.0)
σ_τ_diff = 0.0 .* τ_mean["diff"][1]; fill!(σ_τ_diff, 0.0)

for tt = 1:312
    mean_W_diff .+= (W_mean["diff"][tt] ./ 312)
    mean_τ_diff .+= (τ_mean["diff"][tt] ./ 312)
end

for tt = 1:312
    σ_W_diff .+= (map.(x -> x.^2, W_mean["diff"][tt] .- mean_W_diff) ./ 312)
    σ_τ_diff .+=  (map.(x -> x.^2, τ_mean["diff"][tt] .- mean_τ_diff) ./ 312)
end

covariance = 0.0 .* W_mean["diff"][1]; fill!(covariance, 0.0)
for tt = 1:312
    tmp = map.((x,y) -> x.*y, W_mean["diff"][tt] .- mean_W_diff, τ_mean["diff"][tt] .- mean_τ_diff)
    tmp = tmp ./ 312
    covariance += tmp
end

denom = map.((x,y) -> x.*y, σ_W_diff, σ_τ_diff)
denom = map.((x) -> sqrt.(x), denom)

corrs = map.((x,y) -> x./y, covariance, denom)

fig, ax = plt.subplots(figsize=(14,10), subplot_kw=Dict("projection"=> proj0))
data = 1 .* corrs[:, 43]
nm.maximum((corrs .* PAC_msk)[4, 1])
for ff = 1:5
    cf = ax.pcolormesh(λ[ff], ϕ[ff],  data[ff], transform=projPC, 
    cmap = cmo.balance, rasterized = true, shading = "nearest", vmin = -1, vmax = 1) 
    # push!(CF, cf)
end
fig

corrs_nonans = 1 .* corrs
for ijk in eachindex(corrs_nonans)
    corrs_nonans.f[ijk][isnan.(corrs_nonans.f[ijk])] .= 0.0
end

corrs_zonal = zonal_average(corrs_nonans, cell_volumes)
corrs_zonal[corrs_zonal .== 0.0] .= NaN

fig,axs=plt.subplots( figsize = (17, 6), sharey = true)
cms = []
axs.set_facecolor("black")
levels = -1:0.2:1
ϕ_avg = zonal_average(ϕ, PAC_msk .* area); 
corrs_zonal = corrs_zonal[:, (!isnan).(ϕ_avg)]
ϕ_avg = ϕ_avg[(!isnan).(ϕ_avg), :][:]

CM = axs.contourf(ϕ_avg, z, corrs_zonal,levels = levels, cmap=cmo.balance, 
vmin = -1, vmax = 1, extend = "both")
fig
axs.invert_yaxis()
axs.set_xticks(-40:20:60)
# axs.set_xlim(-34, 60)
fig


lab = string.(abs.(collect(-40:20:60)))
lab = lab .* ["°S", "°S", "", "°N", "°N", "°N"]
axs.set_xticklabels(lab)
# axs.set_title("corr(" * L"\Delta^{\tau} \mathcal{T}, ~" *  labels[i] * ")")
rect = patches.Rectangle((23, 2000), 59, 1000, linewidth=3, edgecolor="black",facecolor="none", alpha = 0.7)
axs.add_patch(rect)
fig


regression_coeffs1 = 0.0 .* W_mean["diff"][1]; fill!(σ_W_diff, 0.0)
regression_coeffs2 = 0.0 .* W_mean["diff"][1]; fill!(σ_W_diff, 0.0)

for tt = 1:312
    tmp1 = map.((x,y) -> x.*y, W_mean["diff"][tt], τ_mean["diff"][tt])
    tmp2 = map.((x,y) -> x.*y, τ_mean["diff"][tt], τ_mean["diff"][tt])
    regression_coeffs1 .+= tmp1
    regression_coeffs2 .+= tmp2
end

α = map.((x,y) -> x./y, regression_coeffs1, regression_coeffs2)

fig, ax = plt.subplots(figsize=(14,10), subplot_kw=Dict("projection"=> proj0))
data = 1 .* α[:, 43]
nm.maximum((data .* PAC_msk)[4, 1])
for ff = 1:5
    cf = ax.pcolormesh(λ[ff], ϕ[ff],  data[ff], transform=projPC, vmin = -100, vmax = 100) 
end
fig


fig, ax = plt.subplots(2, 1, figsize=(14,10), subplot_kw=Dict("projection"=> proj0))
tmp = W_mean["diff"][200]
[tmp.f[ff][iszero.(tmp.f[ff])].=NaN for ff in 1:5]
for ff = 1:5
    cf = ax[1].pcolormesh(λ[ff], ϕ[ff],  α[ff, 43] .* τ_mean["diff"][200][ff], transform=projPC, vmin = -1e-6, vmax = 1e-6, cmap = cmo.balance) 
    cf = ax[2].pcolormesh(λ[ff], ϕ[ff],  tmp[ff, 43], transform=projPC, vmin = -1e-6, vmax = 1e-6, cmap = cmo.balance) 
end
fig

