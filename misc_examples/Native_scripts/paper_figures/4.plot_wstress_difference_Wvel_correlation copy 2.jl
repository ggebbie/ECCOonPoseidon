include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, DSP, PyCall, 
    RollingFunctions, LinearAlgebra, Statistics, FFTW, LaTeXStrings, GLM, DataFrames
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

τ_mean = Dict()
τ_mean["only_wind"]   = get_winds(diagpath, "only_wind", γ)
τ_mean["iter0_bulkformula"]   = get_winds(diagpath, "iter0_bulkformula", γ)
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

W = WEul

ΨCorr = zeros(size(W)[1:2])
Ψpredict = zeros(size(W)[1:2])
for i in 1:size(W, 1), j in 1:size(W, 2)
    y = W[i, j, :][:]
    x = τ[j, :]
    if (!isnan)(sum(y))
        data = DataFrame(X=τ[j, :], Y=Float64.(W[i, j, :][:]))
        data = data .* 1e13
        ols = lm(@formula(Y ~ X), data)
        ΨCorr[i, j] = round(r2(ols); digits=5)
        intercept, slope = coef(ols)
        # y_pred = (slope .* τ[j, :]) .+ (intercept ./ 1e13)
        y_pred = (slope .* τ[j, :])

        Ψpredict[i, j] = mean(round.(y_pred; digits=5))
    else
        ΨCorr[i, j] = NaN
        Ψpredict[i, j] = NaN

    end
end
data = DataFrame(X=tecco, Y=tecco)
ols = lm(@formula(Y ~ X), data)
coef(ols)

fig,axs =plt.subplots(figsize = (17, 6), sharey = true)
axs.set_facecolor("black")
CM = axs.pcolormesh(ϕ_avg[2:end-1], z, ΨCorr, vmin = 0, vmax = 1, cmap = cmo.thermal)
axs.invert_yaxis()
fig.colorbar(CM)
fig


fig,axs =plt.subplots(1, 2, figsize = (17, 6), sharey = true)
axs[1].set_facecolor("black")
axs[2].set_facecolor("black")
vmax = nm.maximum(mean(W, dims = 3)[:, :, 1])
CM = axs[1].pcolormesh(ϕ_avg[2:end-1], z, mean(W, dims = 3)[:, :, 1], vmin = -vmax, vmax = vmax, cmap = cmo.balance)
CM = axs[2].pcolormesh(ϕ_avg[2:end-1], z, Ψpredict, vmin = -vmax, vmax = vmax, cmap = cmo.balance)
axs[1].invert_yaxis()
axs[2].invert_yaxis()

fig.colorbar(CM, ax = axs)
fig


ψ_diff = mean(W, dims = 3)[:, :, 1]
Ψ_diff_reconstruct = 1 .* Ψpredict

ψ_diff[isnan.(ψ_diff)] .= 0.0
Ψ_diff_reconstruct[isnan.(Ψ_diff_reconstruct)] .= 0.0

ψ_diff = reverse(cumsum(reverse(ψ_diff, dims = 2), dims = 2), dims = 2)
Ψ_diff_reconstruct = reverse(cumsum(reverse(Ψ_diff_reconstruct, dims = 2), dims = 2), dims = 2)

Ψ_bounds = 2
levels = collect(-Ψ_bounds:0.1:Ψ_bounds)

fig,axs =plt.subplots(1, 2, figsize = (17, 6), sharey = true)
axs[1].set_facecolor("black")
axs[2].set_facecolor("black")
vmax = nm.maximum(ψ_diff)

CM = axs[1].contourf(ϕ_avg[2:end-1], z, 1e-6 .* mean(ΨEul["only_wind"][:, 2:end-1, :], dims = 3)[:, :, 1], cmap=cmo.delta,levels = levels, 
vmin = -Ψ_bounds, vmax = Ψ_bounds, extend = "both")
CM = axs[2].contourf(ϕ_avg[2:end-1], z, -1e-6 .* Ψ_diff_reconstruct, cmap=cmo.delta,levels = levels, 
vmin = -Ψ_bounds, vmax = Ψ_bounds, extend = "both")
# axs[1].invert_yaxis()
axs[2].invert_yaxis()

fig.colorbar(CM, ax = axs)
fig

