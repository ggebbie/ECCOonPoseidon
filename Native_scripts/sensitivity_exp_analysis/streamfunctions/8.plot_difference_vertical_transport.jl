include("../../../src/intro.jl")

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
# @pyimport cmocean.cm as cmo

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

vars =  ["only_init",  "iter0_bulkformula", "only_kappa", "only_wind"]
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
ΨEulBol["only_kappa"] .-= ΨEulBol["iter0_bulkformula"]
ΨEulBol["only_init"] .-= ΨEulBol["iter0_bulkformula"]

ΨEul["only_wind"] .-= ΨEul["iter0_bulkformula"]
ΨEul["only_kappa"] .-= ΨEul["iter0_bulkformula"]
ΨEul["only_init"] .-= ΨEul["iter0_bulkformula"]

ΨBol["only_wind"] .-= ΨBol["iter0_bulkformula"]
ΨBol["only_kappa"] .-= ΨBol["iter0_bulkformula"]
ΨBol["only_init"] .-= ΨBol["iter0_bulkformula"]

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); 
ϕ_avg = ϕ_avg[(!isnan).(ϕ_avg), :][:]

##############

fig, ax = plt.subplots(1, 3, figsize = (18, 5), sharey = true)
labels = ["WIND", "INIT", "MIXING"]
Ψ_bounds = 0.3
levels = -Ψ_bounds:0.005:Ψ_bounds
cf = []
for (i, expt) in enumerate(["only_wind", "only_init", "only_kappa"])
    ax[3].set_ylim(500, 6000)

    ax[i].set_title(labels[i], fontweight = "bold")
    predi = -(ΨEulBol[expt][:, 1:end-2, :] .- ΨEulBol[expt][:, 3:end, :])
    mean_Ψ =  1e-6 .* mean(predi, dims = 3)[:,:,1]
    CM = ax[i].contourf(ϕ_avg[2:end-1], z, mean_Ψ,levels = levels, cmap=cmo.balance, 
    vmin = -Ψ_bounds, vmax = Ψ_bounds, extend = "both")
    push!(cf, CM)
    ax[i].invert_yaxis()
    ax[i].set_xlabel("Latitude", fontweight = "bold")

end
ax[1].set_ylabel("Depth [m]", fontweight = "bold")
fig.subplots_adjust(wspace = 0.1)
levels = -Ψ_bounds:0.1:Ψ_bounds
cbar = fig.colorbar(cf[1], ax= ax[:], orientation = "vertical", fraction = 0.04, label = L"[\mathbf{\Delta W_{res}}]")
cbar.set_ticks(levels)
fig
fig.savefig(plotsdir("native/sensitivity_exps/3.Pacific_Sens_Exps_ΔW_res.png"), bbox_inches = "tight")

fig, ax = plt.subplots(1, 3, figsize = (18, 5), sharey = true)
cf = []
Ψ_bounds = 0.3
levels = -Ψ_bounds:0.005:Ψ_bounds
for (i, expt) in enumerate(["only_wind", "only_init", "only_kappa"])
    ax[3].set_ylim(500, 6000)

    ax[i].set_title(labels[i], fontweight = "bold")
    predi = -(ΨEul[expt][:, 1:end-2, :] .- ΨEul[expt][:, 3:end, :])
    mean_Ψ =  1e-6 .* mean(predi, dims = 3)[:,:,1]
    CM = ax[i].contourf(ϕ_avg[2:end-1], z, mean_Ψ,levels = levels, cmap=cmo.balance, 
    vmin = -Ψ_bounds, vmax = Ψ_bounds, extend = "both")
    ax[i].invert_yaxis()
    push!(cf, CM)

    ax[i].set_xlabel("Latitude", fontweight = "bold")

end
ax[1].set_ylabel("Depth [m]", fontweight = "bold")
fig.subplots_adjust(wspace = 0.1)
levels = -Ψ_bounds:0.1:Ψ_bounds
cbar = fig.colorbar(cf[1], ax= ax[:], orientation = "vertical", fraction = 0.04, label = L"[\mathbf{\Delta W}]")
cbar.set_ticks(levels)
fig
fig.savefig(plotsdir("native/sensitivity_exps/3.Pacific_Sens_Exps_ΔW.png"), bbox_inches = "tight")

fig, ax = plt.subplots(1, 3, figsize = (18, 5), sharey = true)
cf = []
Ψ_bounds = 0.3
levels = -Ψ_bounds:0.005:Ψ_bounds
for (i, expt) in enumerate(["only_wind", "only_init", "only_kappa"])
    ax[3].set_ylim(500, 6000)

    ax[i].set_title(labels[i], fontweight = "bold")
    predi = -(ΨBol[expt][:, 1:end-2, :] .- ΨBol[expt][:, 3:end, :])
    mean_Ψ =  1e-6 .* mean(predi, dims = 3)[:,:,1]
    CM = ax[i].contourf(ϕ_avg[2:end-1], z, mean_Ψ,levels = levels, cmap=cmo.balance, 
    vmin = -Ψ_bounds, vmax = Ψ_bounds, extend = "both")
    ax[i].invert_yaxis()
    push!(cf, CM)

    ax[i].set_xlabel("Latitude", fontweight = "bold")
end
ax[1].set_ylabel("Depth [m]", fontweight = "bold")
fig.subplots_adjust(wspace = 0.1)
levels = -Ψ_bounds:0.1:Ψ_bounds
cbar = fig.colorbar(cf[1], ax= ax[:], orientation = "vertical", fraction = 0.04, label = L"[\mathbf{\Delta W^* }]")
cbar.set_ticks(levels)
fig.savefig(plotsdir("native/sensitivity_exps/3.Pacific_Sens_Exps_ΔW_bol.png"), bbox_inches = "tight")
fig