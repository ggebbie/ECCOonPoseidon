using Pkg
Pkg.activate(".")
include("../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings
using Statistics
import NaNMath as nm
using PyCall
import PyPlot as plt
using NaNStatistics
include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

include(srcdir("plot_and_dir_config.jl"))
 
ocean_mask = wet_pts(Γ)
region = "PAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "false", include_bering = true)

tecco = 1992+1/24:1/12:2018

vars =  ["only_init", "only_kappa", "only_sfc", "iter129_bulkformula",  "iter0_bulkformula"]
vars =  ["only_buoyancy", "only_wind"]

uplvl = -2.4e3; botlvl = -2.6e3; suffix = "2to3"
lvls = findall( botlvl .<= -z[:].<= uplvl)

diagpath["mean_tau_noadjust_redo_bf"] = vastdiagdir("seasonalclimatology", "run_only_clim_iter0_tau_bf")
diagpath["mean_tau_yesadjust_redo_bf"] = vastdiagdir("seasonalclimatology", "run_only_clim_iter129_tau_bf2")

lon=[i for i=-179.:1:179., j=-89.:1:89.]
lat=[j for i=-179.:1:179., j=-89.:1:89.]
(f,i,j,w)=InterpolationFactors(Γ,vec(lon),vec(lat))

# (lon, lat, f,i,j,w) = load(datadir("0.5deg_interpolation_factors.jld2"), "lon", "lat", "f", "i", "j", "w")

Wres_reg = Dict()

Wres_reg["mean_tau_noadjust_redo_bf"] = zeros(size(lon)..., 312)
Wres_reg["mean_tau_yesadjust_redo_bf"] = zeros(size(lon)..., 312)

exps = ["mean_tau_noadjust_redo_bf", "mean_tau_yesadjust_redo_bf"]
for expname in exps
    fname = datadir("native/" * expname * "_W_residual_2500m.jld2")
    Wres = load(fname)["Wres"]
    print(size(Wres))
    for tt in 1:312
        W_interp = Interpolate(Wres[:, tt] .* 1,f,i,j,w) #interpolate using half-degree resolution
        W_interp = reshape(W_interp, size(lon));
        W_interp[W_interp .== 0.0] .= NaN
        Wres_reg[expname][:, :, tt] .= W_interp

    end
end

Wres_reg["Difference"] = Wres_reg["mean_tau_yesadjust_redo_bf"] .- Wres_reg["mean_tau_noadjust_redo_bf"]
# Wres_reg["Difference"] = (Wres_reg["mean_tau_yesadjust_redo_bf"] .- Wres_reg["mean_tau_noadjust_redo_bf"])

region = "NPAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "false", include_bering = true)
# PAC_msk = wet_pts(Γ)
interp_mask = Interpolate(PAC_msk,f,i,j,w) #interpolate using half-degree resolution
interp_mask = reshape(interp_mask, size(lon));
interp_mask[interp_mask .> 0.2] .= 1
interp_mask[interp_mask .<= 0] .= NaN

fig, ax = plt.subplots(sharey = true)
ax.contourf(lon, lat, nanmean(Wres_reg["Difference"] .* interp_mask, dims = 3)[:, :, 1])
fig



normalize(x) = (x .- mean(x)) ./ std(x)
loc1 = Base.findmin( ((40 .- lat).^2) .+ ((-150 .- lon).^2) )[2]
loc2 = Base.findmin( ((40 .- lat).^2) .+ ((150 .- lon).^2) )[2]

fig, ax = plt.subplots()
ax.plot(normalize(Wres_reg["Difference"][loc1, :]))
ax.plot(normalize(Wres_reg["Difference"][loc2, :]))
fig 

output = ax.xcorr(normalize(Wres_reg["Difference"][loc2, :]), 
                normalize(Wres_reg["Difference"][loc1, :]), 
maxlags=120, usevlines = false, normed=true, alpha=0.0)

fig, ax = plt.subplots(figsize = (6.5, 5))
ax.plot(output[1] ./ 12, output[2], c = "k", marker = "o", markersize = 5)
ax.set_xlabel("lag (years)")
ax.set_ylabel(L"\Delta w^{res}" * " correlation")
ax.set_title("Lagged correlation between\ntwo points in Pacific Ocean\nalong 45°N")
ax.text(0.2, -0.25, "150°W leads", fontsize = 15, 
                  transform=ax.transAxes, horizontalalignment = "center", fontweight = "bold")
ax.text(0.8, -0.25, "150°E leads", fontsize = 15, 
                  transform=ax.transAxes, horizontalalignment = "center", fontweight = "bold")
output[1][Base.findmax(output[2])[2]] / 12
fig.savefig(plotsdir("wind_rewrite/7.w_cross_correlations.png"), bbox_inches = "tight", dpi = 400)

fig

