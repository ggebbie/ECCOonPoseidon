#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../../src/intro.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
Printf, LaTeXStrings, PyCall, GibbsSeaWater, Dierckx, Interpolations

import PyPlot as plt, NaNMath as nm
@pyimport seaborn as sns
@pyimport pandas as pd
@pyimport cmocean.cm as cmo

include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());


include(srcdir("plot_and_dir_config.jl"))

region = "NPAC"
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths);


tecco= 1992+1/24:1/12:2018 # ecco years

p₀ = 2000; nz = 50
P = MeshArray(γ,Float32, nz)
for ijk in eachindex(P)
    P[ijk] .= pstdz[ijk[2]]
end
function get_densities(expname)
    filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
    datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    nt = length(datafilelist_θ[1:312]); nz = 50
    σz_mean = zeros(Float32, nt, 50, 270)
    σ = MeshArray(γ,Float32,nz)
    for tt in 1:nt
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
        θname = datafilelist_θ[tt]

        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        @time θSz = γ.read(diagpath[expname]*θname,MeshArray(γ,Float32,2*nz))
        θz = θSz[:, 1:nz]; Sz = θSz[:, nz+1:end]
        for ijk in eachindex(P)
            σtemp = densityJMD95.(θz.f[ijk],Sz.f[ijk], P[ijk], p₀) #EOS from MITGCM 
            σ.f[ijk] .= σtemp .- 1000
        end 
            
        σz_mean[tt, :, :] .= zonal_average(σ, cell_volumes)

    end
    return σz_mean
end

σz_mean = Dict()

for expname in ["iter0_bulkformula", "only_kappa"]
    σz_mean[expname] = get_densities(expname)
end

ϕ_zonal = zonal_average(ϕ, PAC_msk)

σz_0 = 1 .* permutedims(σz_mean["iter0_bulkformula"], (2, 3, 1))
Δz = -z[2:end] .- -z[1:end-1]
N_0 = σz_0[2:end, :, :] .- σz_0[1:end-1, :, :] ./ Δz
N_0 = 9.81 .* N_0 ./ 1029
N_0 = mean(N_0, dims = 3)[:, :, 1]
zF = (z[2:end] .+ z[1:end-1]) ./ 2

σz_129 = 1 .* permutedims(σz_mean["only_kappa"], (2, 3, 1))
N_129 = σz_129[2:end, :, :] .- σz_129[1:end-1, :, :] ./ Δz
N_129 = 9.81 .* N_129 ./ 1029
N_129 = mean(N_129, dims = 3)[:, :, 1]

fig,ax=plt.subplots(figsize = (10, 7.5))
not_nans = (!isnan).(ϕ_zonal)
var = 100 .* (N_0[:, not_nans] .- N_129[:, not_nans]) ./N_0[:, not_nans]
vmax = nm.maximum(abs.(var)); 
cms = ax.contourf(ϕ_zonal[not_nans], zF,  var, levels = LinRange(0, vmax, 100),
vmin = -vmax, vmax = vmax, cmap = cmo.balance)
ax.invert_yaxis()
fig.colorbar(cms)
fig



σz_time_mean = deepcopy(σz_mean)
[σz_time_mean[key] = mean(σz_mean[key], dims = 1)[:] for key in keys(σz_mean)]
dρdz(ρ, Δz) = (ρ[2:end] .- ρ[1:end-1]) ./ Δz
N(ρ, Δz) = sqrt.(-10 * dρdz(ρ, Δz) / 1029)
[σz_time_mean[key] = N(σz_time_mean[key], Δz) for key in keys(σz_mean)]

key = "only_kappa"
ax.plot(100 * (σz_time_mean[key] .- σz_time_mean["iter0_bulkformula"]) ./ σz_time_mean["iter0_bulkformula"], zF, label = key)
fig
E, F = trend_matrices(tecco)

100 * (sum(σz_time_mean[key] .* Δz) .- sum(σz_time_mean["iter0_bulkformula"] .* Δz)) ./ sum(σz_time_mean["iter0_bulkformula"] .* Δz)
lin_trend(x, t) = (sum(F[2, :] .* x) .* (t .- mean(t))) .+ sum(F[1, :] .* x) 
tecco = collect(tecco)


NPAC_σ = vcat(collect(32:0.25:36.5), collect(36.7:0.01:37), collect(37.03:0.02:37.1))
NPAC_σ = sort(unique(NPAC_σ))
σ_levels = NPAC_σ
vmin, vmax = extrema(σ_levels[34:end])


expts = ["iter0_bulkformula", "only_init", "only_kappa", "only_sfc", "iter129_bulkformula"]
σP_d = Dict()
for (i, expt) in enumerate(expts)
    σP = zeros(312)
    for tt in 1:nt
        σtoP = Spline1D(σz_mean[expt][tt, 35:end], pstdz[35:end]; k = 3)
        σP[tt] = σtoP(36.95)
    end
    σP_d[expt] = σP
end

iter0_trend = lin_trend(σP_d["iter0_bulkformula"], tecco)

fig,ax=plt.subplots(1, sharey = true, figsize = (10, 7.5))
ax.set_title("Pressure Level for σ = 36.95")
for (i, expt) in enumerate(["iter0_bulkformula", "only_init", "only_kappa", "only_sfc", "iter129_bulkformula"])
    ax.plot(tecco[1:nt], σP_d[expt], color = exp_colors[expt], label = plot_labels[expt])

end


ax.set_ylabel("dbar")
fig
#  [ax.set_ylim(2200, 3000) for ax in axs]
# axs[1].invert_yaxis()
fig.savefig(plotsdir("native/sensitivity_exps/NPAC_single_isopycnal_timeseries.png"))


fig,ax=plt.subplots(1, sharey = true, figsize = (10, 7.5))
ax.set_title("Pressure Level for σ = 36.95")
for (i, expt) in enumerate(["iter0_bulkformula", "only_init", "only_kappa", "only_sfc", "iter129_bulkformula"])
    ax.plot(tecco[1:nt], σP_d[expt] .- iter0_trend, color = exp_colors[expt], label = plot_labels[expt])

end


ax.set_ylabel("dbar")
fig
#  [ax.set_ylim(2200, 3000) for ax in axs]
# axs[1].invert_yaxis()
fig.savefig(plotsdir("native/sensitivity_exps/NPAC_single_isopycnal_timeseries_iter0trendremoved.png"))

