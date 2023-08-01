include("../src/intro.jl")
include("../src/experimental.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, LaTeXStrings, PyPlot, PyCall
using .experimental
import NaNMath as nm
include(srcdir("config_exp.jl"))
@pyimport cmocean.cm as cmo

(ϕ,λ) = latlonC(γ); area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());

ocean_mask = wet_pts(Γ)
region = "PAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "full")

#define experimennt
expname = "iter129_bulkformula"

Ψ_exp, ϕ_avg = extract_meridional_Ψ_Mean(expname,diagpath, Γ, γ, PAC_msk)
Ψ_exp_timeseries, ϕ_avg = extract_meridional_Ψ(expname,diagpath, Γ, γ, PAC_msk)

Ψ_exp_timeseries_mean = mean(Ψ_exp_timeseries, dims = 3)[:, :, 1]

#plotting the streamfunction
fig,ax=plt.subplots(1, 2, figsize = (24, 7))
Ψ_bounds = 10.5
levels = -Ψ_bounds:1.5:Ψ_bounds
ax.contourf(ϕ_avg, abs.(z[:]), 1e-6.* Ψ_exp, cmap=cmo.delta,levels = levels, 
vmin = -1.4*Ψ_bounds, vmax = 1.4*Ψ_bounds, extend = "both")

ax.contourf(ϕ_avg, abs.(z[:]), 1e-6.* Ψ_exp_timeseries_mean, cmap=cmo.delta,levels = levels, 
vmin = -1.4*Ψ_bounds, vmax = 1.4*Ψ_bounds, extend = "both")

ax.invert_yaxis()
ax.set_xticks(-40:10:60)
ax.set_xlim(-34, 60)
fig
