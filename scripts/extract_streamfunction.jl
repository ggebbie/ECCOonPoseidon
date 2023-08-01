
include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyPlot, PyCall
using .OHC_helper

include(srcdir("config_exp.jl"))
@pyimport cmocean.cm as cmo

(ϕ,λ) = latlonC(γ)

runpath,diagpath = listexperiments(exprootdir());

ocean_mask = wet_pts(Γ)
region = "PAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "full")
msk = PAC_msk;

#low memory footprint MeshArray(γ, Int8, 180)
#define experimennt
expname = "iter129_bulkformula"

Ψ_exp = extract_meridional_Ψ_Mean_EulBolus(expname,diagpath, Γ, γ, PAC_msk)

#plotting
XΨ=collect(-89.0:89.0); Y=reverse(z); #coordinate variables
fig,ax=plt.subplots(figsize = (12, 7))
Ψ_mean = reverse(Ψ_exp', dims = 1)
Ψ_bounds = 10.5
levels = -Ψ_bounds:1.5:Ψ_bounds
ax.contourf(XΨ, abs.(z[:]), 1e-6.* Ψ_mean, cmap=cmo.delta,levels = levels, 
vmin = -1.4*Ψ_bounds, vmax = 1.4*Ψ_bounds, extend = "both")
ax.invert_yaxis()
ax.set_xticks(-40:10:60)
ax.set_xlim(-34, 60)
fig