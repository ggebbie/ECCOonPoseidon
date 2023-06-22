include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings
using .OHC_helper
using Plots
using ColorSchemes
import NaNMath as nm
import PyPlot as plt
using PyCall

@pyimport cmocean.cm as cm
@pyimport seaborn as sns;
@pyimport pandas as pd;

subtract_dict(D, subkey) = Dict(key => value .- D[subkey] for (key, value) in D)
mask_dict(D, lvls) = Dict(key => value[lvls, :] for (key, value) in D)

colors =  sns.color_palette("deep")[1:4]
sns.set_theme(context = "talk", font_scale = 1.0,
              palette = sns.color_palette("deep"));
sns.set_style("white")
include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
marks = expsymbols()
 
ocean_mask = wet_pts(Γ)
region = "PAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "full")
msk = PAC_msk;
tecco = 1992+1/24:1/12:2018

X=collect(-89.0:89.0); Y=reverse(z); #coordinate variables

Ψ_exp_mean = Dict();
Ψ_exp_labels = [L"\Psi^{\Delta F, \Delta T }", L"\Psi^{\Delta F}",
                L"\Psi^{\Delta T}", L"\Psi^{0}"]

@time for (key,values) in shortnames
    expname = key
    println(expname)
    Ψ_exp_mean[expname] = extract_meridionalΨ̄(expname,diagpath, Γ, γ, msk)
end
jldsave(datadir("Ψmean_"*region*".jld2"); Ψ_exp_mean)

# subtract_dict(D, subkey) = Dict(key => value .- D[subkey] for (key, value) in D)
# mask_dict(D, lvls) = Dict(key => value[lvls, :] for (key, value) in D)

# Ψ_exp_mean = Dict(key => reverse((value[2, :, :] .- value[1, :, :])', dims = 1) for (key, value) in Ψ_exp_mean)
# # Ψ_exp_mean = mask_dict(Ψ_exp_mean, lvls)
# Ψ_exp_mean_anom0 = subtract_dict(Ψ_exp_mean, "iter0_bulkformula") 

# pygui(true)
# fig,axs=plt.subplots(4,1, figsize = (7.5, 15), sharex = "col")
# Ψ_bounds = 1e-6.* maximum(abs.(extrema(Ψ_exp_mean)))
# Ψlevels = LinRange(-Ψ_bounds, Ψ_bounds, 35)
# cs = Vector{Any}(missing, 1)

# @time for (i, expname) in enumerate(keys(shortnames))
#     Ψ_exp = Ψ_exp_mean[expname]
#     cs[1] = axs[i].contourf(X, Y, reverse(1e-6.* Ψ_exp, dims = 1), 
#                     levels = Ψlevels, cmap = cm.delta, inline = true)
#     c = axs[i].contour(X, Y, reverse(1e-6.* Ψ_exp, dims = 1), 
#                     levels = Ψlevels, colors = "k")
#     axs[i].clabel(c, fontsize=15, inline=true, fmt = "%0.2f")

#     axs[i].set_title(Ψ_exp_labels[i])
# end
# [ax.set_xlim(-35, 65) for ax in axs]
# [ax.set_xticks(-40:10:60) for ax in axs]
# fig.tight_layout()
