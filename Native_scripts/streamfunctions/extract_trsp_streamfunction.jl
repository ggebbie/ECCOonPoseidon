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
# Ψ_exp_labels = [L"\Psi^{\Delta F, \Delta T }", L"\Psi^{\Delta F}",
#                 L"\Psi^{\Delta T}", L"\Psi^{0}"]

@time for (key,values) in shortnames
    expname = key
    println(expname)
    @time Ψ_exp_mean[expname] = extract_meridionalΨ(expname,diagpath, Γ, γ, msk)
end
jldsave(datadir("ΨwBolusmean_"*region*".jld2"); Ψ_exp_mean)