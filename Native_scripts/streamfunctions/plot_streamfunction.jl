include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings
using .OHC_helper
import NaNMath as nm
import PyPlot as plt
using PyCall

include(srcdir("config_exp.jl"))
@pyimport cmocean.cm as cm
@pyimport seaborn as sns;
@pyimport pandas as pd;

subtract_dict(D, subkey) = Dict(key => value .- D[subkey] for (key, value) in D)
mask_dict(D, lvls) = Dict(key => value[lvls, :] for (key, value) in D)

colors =  sns.color_palette("deep")[1:4]
labels = ["ECCO", "ECCO Forcing", "ECCO Init. Cond. and κ", "CTRL"]
sns.set_theme(context = "poster", style = "white", font_scale = 1.0,
              palette = sns.color_palette("deep"));

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
marks = expsymbols()
 
ocean_mask = wet_pts(Γ)
region = "PAC"; 
tecco = collect(1992+1/24:1/12:2018)
#first 3 yrs tecco[1:36]
#last 3 yrs tecco[end-36+1:end]
vcat(1:36,312-36+1:312)
XΨ=collect(-89.0:89.0); Y=reverse(z); #coordinate variables
Xθ = OHC_helper.ma_zonal_sum(ϕ .* area) ./ OHC_helper.ma_zonal_sum(area)

insuffix = "sfctobot"
outsuffix = "1tobot"
uplvl = -1e3; botlvl = -Inf
lvls = findall( botlvl .<= z[:].<= uplvl)
Y_mask = reverse(z[lvls])
lat, lon = (ϕ .- (-30), λ .- (-140))
findmin((lat).^2 + (lon).^2)
λ[4][20932], ϕ[4][20932]

Ψ_exp_mean = jldopen(datadir("Ψmean_"*region*".jld2"))["Ψ_exp_mean"]
Ψ_exp_mean = Dict(key => reverse(value', dims = 1) for (key, value) in Ψ_exp_mean)

XΨ_mask = XΨ[XΨ .> -40]

cp = Vector{Any}(missing, 1)  
# anom_label = [sym * " minus CTRL" for sym in labels]
for expname in keys(shortnames)
    fig,ax=plt.subplots(1,1, figsize = (20, 12.5))
    Ψ_bounds = round(1e-6.* maximum(abs.(extrema(Ψ_exp_mean))))
    Ψ_bounds = 9
    levels = -Ψ_bounds:1.5:Ψ_bounds
    Ψ_exp = Ψ_exp_mean[expname][:, XΨ .> -40]
    # cp[1] = ax[i].contourf(Xθ, Y_mask, reverse(θ_zonal_exp, dims = 1), vmin = -θ_bounds, vmax = θ_bounds,
    # cmap = cm.balance)
    cs = ax.contourf(XΨ_mask, abs.(z[:]), 1e-6.* Ψ_exp, cmap=cm.delta,levels = levels, 
    vmin = -Ψ_bounds-4, vmax = Ψ_bounds+4, extend = "both")
    cs = ax.contour(XΨ_mask, abs.(z[:]), 1e-6.* Ψ_exp, colors="k",levels = levels, 
    vmin = -Ψ_bounds, vmax = Ψ_bounds)
    ax.clabel(cs, fontsize=20, inline=true, fmt = "%.1f")
    ax.invert_yaxis()
    ax.set_xticks(-40:10:60)
    ax.set_xlim(-39, 60)
    ax.set_title("Pacific Streamfunction in ECCO")
    fig.savefig(plotsdir("native/generals/streamfunction" * expname *  region * ".png"), dpi = 400)
end