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
sns.set_theme(context = "notebook", style = "ticks",
              palette = sns.color_palette("colorblind"));

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
 
ocean_mask = OHC_helper.wet_pts(Γ)
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

Ψ_exp_mean = jldopen(datadir("ΨwBolusmean_"*region*".jld2"))["Ψ_exp_mean"]
Ψ_exp_mean = Dict(key => reverse(value', dims = 1) for (key, value) in Ψ_exp_mean)

XΨ_mask = XΨ[XΨ .> -40]
@pyimport matplotlib.patches as patches

# θ_bounds = maximum(abs.(extrema(θ_zonal_mean_anom0)))
# Ψ_bounds = 1e-6.* maximum(abs.(extrema(Ψ_exp_mean_anom0)))
# Ψlevels = LinRange(-Ψ_bounds, Ψ_bounds, 13)
# θlevels = LinRange(-θ_bounds, θ_bounds, 25)
Ψ_exp_mean
cp = Vector{Any}(missing, 1)  
# anom_label = [sym * " minus CTRL" for sym in labels]
fig,axes=plt.subplots(4,1, figsize = (12, 28))
alabels = ["Iteration 129", "Iteration 0", "Initial 129", "Initial 0"]
clevels = ["-9.0", "-7.5", "-1.5"]
for (i, expname) in enumerate(["iter129_bulkformula", "iter0_bulkformula", "seasonalclimatology", "seasonalclimatology_iter0"])
    ax = axes[i]
    ax.set_facecolor("black")
    Ψ_bounds = round(1e-6 .* maximum(abs.(extrema(Ψ_exp_mean))))
    Ψ_bounds = 10.5
    levels = -Ψ_bounds:1.5:Ψ_bounds
    # levels = sort(vcat(levels, [-0.5, 0.5]))
    Ψ_exp = Ψ_exp_mean[expname][:, XΨ .> -40]
    ax.contourf(XΨ_mask, abs.(z[:]), 1e-6.* Ψ_exp, cmap=cm.delta,levels = levels, 
    vmin = -1.4*Ψ_bounds, vmax = 1.4*Ψ_bounds, extend = "both")
    cs2 = ax.contour(XΨ_mask, abs.(z[30:46]), 1e-6.* Ψ_exp[30:46, :], colors="k",levels = levels, 
    vmin = -Ψ_bounds, vmax = Ψ_bounds)
    ax.contour(XΨ_mask, abs.(z[1:30]), 1e-6.* Ψ_exp[1:30, :], colors="k",levels = levels, 
    vmin = -Ψ_bounds, vmax = Ψ_bounds)
    ax.contour(XΨ_mask, abs.(z[46:end]), 1e-6.* Ψ_exp[46:end, :], colors="k",levels = levels, 
    vmin = -Ψ_bounds, vmax = Ψ_bounds)
    labels = ax.clabel(cs2, fontsize=20, inline=true, fmt = "%.1f", 
    inline_spacing = 12, rightside_up = true, use_clabeltext = true)
    for label in labels
        text = label.get_text()
        label.set_fontsize(20)  # Adjust the fontsize as desired
        if text in clevels
            label.set_fontweight("bold")
            label.set_fontsize(22)  # Adjust the fontsize as desired
            label.set_zorder(100)    # Set the zorder value
        end
    end
    rect = patches.Rectangle((23, 2000), 59, 1000, linewidth=3, edgecolor="red",facecolor="none", alpha = 0.7)
    ax.add_patch(rect)
    rect = patches.Rectangle((2.5, 4000), 0.1, 1500, linewidth=3, edgecolor="yellow",facecolor="none", alpha = 0.7)
    ax.add_patch(rect)
    ax.invert_yaxis()
    ax.set_xticks(-40:10:60)
    ax.set_xlim(-34, 60)
    ax.set_title(alabels[i])
    ax.set_ylabel("Depth [m]")
    lab = string.(abs.(collect(-40:10:60)))
    lab = lab .* ["°S", "°S", "°S", "°S", "", "°N", "°N", "°N", "°N", "°N", "°N"]
    ax.set_xticklabels(lab)
end

fig
fig.tight_layout()
fig
fig.savefig(plotsdir("native/generals/ΨwBolus_iter129_iter0_" * region * ".png"), dpi = 400, 
bbox_inches = "tight")

Ψ_exp = jldopen(datadir("Ψtimeseries_"*region*".jld2"))["Ψ_exp"]
Ψ_exp = Dict(key => reverse(permutedims(value, (2, 1, 3)), dims = 1) for (key, value) in Ψ_exp)
Ψ_exp = Dict(key => -mean(value[44:46, 80:90, :], dims = [1,2]) for (key, value) in Ψ_exp) #draws a circle around the maximum overturning
fig,axes=plt.subplots(1, 3, figsize = (10, 8), sharex = true, sharey = true)
for (i, expname) in enumerate(["iter129_bulkformula", "iter0_bulkformula", "seasonalclimatology"])
    axes[i].plot(tecco, 1e-6 .* Ψ_exp[expname][1, 1, :], label = expname)
    axes[i].set_title(expname)
end
fig
# # Create a Rectangle patch
# fig,ax=plt.subplots(1,1, figsize = (12, 7))


# # Add the patch to the Axes
# ax.invert_yaxis()
# ax.set_xticks(-40:10:60)
# ax.set_yticks(1000:500:4000)
# ax.set_ylim(1000, 4000)
# ax.set_xlim(-34, 60)



fig