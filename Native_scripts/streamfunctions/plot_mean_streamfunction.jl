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
@pyimport cmocean.cm as cmo
@pyimport seaborn as sns;
@pyimport matplotlib.patches as patches

sns.set_theme(context = "notebook", style = "ticks",
              palette = sns.color_palette("colorblind"));


(ϕ,λ) = latlonC(γ)

ocean_mask = OHC_helper.wet_pts(Γ)
region = "PAC"; 

XΨ=collect(-89.0:89.0); Y=reverse(z); #coordinate variables
Xθ = OHC_helper.ma_zonal_sum(ϕ .* area) ./ OHC_helper.ma_zonal_sum(area)

uplvl = -1e3; botlvl = -Inf
lvls = findall( botlvl .<= z[:].<= uplvl)
Y_mask = reverse(z[lvls])
XΨ_mask = XΨ[XΨ .> -40]
cp = Vector{Any}(missing, 1)  
# alabels = ["Iteration 129", "Iteration 0", "Initial 129", "Initial 0"]
clevels = ["-9.0", "-7.5", "-1.5"]
for (i, expname) in enumerate(myexps)
    fig,ax=plt.subplots(figsize = (12, 7))
    Ψ_exp = jldopen(datadir("ΨwBolustimeseries_"*region*"_" * expname *".jld2"))["Ψ_exp"]
    Ψ_mean = mean(Ψ_exp, dims = 3)[:, :, 1]
    Ψ_mean = reverse(Ψ_mean', dims = 1)
    ax.set_facecolor("black")
    Ψ_bounds = round(1e-6 .* maximum(abs.(extrema(Ψ_exp_mean))))
    Ψ_bounds = 10.5
    levels = -Ψ_bounds:1.5:Ψ_bounds
    Ψ_exp = Ψ_mean[:, XΨ .> -40]
    ax.contourf(XΨ_mask, abs.(z[:]), 1e-6.* Ψ_exp, cmap=cmo.delta,levels = levels, 
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
    ax.set_title(expname)
    ax.set_ylabel("Depth [m]")
    lab = string.(abs.(collect(-40:10:60)))
    lab = lab .* ["°S", "°S", "°S", "°S", "", "°N", "°N", "°N", "°N", "°N", "°N"]
    ax.set_xticklabels(lab)
    fig.savefig(plotsdir("native/generals/ΨwBolus_"*expname * "_" * region * ".png"), 
                dpi = 400, bbox_inches = "tight")
end

