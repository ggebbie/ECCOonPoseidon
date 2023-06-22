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

include(srcdir("config_exp.jl"))
@pyimport cmocean.cm as cm
@pyimport seaborn as sns;
@pyimport pandas as pd;

subtract_dict(D, subkey) = Dict(key => value .- D[subkey] for (key, value) in D)
mask_dict(D, lvls) = Dict(key => value[lvls, :] for (key, value) in D)

colors =  sns.color_palette("deep")[1:4]
labels = ["ECCO", "ECCO Forcing", "ECCO Init. Cond. and κ", "CTRL"]
sns.set_theme(context = "talk", style = "white", font_scale = 1.0,
              palette = sns.color_palette("deep"));

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
marks = expsymbols()
 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
                            region, include_bering = false)
msk = PAC_msk;
(msk == ocean_mask) && (region = "GLOB")
uplvl = Inf; botlvl = -Inf
lvls = findall( botlvl .<= z[:].<= uplvl)

cell_depths = get_cell_depths(msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths);
smush_depths = smush(cell_depths); smush_depths[findall(smush_depths .== 0)] = Inf
inv_depths = 1 ./ smush_depths

tecco = collect(1992+1/24:1/12:2018)
#first 3 yrs tecco[1:36]
#last 3 yrs tecco[end-36+1:end]
# vcat(1:36,312-36+1:312)
XΨ=collect(-89.0:89.0); Y=reverse(z); #coordinate variables
Xθ = OHC_helper.ma_zonal_sum(ϕ .* area) ./ OHC_helper.ma_zonal_sum(area)

insuffix = "sfctobot"

Xθ = OHC_helper.ma_zonal_sum(ϕ .* area) ./ OHC_helper.ma_zonal_sum(area)

fig,ax=plt.subplots(1,1, figsize = (20, 12.5))
i = 1
expname = "iter0_bulkformula"


cs = ax.contourf(Xθ, abs.(Y), 0.0 .* reverse(θ_zonal_mean[expname], dims = 1), cmap="RdBu_r", vmin = -1, vmax = 1)
ax.invert_yaxis()
ax.set_facecolor("#964B00")
# ax[i].set_ylim(-3500, -1500)
ax.set_xticks(-40:10:70)
ax.set_title(anom_label[i])
ax.set_ylabel("Depth [m]")
ax.set_xlabel("Latitude")
ax.set_xlim(-50, 70)
fig
fig.suptitle("<ψ> (Solid); <θ> (Colored) in North Pacific")

fig.tight_layout()

fig
# fig.savefig(plotsdir() * "/GCC_Plots/TrspNθAnomiter0_" * region * "_" * outsuffix * ".png",
# dpi = 1000)
