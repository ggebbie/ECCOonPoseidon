#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
using .OHC_helper
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");
@pyimport seaborn as sns;

colors =  sns.color_palette("deep")[1:4]
exp_colors = colors[[1, 4, 3, 2]]
sns.set_theme(context = "talk", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));
(ϕ,λ) = latlonC(γ)
area = readarea(γ)
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

ocean_mask = OHC_helper.wet_pts(Γ)
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));
H = OHC_helper.smush(cell_depths, γ); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H

uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= z[:].<= uplvl)
tecco = 1992+1/24:1/12:2018

runpath,diagpath = listexperiments(exprootdir());
fig, axess = plt.subplots(2, 2, figsize = (10, 10), sharex = true, sharey = true)
[ax.set_ylabel("[cK]") for ax in axess]
[ax.set_xlabel("time") for ax in axess]
colors = cat("black", sns.color_palette()[5:6], sns.color_palette()[end], "grey", dims = 1)
alabels = ["Iteration 129", "Iteration 0"]
for (i, expname) in enumerate(["iter129_bulkformula", "iter0_bulkformula"])
    fname = datadir("native/" * expname * region * "_THETA_budget_ref_Reynolds_" * suffix * ".jld2")
    vars = load(fname)["dθ"]
    axes = axess[:, i]
    print(expname)
    
    wθTop = vars["Wbθb_Top"] .+ vars["Wbθp_Top"] .+ vars["Wpθb_Top"] .+ vars["Wpθp_Top"]
    wθBot = vars["Wbθb_Bottom"] .+ vars["Wbθp_Bottom"] .+ vars["Wpθb_Bottom"] .+ vars["Wpθp_Bottom"]
    nt = length(wθBot)
    lw = 2.5
    integrate(start, x) = cumsum([start, x...])[1:end-1]
    # axes[2].set_title("Bottom Boundary")

    # axes[1].set_title("Top Boundary")
    axes[1].text(0.4, 0.92, alabels[i], transform=axes[1].transAxes, fontweight = "bold")
    axes[1].plot(tecco[1:nt], -100 .* integrate(0, wθTop) .* 2.628e+6, label = L"\mathcal{A}_T (w \theta)", alpha = 0.5, linewidth = 2*lw, c = exp_colors[i], zorder = 0)
    axes[1].plot(tecco[1:nt], -100 .* integrate(0, vars["Wbθb_Top"]) .* 2.628e+6, alpha = 1, label = L"\mathcal{A}_T (\overline{w} \overline{\theta})", linewidth = lw, c = exp_colors[i])
    axes[1].plot(tecco[1:nt], -100 .* integrate(0, vars["Wbθp_Top"]) .* 2.628e+6, linewidth = lw, c = colors[2], alpha = 0.3, linestyle = "--")
    axes[1].plot(tecco[1:nt], -100 .* integrate(0, vars["Wpθb_Top"]) .* 2.628e+6, linewidth = lw, c = colors[3], alpha = 0.3, linestyle = "--")
    axes[1].plot(tecco[1:nt], -100 .* integrate(0, vars["Wpθp_Top"]) .* 2.628e+6, linewidth = lw, c = colors[4], alpha = 0.3, linestyle = "--")
    axes[1].legend(markerscale = 0.1, ncols = 1, columnspacing = 0.5, frameon = false)

    axes[2].text(0.43, 0.92, alabels[i], transform=axes[2].transAxes, fontweight = "bold")
    axes[2].plot(tecco[1:nt], 100 .* integrate(0, wθBot) .* 2.628e+6, label = L"\mathcal{A}_B (w \theta)", alpha = 0.5, linewidth = 2*lw, c = exp_colors[i], zorder = 0)
    axes[2].plot(tecco[1:nt], 100 .* integrate(0, vars["Wbθb_Bottom"]) .* 2.628e+6, alpha =1 , label = L"\mathcal{A}_B (\overline{w} \overline{\theta})", linewidth = lw,  c = exp_colors[i])
    axes[2].plot(tecco[1:nt], 100 .* integrate(0, vars["Wbθp_Bottom"]) .* 2.628e+6, linewidth = lw, c = colors[2], alpha = 0.3, linestyle = "--")
    axes[2].plot(tecco[1:nt], 100 .* integrate(0, vars["Wpθb_Bottom"]) .* 2.628e+6, linewidth = lw, c = colors[3], alpha = 0.3, linestyle = "--")
    axes[2].plot(tecco[1:nt], 100 .* integrate(0, vars["Wpθp_Bottom"]) .* 2.628e+6, linewidth = lw, c = colors[4], alpha = 0.3, linestyle = "--")
    axes[2].legend(markerscale = 0.1, ncols = 1, columnspacing = 0.5, frameon = false)

end
fig.tight_layout()
fig
# sns.move_legend(axes[2], "lower center",  bbox_to_anchor=(0, -0.25), ncol=5, frameon=true, borderaxespad=0.)
fname = region * "_THETA_INT_budget_ref_reynolds_iter129iter0_hidden" * suffix * ".png"
fig.savefig(plotsdir("native/generals/" * fname), dpi = 400)
# fig


# fig, axes = plt.subplots(2, 2, figsize = (20, 5), sharex = true)
# axes[1, 1].set_title("Meridional Heat Flux")
# axes[1, 2].set_title("Vertical Heat Flux")
# axes[1].set_ylabel("[cK / century]")
# axes[2].set_ylabel("[cK / century]")

# [ax.set_xlabel("time") for ax in axes[2, :]]
# nt = length(tot_flux)
# println(tot_flux)
# axes[1, 1].plot(tecco[1:nt], vars["VθSouth"] .* (3.154e+9 * 100), label = L"Vθ_{in}", c = colors[i])
# axes[1, 1].plot(tecco[1:nt], -vars["VθNorth"].* (3.154e+9 * 100), label = L"-Vθ_{out}", c = colors[i], alpha = 0.5)
# axes[2, 1].plot(tecco[1:nt], (vars["VθSouth"] .-vars["VθNorth"]) .* (3.154e+9 * 100), label = "∇(Vθ)", c = "k", alpha = 0.7)

# axes[1, 2].plot(tecco[1:nt], vars["wθBot"] .* (3.154e+9 * 100), label = L"Wθ_{in}", c = colors[i])
# axes[1, 2].plot(tecco[1:nt], -vars["wθTop"] .* (3.154e+9 * 100), label = L"-Wθ_{out}", c = colors[i], alpha = 0.5)
# axes[2, 2].plot(tecco[1:nt], (vars["wθBot"].-vars["wθTop"]) .* (3.154e+9 * 100), label = L"∇(Wθ)", c = "k", alpha = 0.7)

# # axes[2, 2].legend(ncol = 3)
# [ax.legend() for ax in axes]
# # sns.move_legend(axes[3], "lower center",  bbox_to_anchor=(0, -0.34), ncol=4, frameon=true, borderaxespad=0.)
# fig