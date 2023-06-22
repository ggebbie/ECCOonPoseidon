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
=
alabels = ["Iteration 129", "Iteration 0"]
for (i, expname) in enumerate(["iter129_bulkformula", "iter0_bulkformula"])
    fname = datadir("native/" * expname * region * "_THETA_budget_ref_Reynolds_" * suffix * ".jld2")
    vars = load(fname)["dθ"]
    axes = axess[:, i]
    println(expname)
    
    wθTop = vars["Wbθb_Top"] .+ vars["Wbθp_Top"] .+ vars["Wpθb_Top"] .+ vars["Wpθp_Top"]
    wθBot = vars["Wbθb_Bottom"] .+ vars["Wbθp_Bottom"] .+ vars["Wpθb_Bottom"] .+ vars["Wpθp_Bottom"]
    nt = length(wθBot)
    lw = 2.5
    integrate(start, x) = cumsum([start, x...])[1:end-1]
    # axes[2].set_title("Bottom Boundary")
    println("Top")
    println(1 - (var(integrate(0, wθTop .- vars["Wbθb_Top"])) / var(integrate(0, wθTop))))
    println("Bot")
    println(1 - (var(integrate(0, wθBot .- vars["Wbθb_Bottom"])) / var(integrate(0, wθBot))))

end



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