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
vars = ["iter129_bulkformula",  "iter0_bulkformula", "seasonalclimatology"]
i = 1; expname = "iter129_bulkformula"
fname = datadir("native/" * expname * region * "_THETA_budget_ref_TEST" * suffix * ".jld2")
vars = load(fname)["dθ"]

fig, axes = plt.subplots(2, 2, figsize = (20, 5), sharex = true, sharey = "row")
axes[1, 1].set_title("Meridional Heat Flux")
axes[1, 2].set_title("Vertical Heat Flux")
axes[1].set_ylabel("[cK]")
[ax.set_xlabel("time") for ax in axes[2, :]]
tot_flux = vars["GTH"] .+ vars["VθSouth"] .- vars["wθTop"] .+ vars["wθBot"]  .+ vars["κxyθ"] .+ vars["κzθ"] .- vars["VθNorth"]
nt = length(tot_flux)
println(tot_flux)
axes[1, 1].plot(tecco[1:nt], vars["VθSouth"], label = L"Vθ_{in}", c = colors[i])
axes[1, 1].plot(tecco[1:nt], -vars["VθNorth"], label = L"-Vθ_{out}", c = colors[i], alpha = 0.5)
axes[2, 1].plot(tecco[1:nt], vars["VθSouth"] .-vars["VθNorth"], label = "∇(Vθ)", c = "k", alpha = 0.7)
axes[2, 1].plot(tecco[1:nt], vars["VθConv"], label = "True ∇(Vθ)", c = "red", alpha = 0.7)

axes[1, 2].plot(tecco[1:nt], vars["wθBot"], label = L"Wθ_{in}", c = colors[i])
axes[1, 2].plot(tecco[1:nt], -vars["wθTop"], label = L"-Wθ_{out}", c = colors[i], alpha = 0.5)
axes[2, 2].plot(tecco[1:nt], vars["wθBot"].-vars["wθTop"], label = L"∇(Wθ)", c = "k", alpha = 0.7)
axes[2, 2].plot(tecco[1:nt], vars["WθConv"], label = "True ∇(Wθ)", c = "red", alpha = 0.7)
fig.suptitle(expname)
# axes[2, 2].legend(ncol = 3)
[ax.legend() for ax in axes]
# sns.move_legend(axes[3], "lower center",  bbox_to_anchor=(0, -0.34), ncol=4, frameon=true, borderaxespad=0.)
fig 

fig, axes = plt.subplots(1, figsize = (20, 5))
axes.plot(tecco[1:nt], cumsum(tot_flux) .* 2.628e+6, c="red", label = "approx")
axes.plot(tecco[1:nt], vars["θ"] .- vars["θ"][1], c="black", label = "true")
axes.legend()
fig
