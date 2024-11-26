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

i = 1; expname = "iter129_bulkformula"
fname = datadir("native/" * expname * region * "_THETA_budget_ref_" * suffix * ".jld2")
vars = load(fname)["dθ"]

GTFandκ = vars["κxyθ"] .+ vars["GTH"] .+ vars["κzθ"]
tot_flux =GTFandκ .+ vars["VθSouth"] .- vars["wθTop"] .+ vars["wθBot"] .- vars["VθNorth"]

fig, axes = plt.subplots(1, 4, figsize = (15, 5), sharex = true)
axes[1].set_title("Integrated Meridional Heat Flux"); axes[2].set_title("Integrated Vertical Heat Flux")
axes[3].set_title("Integrated Diffusive and \n Geothermal Heat Fluxes")
axes[4].set_title("Net θ"); 
[ax.set_ylabel("[cK]") for ax in axes]
[ax.set_xlabel("time") for ax in axes]

integrate(start, x) = cumsum([start, x...])[1:end-1]

nt = length(tot_flux)
lw = 2
axes[1].plot(tecco[1:nt], 100 .* integrate(0, vars["VθSouth"]) .* 2.628e+6, label = "Southern Boundary", c = colors[i], linewidth = lw)
axes[1].plot(tecco[1:nt], 100 .* integrate(0, -vars["VθNorth"]) .* 2.628e+6, label = "Northern Boundary", c = colors[i], alpha = 0.5, linewidth = lw)
axes[2].plot(tecco[1:nt], 100 .* integrate(0, vars["wθBot"]) .* 2.628e+6, label = "Bottom Boundary", c = colors[i], linewidth = lw)
axes[2].plot(tecco[1:nt], 100 .* integrate(0, -vars["wθTop"]) .* 2.628e+6, label = "Top Boundary", c = colors[i], alpha = 0.5, linewidth = lw)
axes[3].plot(tecco[1:nt], 100 .* integrate(0, GTFandκ) .* 2.628e+6, c = colors[i], linewidth = lw)
axes[4].plot(tecco[1:nt], 100 .* (vars["θ"] .- vars["θ"][1] ), label = "true θ", c = "k", alpha = 0.3, linewidth = 4*lw)
axes[4].plot(tecco[1:nt], 100 .* integrate(0, vcat(0, tot_flux[1:end-1])) .* 2.628e+6, label = "reconstruction", c = colors[i], alpha = 1, linewidth = lw)
[ax.legend(loc = "lower left") for ax in axes]
fig.suptitle("North Pacific Heat Budget for " * shortnames[expname] * " [z = 2-3km]")
fig.tight_layout()
fig
fname = expname * region * "_THETA_INT_budget_ref_" * suffix * ".png"
fig.savefig(plotsdir("native/" * fname), dpi = 400)


fig, axes = plt.subplots(1, 2, figsize = (20, 5), sharex = true, sharey = true)
axes[2].set_title("Vertical Heat Flux")
axes[1].set_title("Vertical Heat Flux")
axes[1].set_ylabel("[cK / century]")
axes[2].set_ylabel("[cK / century]")
[ax.set_xlabel("time") for ax in axes[2, :]]
nt = length(tot_flux)

# axes[1].plot(tecco[1:nt], vars["VθSouth"] .* (3.154e+9 * 100), label = "Southern Boundary", c = colors[i]);
# axes[1].plot(tecco[1:nt], -vars["VθNorth"].* (3.154e+9 * 100), label = "Northern Boundary", c = colors[i], alpha = 0.5)

axes[1].plot(tecco[1:nt], vars["wθBot"] .* (3.154e+9 * 100), label = "Bottom Boundary", c = colors[i])
axes[2].plot(tecco[1:nt], -vars["wθTop"] .* (3.154e+9 * 100), label = "Top Boundary", c = colors[i], alpha = 0.5)
[ax.legend() for ax in axes]
# axes[1].set_xlim(2000, 2003)
fig

