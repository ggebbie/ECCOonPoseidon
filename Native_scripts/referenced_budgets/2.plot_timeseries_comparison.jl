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
sns.set_theme(context = "talk", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));
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

fig, axes = plt.subplots(1, 1, figsize = (5, 5), sharex = true)
axes.set_title("Net θ"); axes.set_ylabel("[cK]"); axes.set_xlabel("time")
lw = 2
for (i, expname) in enumerate(["iter129_bulkformula", "iter0_bulkformula", "seasonalclimatology", "seasonalclimatology_iter0"])
    fname = datadir("native/" * expname * region * "_THETA_budget_ref_" * suffix * ".jld2")
    vars = load(fname)["dθ"]
    axes.plot(tecco, 100 .* (vars["θ"] .- vars["θ"][1] ), label = expname, c = colors[i], alpha = 0.7, linewidth = lw)
end
axes.legend(loc = "lower left")
fig.suptitle("Equatorial/South Pacific Temperature [z = 2-3km]")
fig.tight_layout()
fname = region * "_THETA_INT_ref_" * suffix * ".png"
fig.savefig(plotsdir("native/" * fname), dpi = 400)
fig