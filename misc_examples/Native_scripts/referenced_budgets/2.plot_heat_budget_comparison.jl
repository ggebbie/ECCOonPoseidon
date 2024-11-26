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
colors = colors[[1, 4, 3, 2]]
#blue, red, #green, orange

sns.set_theme(context = "notebook", style = "ticks", font_scale = 1.1,
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

fig, axes = plt.subplots(1, 4, figsize = (16, 5), sharex = true, sharey = true)
axes[1].set_title("Mid-Depth North Pacific \n Temperature")
axes[2].set_title("Advective Heat Flux \n Convergence Contribution")
axes[3].set_title("Diffusive Heat Flux \n Convergence Contribution")
axes[4].set_title("Geothermal Heat Flux \n Contribution")

[ax.set_ylabel("[cK]") for ax in axes]
[ax.set_xlabel("time") for ax in axes]
fig.tight_layout()
integrate(start, x) = cumsum([start, x...])[1:end-1]

lw = 2.5; α = 0.8
labels = ["Iteration 129", "Iteration 0"]
zorders = [1, 0]
for (i, expname) in enumerate(["iter129_bulkformula", "iter0_bulkformula"])
    fname = datadir("native/" * expname * region * "_THETA_budget_ref_" * suffix * ".jld2")
    vars = load(fname)["dθ"]
    GTFandκ = vars["κxyθ"] .+ vars["GTH"] .+ vars["κzθ"]
    uθ = integrate(0, vars["VθSouth"])
    ∇wθ = integrate(0, vars["wθBot"] .- vars["wθTop"])
    thetaend = 100 .*(vars["θ"] .- vars["θ"][1])
    println(thetaend[end])
    println((100 .* (uθ .+ ∇wθ) .* 2.628e+6)[end])

    axes[1].plot(tecco, 100 .* (vars["θ"] .- vars["θ"][1] ), label = labels[i], c = colors[i], alpha =α,  linewidth = lw, zorder = zorders[i])
    axes[2].plot(tecco, 100 .* (uθ .+ ∇wθ) .* 2.628e+6 , label = labels[i], c = colors[i], alpha = α, linewidth = lw, zorder = zorders[i])
    axes[3].plot(tecco, 100 .* integrate(0, vars["κxyθ"] .+ vars["κzθ"]).* 2.628e+6, label = labels[i], c = colors[i], alpha = α, linewidth = lw, zorder = zorders[i])
    axes[4].plot(tecco, 100 .* integrate(0, vars["GTH"]).* 2.628e+6, label = labels[i], c = colors[i], alpha = α, linewidth = lw, zorder = zorders[i])

end
[ax.legend(loc = "lower center", frameon = false) for ax in axes]
fname = region * "_THETA_INT_ref_iter129anditer0" * suffix * ".png"
fig.savefig(plotsdir("native/generals/" * fname), dpi = 400, bbox_inches = "tight")
fig