include("../../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson,
    BenchmarkTools, LaTeXStrings, PyCall
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");

#blue, red, #green, orange

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)
region = "NPAC"; 

uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= z[:].<= uplvl)
tecco = 1992+1/24:1/12:2018

include(srcdir("plot_and_dir_config.jl"))

lw = 2.5; α = 0.8

exps =  ["iter129_bulkformula", "iter0_bulkformula", "only_wind", "only_init", "only_kappa", "only_buoyancy"]
advection = Dict(); diffustion = Dict(); GTF = Dict()
for (i, expname) in enumerate(exps)
    fname = datadir("native/" * expname * region * "_THETA_budget_ref_with_Bolus" * suffix * ".jld2")
    vars = load(fname)["dθ"]
    uθ = vars["VθSouth"]
    ∇wθ = vars["wθBot"] .- vars["wθTop"]
    oneyear = 365 * 86400 
    ck = 100
    advection[expname] = (uθ .+ ∇wθ) .* (ck * oneyear)
    diffustion[expname] = (vars["κxyθ"] .+ vars["κzθ"]) .* (ck * oneyear)
    GTF[expname] = (vars["GTH"]) .* (ck * oneyear)
end
effect_exp = ["only_init", "only_kappa", "only_wind", "only_buoyancy"]
[advection[expt] .-= advection["iter0_bulkformula"] for expt in effect_exp]
[diffustion[expt] .-= diffustion["iter0_bulkformula"] for expt in effect_exp]
[GTF[expt] .-= GTF["iter0_bulkformula"] for expt in effect_exp]


fig, axes = plt.subplots(1, 3, figsize = (16, 5), sharex = true, sharey = true)
axes[1].set_title("Advective Heat Flux \n Convergence Contribution")
axes[2].set_title("Diffusive Heat Flux \n Convergence Contribution")
axes[3].set_title("Geothermal Heat Flux \n Contribution")

[ax.set_ylabel("[cK]") for ax in axes]
[ax.set_xlabel("time") for ax in axes]
fig.tight_layout()
integrate(start, x) = cumsum([start, x...])[1:end-1]

for (i, expname) in enumerate(effect_exp)
    println(expname)

    axes[1].plot(tecco, advection[expname], label = plot_labels[expname], c = exp_colors[expname], linewidth = lw)
    axes[2].plot(tecco, diffustion[expname], label = plot_labels[expname], c = exp_colors[expname])
    axes[3].plot(tecco, GTF[expname], label = plot_labels[expname], c = exp_colors[expname], linewidth = lw)
end
fig
[ax.legend(loc = "lower center", frameon = false) for ax in [axes[2]]]
fname = region * "_THETA_INT_ref_iter129anditer0" * suffix * ".png"
# fig.savefig(plotsdir("native/generals/" * fname), dpi = 400, bbox_inches = "tight")
fig

fig.savefig(plotsdir("native/sensitivity_exps/1.HeatBudgetFluxes.png"))
