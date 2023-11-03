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
sns.set_theme(context = "notebook", style = "ticks",
              palette = colors, rc = custom_params);

fig, axes = plt.subplots(1, 2, figsize = (11, 5), sharex = true, sharey = true)
axes[1].set_title("Horizontal Advective Heat Flux \n Convergence Contribution")
axes[2].set_title("Vertical Advective Heat Flux \n Convergence Contribution")
axes[1].set_ylabel("[cK]")
[ax.set_xlabel("time") for ax in axes]
fig.tight_layout()
integrate(start, x) = cumsum([start, x...])[1:end-1]

lw = 2.5; α = 0.8

exps =  ["iter0_bulkformula", "only_init", "only_kappa", "only_wind"]
for (i, expname) in enumerate(exps)
    println(expname)
    fname = datadir("native/" * expname * region * "_THETA_budget_ref_with_Bolus" * suffix * ".jld2")
    vars = load(fname)["dθ"]
    GTFandκ = vars["κxyθ"] .+ vars["GTH"] .+ vars["κzθ"]
    uθ = integrate(0, vars["VθSouth"])
    ∇wθ = integrate(0, vars["wθBot"] .- vars["wθTop"])
    thetaend = 100 .*(vars["θ"] .- vars["θ"][1])
    println(thetaend[end])
    println((100 .* (uθ .+ ∇wθ) .* 2.628e+6)[end])

    axes[1].plot(tecco, 100 .* uθ .* 2.628e+6 , label = plot_labels[expname], c = exp_colors[expname], linewidth = lw)
    axes[2].plot(tecco, 100 .* ∇wθ .* 2.628e+6, label = plot_labels[expname], c = exp_colors[expname], linewidth = lw)
end
fig
[ax.legend(loc = "lower left", frameon = false) for ax in axes]
fname = region * "_THETA_INT_ref_iter129anditer0" * suffix * ".png"
# fig.savefig(plotsdir("native/generals/" * fname), dpi = 400, bbox_inches = "tight")
fig.savefig(plotsdir("native/sensitivity_exps/2.HeatBudgetVertHorizAdvection.png"))

fig

# fig.savefig(plotsdir("native/sensitivity_exps/HeatBudget.png"))
