include("../../src/intro.jl")

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
sns.set_theme(context = "talk", style = "ticks",
              palette = colors, rc = custom_params);

exps =  ["iter129_bulkformula", "iter0_bulkformula"]
Hadvection = Dict(); diffusion = Dict(); GTF = Dict()
Vadvection = Dict()

for (i, expname) in enumerate(exps)
    fname = datadir("native/" * expname * region * "_THETA_budget_ref_with_Bolus" * suffix * ".jld2")
    vars = load(fname)["dθ"]
    uθ = vars["VθSouth"]
    ∇wθ = vars["wθBot"] .- vars["wθTop"]
    Hadvection[expname] = (uθ)
    Vadvection[expname] = (∇wθ)

    diffusion[expname] = (vars["κxyθ"] .+ vars["κzθ"])
    GTF[expname] = (vars["GTH"])
end

fig, axes = plt.subplots(1, 4, figsize = (16, 5), sharex = true, sharey = true)
axes[1].set_title(L"\mathbf{A}_Z")
axes[2].set_title(L"\mathbf{A}_H")
axes[3].set_title(L"\mathbf{F}_\kappa")
axes[4].set_title(L"\mathbf{F}_{geo}")
axes[1].set_ylabel("[cK]", fontweight = "bold")
[ax.set_xlabel("time", fontweight = "bold") for ax in axes]
fig.tight_layout()

integrate(start, x) = cumsum([start, x...])[1:end-1] .* (100 * 2.628e+6)

lw = 2.5; α = 0.8

exps =  ["iter0_bulkformula", "iter129_bulkformula"]
for (i, expname) in enumerate(exps)
    println(expname)
    axes[1].plot(tecco, integrate(0, Vadvection[expname]), label = plot_labels[expname], c = exp_colors[expname], linewidth = lw)
    println("Vertical")
    println(round(integrate(0, Vadvection[expname])[end], digits = 2))
    axes[2].plot(tecco, integrate(0, Hadvection[expname]), label = plot_labels[expname], c = exp_colors[expname], linewidth = lw)
    println("Horizontal")
    println(round(integrate(0, Hadvection[expname])[end], digits = 2))
    axes[3].plot(tecco, integrate(0, diffusion[expname]), label = plot_labels[expname], c = exp_colors[expname])
    println("Diff")
    println(round(integrate(0, diffusion[expname])[end], digits = 2))
    axes[4].plot(tecco, integrate(0, GTF[expname]), label = plot_labels[expname], c = exp_colors[expname], linewidth = lw)
    println("Geo")
    println(round(integrate(0, GTF[expname])[end], digits = 2))

end
axes[2].legend(loc = "lower left", frameon = false, bbox_to_anchor = (0.1, -0.4), ncols = 5)
fig.subplots_adjust(wspace = 0.1)
fig
fig.savefig(plotsdir("native/paper_figures/1.HeatBudgetTimeSeries_simple.png"), bbox_inches = "tight")
