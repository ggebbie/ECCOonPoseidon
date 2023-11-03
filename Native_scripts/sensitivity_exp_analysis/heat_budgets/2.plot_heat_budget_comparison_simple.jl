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
sns.set_theme(context = "talk", style = "ticks",
              palette = colors, rc = custom_params);

exps =  ["iter129_bulkformula", "iter0_bulkformula"]
Hadvection = Dict(); diffustion = Dict(); GTF = Dict()
Vadvection = Dict()
for (i, expname) in enumerate(exps)
    fname = datadir("native/" * expname * region * "_THETA_budget_ref_with_Bolus" * suffix * ".jld2")
    vars = load(fname)["dθ"]
    uθ = vars["VθSouth"]
    ∇wθ = vars["wθBot"] .- vars["wθTop"]
    Hadvection[expname] = (uθ)
    Vadvection[expname] = (∇wθ)

    diffustion[expname] = (vars["κxyθ"] .+ vars["κzθ"])
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
    
    axes[1].plot(tecco, integrate(0, Vadvection[expname]), label = plot_labels[expname], c = exp_colors[expname], linewidth = lw)
    axes[2].plot(tecco, integrate(0, Hadvection[expname]), label = plot_labels[expname], c = exp_colors[expname], linewidth = lw)
    axes[3].plot(tecco, integrate(0, diffustion[expname]), label = plot_labels[expname], c = exp_colors[expname])
    axes[4].plot(tecco, integrate(0, GTF[expname]), label = plot_labels[expname], c = exp_colors[expname], linewidth = lw)

end
axes[2].legend(loc = "lower left", frameon = false, bbox_to_anchor = (-0.5, -0.4), ncols = 5)
fig.subplots_adjust(wspace = 0.1)
fig
fig.savefig(plotsdir("native/sensitivity_exps/1.HeatBudgetTimeSeries_simple.png"), bbox_inches = "tight")


fig, axes = plt.subplots(1, 4, figsize = (16, 5), sharex = true, sharey = true)
axes[1].set_title(L"\Delta \mathbf{A}_Z")
axes[2].set_title(L"\Delta \mathbf{A}_H")
axes[3].set_title(L"\Delta \mathbf{F}_\kappa")
axes[4].set_title(L"\Delta \mathbf{F}_{geo}")

axes[1].set_ylabel("[cK]", fontweight = "bold")
[ax.set_xlabel("time", fontweight = "bold") for ax in axes]
fig.tight_layout()
integrate(start, x) = cumsum([start, x...])[1:end-1] .* (100 * 2.628e+6)

lw = 2.5; α = 0.8

exps =  ["only_init", "only_kappa", "only_wind", "only_buoyancy"]
for (i, expname) in enumerate(exps)
    
    axes[1].plot(tecco, integrate(0, Vadvection[expname] .- Vadvection["iter0_bulkformula"]), label = plot_labels[expname], c = exp_colors[expname], linewidth = lw)
    axes[2].plot(tecco, integrate(0, Hadvection[expname] .- Hadvection["iter0_bulkformula"]), label = plot_labels[expname], c = exp_colors[expname], linewidth = lw)
    axes[3].plot(tecco, integrate(0, diffustion[expname] .- diffustion["iter0_bulkformula"]), label = plot_labels[expname], c = exp_colors[expname])
    axes[4].plot(tecco, integrate(0, GTF[expname] .- GTF["iter0_bulkformula"]), label = plot_labels[expname], c = exp_colors[expname], linewidth = lw)

end
axes[2].legend(loc = "lower left", frameon = false, bbox_to_anchor = (-0.1, -0.4), ncols = 5)
fig.subplots_adjust(wspace = 0.1)

fig
fig.savefig(plotsdir("native/sensitivity_exps/1.ΔHeatBudgetTimeSeries.png"), bbox_inches = "tight")


