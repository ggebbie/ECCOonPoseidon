include("../../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");
@pyimport seaborn as sns;

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)
region = "NPAC"; 

uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= -z[:].<= uplvl)
tecco = 1992+1/24:1/12:2018

vars =   ["iter0_bulkformula", "only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]

wθ_resid = Dict(); wθ_eul = Dict(); wθ_bol = Dict()

for (i, expname) in enumerate(vars)
    fname = datadir("native/" * expname * region * "_THETA_budget_Wθ_eul_bol_" * suffix * ".jld2")
    print(expname)
    vars = jldopen(fname)["dθ"]

    wθfull = vars["wθBot"] .- vars["wθTop"]
    wθbol= vars["wθBot_bol"] .- vars["wθTop_bol"]
    wθ_resid[expname] = 1.0 .* wθfull * (100 * 86400 * 365 * 10)
    wθ_bol[expname] = 1.0 .* wθbol * (100 * 86400 * 365 * 10)

    wθ_eul[expname] = wθ_resid[expname] .- wθ_bol[expname]

end


effect_exp = ["only_wind", "only_kappa", "only_init"]

[wθ_bol[expt] .-= wθ_bol["iter0_bulkformula"] for expt in effect_exp]
[wθ_eul[expt] .-= wθ_eul["iter0_bulkformula"] for expt in effect_exp]

# wθ_resid["SUM"] = wθ_resid["only_init"] .+ wθ_resid["only_kappa"] .+ wθ_resid["only_sfc"] .- (2 .*wθ_resid["iter0_bulkformula"])
# wθ_eul["SUM"] = wθ_eul["only_init"] .+ wθ_eul["only_kappa"] .+ wθ_eul["only_sfc"] .- (2 .*wθ_eul["iter0_bulkformula"])
# wθ_bol["SUM"] = wθ_bol["only_init"] .+ wθ_bol["only_kappa"] .+ wθ_bol["only_sfc"] .- (2 .*wθ_bol["iter0_bulkformula"])

fig, axes = plt.subplots(2, 1, figsize = (7., 10))
[ax.set_xlabel("time", fontweight = "bold") for ax in axes]
exps =  ["iter0_bulkformula", "iter129_bulkformula", "only_init", "only_kappa", "only_wind"]

axes[1].set_title("Eulerian Vertical Advection Convergence", fontweight = "bold")
axes[2].set_title("Bolus Vertical Advection Convergence", fontweight = "bold")

for (i, expname) in enumerate(effect_exp)
    # println(expname)
    lw = 3

    axes[1].plot(tecco, wθ_eul[expname], label = plot_labels_effects[expname], 
    linewidth = lw, color = exp_colors[expname])
    axes[2].plot(tecco, wθ_bol[expname], label = plot_labels_effects[expname], 
    linewidth = lw, color = exp_colors[expname])


    # axes[1].legend(markerscale = 0.1, ncols = 1, columnspacing = 0.5, frameon = false)
    # axes[3].legend(markerscale = 0.1, ncols = 1, columnspacing = 0.5, frameon = false)

end
axes[1].legend(markerscale = 0.1, ncols = 1, columnspacing = 0.5, frameon = false)

axes[2].legend(markerscale = 0.1, ncols = 1, columnspacing = 0.5, frameon = false)

axes[1].set_ylabel("[cK per decade]", fontweight = "bold")
axes[2].set_ylabel("[cK per decade]", fontweight = "bold")

axes[1].set_ylim(-2.3, 2.3)
axes[2].set_ylim(-0.5, 0.5)

# axes[2].set_ylim(-0.03, 0.03)
fig.savefig(plotsdir("native/sensitivity_exps/4.HeatBudgetVertEulBolAdvectionFluxes.png"))

fig