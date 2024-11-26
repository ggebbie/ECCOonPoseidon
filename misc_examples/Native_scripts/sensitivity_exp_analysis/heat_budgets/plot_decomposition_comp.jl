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

runpath,diagpath = listexperiments(exprootdir());
fig, axes = plt.subplots(1, 2, figsize = (10, 5), sharex = true, sharey = true)
[ax.set_ylabel("[cK]") for ax in axes]
[ax.set_xlabel("time") for ax in axes]
vars =  ["only_init", "only_kappa", "only_sfc", "iter129_bulkformula",  "iter0_bulkformula", "only_buoyancy", "only_wind"]

# alabels = ["Iteration 129", "Iteration 0"]
for (i, expname) in enumerate(vars)

    #load in the decomposition
    fname = datadir("native/" * expname * region * "_THETA_budget_Wθ_eul_bol_" * suffix * ".jld2")
    print(expname)
    vars = jldopen(fname)["dθ"]
    wθest = vars["wθBot"] .- vars["wθTop"]
    fname = datadir("native/" * expname * region * "_THETA_budget_ref_with_Bolus" * suffix * ".jld2")
    vars = jldopen(fname)["dθ"]
    wθtrue = vars["wθBot"] .- vars["wθTop"]

    nt = length(wθest)
    lw = 2.5
    integrate(start, x) = cumsum([start, x...])[1:end-1]
    # axes[2].set_title("Bottom Boundary")

    # axes[1].set_title("Top Boundary")
    # axes[1].text(0.6, 0.92, alabels[i], transform=axes[1].transAxes, fontweight = "bold")
    axes[1].plot(tecco[1:nt], 100 .* integrate(0, wθest) .* 2.628e+6, label = expname, alpha = 0.8, linewidth = lw, zorder = 0)
    axes[2].plot(tecco[1:nt], 100 .* integrate(0, wθtrue[1:nt]) .* 2.628e+6, alpha = 1, label = expname, linewidth = lw)
    
    axes[1].set_title("True Vertical Advection Convergence")
    axes[2].set_title("Estimated Vertical Advection Convergence")

    axes[1].legend(markerscale = 0.1, ncols = 1, columnspacing = 0.5, frameon = false)
    axes[2].legend(markerscale = 0.1, ncols = 1, columnspacing = 0.5, frameon = false)

end
fig.tight_layout()
fig