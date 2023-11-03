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

exps =  ["iter0_bulkformula", "only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]
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

fig, axes = plt.subplots(1, 2, figsize = (12, 7))
axes[1].set_title(L"\mathbf{A}_Z")
axes[2].set_title(L"\Delta \mathbf{A}_Z")

axes[1].set_ylabel("[cK]", fontweight = "bold")
[ax.set_xlabel("time", fontweight = "bold") for ax in axes]
fig.tight_layout()

integrate(start, x) = cumsum([start, x...])[1:end-1] .* (100 * 2.628e+6)

lw = 2.5; α = 0.8


sens_exp = ["only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]
ΔVadvection = deepcopy(Vadvection)
[ΔVadvection[nexp] .-= Vadvection["iter0_bulkformula"] for nexp in sens_exp]
ΔVadvection["SUM"] = 0 .* Vadvection["iter0_bulkformula"]
[ΔVadvection["SUM"] .+= ΔVadvection[nexp] for nexp in sens_exp[1:end-1]]

exps =  ["iter0_bulkformula", "only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]
plot_labels["SUM"] = "SUM"
exp_colors["SUM"] = "grey"
for (i, expname) in enumerate(exps)
    axes[1].plot(tecco, integrate(0, Vadvection[expname]), label = plot_labels[expname], c = exp_colors[expname], linewidth = lw)
    if expname != "iter0_bulkformula"
        axes[2].plot(tecco, integrate(0, ΔVadvection[expname]), label = plot_labels[expname], c = exp_colors[expname], linewidth = lw)
    end
end
axes[2].plot(tecco, integrate(0, ΔVadvection["SUM"]), label = plot_labels["SUM"], c = exp_colors["SUM"], linewidth = lw)
axes[1].legend(loc = "lower left", frameon = false, bbox_to_anchor = (-0.5, -0.4), ncols = 5)
fig.subplots_adjust(wspace = 0.5)
fig
fig.savefig(plotsdir("native/sensitivity_exps/1.HeatBudgetTimeSeries_simple.png"), bbox_inches = "tight")
