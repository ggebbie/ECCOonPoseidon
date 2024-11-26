include("../../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))

# cmo = pyimport("cmocean.cm");
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

Δwθ = Dict()
wΔθ = Dict()
ΔwΔθ = Dict()

integrate(start, x) = cumsum([start, x...])[1:end-1]
include(srcdir("plot_and_dir_config.jl"))

sns.set_theme(context = "talk", style = "ticks",
              palette = colors, rc = custom_params);

for (i, expname) in enumerate(vars)
    fname = datadir("native/" * expname * region * "_THETA_budget_Wθ_eul_bol_ΔDecomp" * suffix * ".jld2")
    print(expname)
    vars = jldopen(fname)["dθ"]

    wθfull = vars["wθBot"] .- vars["wθTop"]
    wθ_resid[expname] = 100 .* integrate(0, wθfull).* 2.628e+6
    Δwθ[expname] = 100 .* integrate(0, vars["ΔW_θBot"] .- vars["ΔW_θTop"]).* 2.628e+6
    wΔθ[expname] = 100 .* integrate(0, vars["W_ΔθBot"] .- vars["W_ΔθTop"]).* 2.628e+6
    ΔwΔθ[expname] = 100 .* integrate(0, vars["ΔW_ΔθBot"] .- vars["ΔW_ΔθTop"]).* 2.628e+6

end

vars =   [ "only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]

[wθ_resid[nexpt] .-= wθ_resid["iter0_bulkformula"] for nexpt in vars]


fig, axs = plt.subplots(3, 1, figsize = ( 12, 14))
[ax.set_xlabel("time", fontweight = "bold") for ax in axs]

vars =   ["only_init", "only_kappa", "only_wind"]

for (i, expname) in enumerate(vars)
    # println(expname)
    lw = 2.5

    axs[i].plot(tecco, wθ_resid[expname], label = L"\Delta \mathbf{A}_Z (w_{res}, \theta^\#)", 
    linewidth = lw, color = exp_colors[expname])
    axs[i].plot(tecco, Δwθ[expname], label = L" \mathbf{A}_Z (\Delta w_{res}, \theta^\#)", 
    linewidth = lw, color = 0.5 .* exp_colors[expname], alpha = 0.9)
    axs[i].plot(tecco, wΔθ[expname], label = L" \mathbf{A}_Z (w_{res}, \Delta \theta^\#)", 
    linewidth = lw, color = exp_colors[expname], alpha = 0.6)
    axs[i].plot(tecco, -ΔwΔθ[expname], label = L" \mathbf{A}_Z (\Delta w_{res}, \Delta \theta^\#)", 
    linewidth = lw, color = exp_colors[expname], alpha = 0.3)
    axs[i].legend(frameon = false, ncols = 4, loc = "lower center",  bbox_to_anchor=(0.5, 1))
    axs[i].axhline(0, c = "k", linestyle = "--")

    # axes[3].legend(markerscale = 0.1, ncols = 1, columnspacing = 0.5, frameon = false)

end
[ax.set_ylabel("[cK]", fontweight = "bold") for ax in axs]
fig.subplots_adjust(hspace = 0.5)

fig.savefig(plotsdir("native/sensitivity_exps/2.HeatBudgetVertΔDecomp.png"), bbox_inches = "tight")

fig