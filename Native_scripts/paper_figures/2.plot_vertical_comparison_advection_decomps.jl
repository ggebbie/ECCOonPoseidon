include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))
include(srcdir("plot_and_dir_config.jl"))

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
    wθ_resid[expname] = 100 .* integrate(0, wθfull).* 2.628f+6
    Δwθ[expname] = 100 .* integrate(0, vars["ΔW_θBot"] .- vars["ΔW_θTop"]).* 2.628f+6
    wΔθ[expname] = 100 .* integrate(0, vars["W_ΔθBot"] .- vars["W_ΔθTop"]).* 2.628f+6
    ΔwΔθ[expname] = 100 .* integrate(0, vars["ΔW_ΔθBot"] .- vars["ΔW_ΔθTop"]).* 2.628f+6

end
wθ_resid["iter129_bulkformula"] .-= wθ_resid["iter0_bulkformula"]

vars =   [ "only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]

fig, axs = plt.subplots(figsize = ( 7, 7))
axs.set_xlabel("time", fontweight = "bold")

expname =   "iter129_bulkformula"
lw = 2.5

axs.plot(tecco, wθ_resid[expname], label = L"\Delta \mathbf{A}_Z (w^{res}, \theta^\#)", 
linewidth = lw, color = exp_colors[expname])
axs.plot(tecco, Δwθ[expname] .-ΔwΔθ[expname] , label = L" \mathbf{A}_Z (\Delta w^{res}, \theta_0^\#)", 
linewidth = 3, color = 0.5 .* exp_colors[expname], alpha = 0.8)
axs.plot(tecco, wΔθ[expname] .-ΔwΔθ[expname] , label = L" \mathbf{A}_Z (w^{res}_0, \Delta \theta^\#)", 
linewidth = lw, color = exp_colors[expname], alpha = 0.6)
axs.plot(tecco, ΔwΔθ[expname], label = L" \mathbf{A}_Z (\Delta w^{res}, \Delta \theta^\#)", 
linewidth = lw, color = exp_colors[expname], alpha = 0.3)
axs.legend(frameon = false, ncols = 2, loc = "lower center",  bbox_to_anchor=(0.5, -0.35))
axs.axhline(0, c = "k", linestyle = "--", alpha = 0.2)

axs.set_ylabel("[cK]", fontweight = "bold")
fig.subplots_adjust(hspace = 0.5)
fig
fig.savefig(plotsdir("native/paper_figures/2.HeatBudgetVertΔDecomp.png"), bbox_inches = "tight")

fig


vars =   [ "only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]

Δwθ["SUM"] = sum([Δwθ[nexpt]  for nexpt in vars[1:end-1]])
Δwθ["eps"] = Δwθ["iter129_bulkformula"] .- Δwθ["SUM"]

ΔwΔθ["SUM"] = sum([ΔwΔθ[nexpt]  for nexpt in vars[1:end-1]])
ΔwΔθ["eps"] = ΔwΔθ["iter129_bulkformula"] .- ΔwΔθ["SUM"]

vars = vcat(vars,"SUM")
plot_labels["eps"] = "ϵ"
exp_colors["eps"] = "grey"
fig, axs = plt.subplots(1, figsize = ( 7, 7))
axs.set_xlabel("time", fontweight = "bold")
vars =   [ "only_init", "only_kappa", "only_wind", "only_buoyancy", "eps"]

axs.set_ylabel(L" \mathbf{A_Z (\Delta w^{res}, \theta_0^\#)}" * " [cK]", fontweight = "bold") 
for (i, expname) in enumerate(vars)
    println(expname)
    lw = 2.5
    tmp = Δwθ[expname].-ΔwΔθ[expname]
    println(tmp[end])

    axs.plot(tecco, tmp, label = plot_labels[expname], 
    linewidth = lw, color = exp_colors[expname])
    axs.legend(frameon = false, ncols = 3, loc = "lower center",  bbox_to_anchor=(0.5, 1))
    axs.axhline(0, c = "k", linestyle = "--")

    # axes[3].legend(markerscale = 0.1, ncols = 1, columnspacing = 0.5, frameon = false)

end
fig
fig.subplots_adjust(hspace = 0.5)

fig.savefig(plotsdir("native/paper_figures/2.HeatBudgetVertΔWDecomp.png"), bbox_inches = "tight")

fig


