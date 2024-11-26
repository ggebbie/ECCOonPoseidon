include("../../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson,
    BenchmarkTools, LaTeXStrings, PyCall, DataFrames
import NaNMath as nm
import PyPlot as plt
import NumericalIntegration

include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");
gs = pyimport("matplotlib.gridspec");
# gridspec = gs.GridSpec
#blue, red, #green, orangeDataFrames

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)
region = "NPAC"; 

uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= z[:].<= uplvl)
tecco = 1992+1/24:1/12:2018

include(srcdir("plot_and_dir_config.jl"))

exps =  ["iter0_bulkformula", "iter129_bulkformula", "only_buoyancy", 
"only_init", "only_kappa", "only_wind"]
Hadvection = Dict(); diffusion = Dict(); GTF = Dict()
Vadvection = Dict()
ext_advection = Dict()
int_advection = Dict()

temps = Dict()
for (i, expname) in enumerate(exps)
    fname = datadir("native/" * expname * region * "_THETA_budget_ref_with_Bolus_wextra" * suffix * ".jld2")
    vars = load(fname)["dθ"]
    temps[expname] = vars["θ"]
    uθ = vars["VθSouth"]
    ∇wθ = vars["wθBot"] .- vars["wθTop"]
    Hadvection[expname] = (uθ)
    Vadvection[expname] = (∇wθ)
    ext_advection[expname] = uθ .+ ∇wθ
    int_advection[expname] = vars["wθ"] .+ vars["uvθ"] .- ext_advection[expname]

    diffusion[expname] = (vars["κxyθ"] .+ vars["κzθ"])
    GTF[expname] = (vars["GTH"])
end

integrate(start, x) = cumsum([start, x...])[1:end-1] .* (100 * 2.628e+6)

sns.set_style("darkgrid", Dict("axes.facecolor" => ".95", "grid.color" => "0.5"))

fig, axs = plt.subplots(1, 2, figsize = (12, 6), sharey = true)
axs[1].set_title("Mid-depth North Pacific\nTemperature Anomaly\n" * L"\theta'")
axs[2].set_title("Heat Budget")
[ax.set_xlabel("time") for ax in axs]
# fig.tight_layout()


lw = 2.5; α = 0.8

exps =  ["iter0_bulkformula", "iter129_bulkformula", "only_buoyancy", 
"only_init", "only_kappa", "only_wind"]
plot_labels_list = ["ECCO V4r4 First-Guess\n(Iteration 0)", "ECCO V4r4 Official Release\n(Iteration 129)", 
"Adjusted\nBuoyancy", "Adjusted\nI.C.", "Adjusted\nMixing", 
"Adjusted\nWind"]
plot_labels["Difference"] = "Difference"
exp_colors["Difference"] = "k"
fig.tight_layout()
alphas = [0.9, 0.9, 0.4, 0.4, 0.4, 0.9]
for (i, expname) in enumerate(exps)
    println(expname)

    axs[1].plot(tecco, 100 .* (temps[expname] .- temps[expname][1]), label = plot_labels_list[i], c = exp_colors[expname], linewidth = lw, alpha = alphas[i])
    axs[2].plot(tecco, integrate(0, int_advection[expname]), label = plot_labels_list[i], c = exp_colors[expname], linewidth = lw, alpha = alphas[i])
    axs[2].plot(tecco, integrate(0, diffusion[expname]), label = plot_labels_list[i], c = exp_colors[expname], alpha = alphas[i])
    axs[2].plot(tecco, integrate(0, GTF[expname]), label = plot_labels_list[i], c = exp_colors[expname], linewidth = lw, alpha = alphas[i])
    axs[2].plot(tecco, integrate(0, ext_advection[expname]), label = plot_labels_list[i], c = exp_colors[expname], linewidth = lw, zorder = 10, alpha = alphas[i])    

end
axs[1].tick_params(which="both", bottom=true, left = true)
axs[2].tick_params(which="both", bottom=true)

fig_labs = uppercase.(["a", "b", "c", "d", "e"])
for (i, a) in enumerate(axs)
    a.annotate(fig_labs[i], (0.05, 0.05), fontsize = 20, 
    xycoords="axes fraction", fontweight = "bold")
end

axs[2].annotate("Advection\n" * L"\mathbf{A}^{}", (0.5, 0.15), fontsize = 12.5, 
xycoords="axes fraction", ha="center", fontweight = "normal")
axs[2].annotate("Diffusion\n"* L"\mathbf{F}_\kappa", (0.5, 0.85), fontsize = 12.5, 
xycoords="axes fraction", ha="center", fontweight = "normal")
axs[2].annotate("Residual\n" * L"\mathbf{R}", (0.8, 0.57), fontsize = 12.5, 
xycoords="axes fraction", ha="center", fontweight = "normal")
axs[2].annotate("Geothermal\n" * L"\mathbf{F}_{geo}", (0.81, 0.75), fontsize = 12.5, 
xycoords="axes fraction", ha="center", fontweight = "normal")
fig.subplots_adjust(wspace = 0.05)
[a.grid(alpha = 0.2) for a in axs]
axs[1].legend(loc = "lower left", frameon = false, bbox_to_anchor = (0.05, 0-.7), ncols = 3)

fig.savefig(plotsdir("wind_rewrite/1.1.HeatBudgetTimeSeries_complex_all_2panels.png"), bbox_inches = "tight", dpi = 400)
fig