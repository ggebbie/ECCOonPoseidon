include("../../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson,
    BenchmarkTools, LaTeXStrings, PyCall, DataFrames
import NaNMath as nm
import PyPlot as plt
import NumericalIntegration

cumul_integrate = NumericalIntegration.cumul_integrate

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

diagpath["mean_tau_noadjust_redo_bf"] = vastdiagdir("seasonalclimatology", "run_only_clim_iter0_tau_bf")
diagpath["mean_tau_yesadjust_redo_bf"] = vastdiagdir("seasonalclimatology", "run_only_clim_iter129_tau_bf2")

exps =  ["iter0_bulkformula", "only_wind", "mean_tau_noadjust_redo_bf", "mean_tau_yesadjust_redo_bf"]
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
function integration(t, x)
    int_x = 3.154e+7* cumul_integrate(tecco, x)
    int_x .-= int_x[1]
    return 100 * int_x
end


lw = 2.5; α = 0.8

sns.set_style("darkgrid", Dict("axes.facecolor" => ".95", "grid.color" => "0.5"))

fig, axs = plt.subplots(1, 2, figsize = (12, 6), sharey = true)
axs[1].set_title("Mid-depth North Pacific\nTemperature Anomaly\n" * L"\theta'")
axs[2].set_title("Heat Budget")
[ax.set_xlabel("time") for ax in axs]

# exps =  ["iter0_bulkformula", "only_wind", "diffF", "mean_tau_noadjust_redo_bf", "mean_tau_yesadjust_redo_bf", "diffC"]
# plot_labels_list = ["Iteration 0", "Iteration 129", "Difference", "Iteration 0 (Climatological Winds)", 
# "Iteration 129 (Climatological Winds)", "Difference"]

exps =  ["iter0_bulkformula", "only_wind", "mean_tau_noadjust_redo_bf", "mean_tau_yesadjust_redo_bf"]
plot_labels_list = ["ECCO V4r4 First-Guess\n(Iteration 0)", "Adjusted Wind", 
"Seasonal First-Guess Wind", 
"Seasonal Adjusted Wind"]

exp_colors["mean_tau_noadjust_redo_bf"] = exp_colors["iter0_bulkformula"] 
exp_colors["mean_tau_yesadjust_redo_bf"] = exp_colors["only_wind"] 
exp_colors["diffC"] = "k"
exp_colors["diffF"] = "k"
fig.tight_layout()
E, F = trend_matrices(tecco)
trend(x) = (F * x[:])[2]

# alphas = [0.4, 0.4, 0.4, 0.9, 0.9, 0.9]
alphas = [0.4, 0.4, 0.9, 0.9]

temps["diffC"] = temps["mean_tau_yesadjust_redo_bf"] .- temps["mean_tau_noadjust_redo_bf"]
temps["diffF"] = temps["only_wind"] .- temps["iter0_bulkformula"]

int_advection["diffC"] = int_advection["mean_tau_yesadjust_redo_bf"] .- int_advection["mean_tau_noadjust_redo_bf"]
int_advection["diffF"] = int_advection["only_wind"] .- int_advection["iter0_bulkformula"]

diffusion["diffC"] = diffusion["mean_tau_yesadjust_redo_bf"] .- diffusion["mean_tau_noadjust_redo_bf"]
diffusion["diffF"] = diffusion["only_wind"] .- diffusion["iter0_bulkformula"]

GTF["diffC"] = GTF["mean_tau_yesadjust_redo_bf"] .- GTF["mean_tau_noadjust_redo_bf"]
GTF["diffF"] = GTF["only_wind"] .- GTF["iter0_bulkformula"]

ext_advection["diffC"] = ext_advection["mean_tau_yesadjust_redo_bf"] .- ext_advection["mean_tau_noadjust_redo_bf"]
ext_advection["diffF"] = ext_advection["only_wind"] .- ext_advection["iter0_bulkformula"]
println("advection trend:", 100 * trend(temps["diffC"]))
println("advection trend:", 100 * trend(temps["diffF"] ))

for (i, expname) in enumerate(exps)
    println(expname)
    axs[1].plot(tecco, 100 .* (temps[expname] .- temps[expname][1]), label = plot_labels_list[i], c = exp_colors[expname], linewidth = lw, alpha = alphas[i])
    println("total trend:", 100 * 100 * trend((temps[expname] .- temps[expname][1])))

    axs[2].plot(tecco, integration(0, int_advection[expname]), label = plot_labels_list[i], c = exp_colors[expname], linewidth = lw, alpha = alphas[i])
    println("residual trend:", 100 * trend(integration(0, int_advection[expname])))

    axs[2].plot(tecco, integration(0, diffusion[expname]), label = plot_labels_list[i], c = exp_colors[expname], alpha = alphas[i])
    println("diffusion trend:", 100 * trend(integration(0, diffusion[expname])))

    axs[2].plot(tecco, integration(0, GTF[expname]), label = plot_labels_list[i], c = exp_colors[expname], linewidth = lw, alpha = alphas[i])
    println("gtf trend:", 100 * trend(integration(0, GTF[expname])))

    axs[2].plot(tecco, integration(0, ext_advection[expname]), label = plot_labels_list[i], c = exp_colors[expname], linewidth = lw, zorder = 10, alpha = alphas[i])    
    println("advection trend:", 100 * trend(integration(0, ext_advection[expname])))
    println(" ")

end
axs[1].tick_params(which="both", bottom=true, left = true)
axs[2].tick_params(which="both", bottom=true)

fig_labs = uppercase.(["a", "b", "c", "d", "e"])
for (i, a) in enumerate(axs)
    a.annotate(fig_labs[i], (0.05, 0.05), fontsize = 20, 
    xycoords="axes fraction", fontweight = "bold")
end

fig
axs[2].annotate("Advection\n" * L"\mathbf{A}^{}", (0.45, 0.15), fontsize = 12.5, 
xycoords="axes fraction", ha="center", fontweight = "normal")
axs[2].annotate("Diffusion\n"* L"\mathbf{F}_\kappa", (0.5, 0.85), fontsize = 12.5, 
xycoords="axes fraction", ha="center", fontweight = "normal")
axs[2].annotate("Residual\n" * L"\mathbf{R}", (0.8, 0.55), fontsize = 12.5, 
xycoords="axes fraction", ha="center", fontweight = "normal")
axs[2].annotate("Geothermal\n" * L"\mathbf{F}_{geo}", (0.81, 0.75), fontsize = 12.5, 
xycoords="axes fraction", ha="center", fontweight = "normal")
fig.subplots_adjust(wspace = 0.05)
[a.grid(alpha = 0.5) for a in axs]
axs[1].legend(loc = "lower left", frameon = false, bbox_to_anchor = (0.1, 0-.55), ncols = 2)

fig.savefig(plotsdir("wind_rewrite/1.2.HeatBudgetTimeSeries_complex_all_2panels_clim.png"), bbox_inches = "tight", dpi = 400)
fig