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

exps =  ["mean_tau_noadjust_redo_bf", "mean_tau_yesadjust_redo_bf"]

#define dictionaries

Hadvection = Dict(); 
diffusion = Dict(); 
GTF = Dict()
Vadvection = Dict()
Topadvection = Dict()
Botadvection = Dict()
temps = Dict()

#fill dictionaries
for (i, expname) in enumerate(exps)
    fname = datadir("native/" * expname * region * "_THETA_budget_ref_with_Bolus_wextra" * suffix * ".jld2")
    vars = load(fname)["dθ"]
    temps[expname] = vars["θ"]
    uθ = vars["VθSouth"]
    ∇wθ = vars["wθBot"] .- vars["wθTop"]
    Hadvection[expname] = (uθ)
    Vadvection[expname] = (∇wθ)
    Botadvection[expname] = vars["wθBot"]
    Topadvection[expname] = -vars["wθTop"]

    diffusion[expname] = (vars["κxyθ"] .+ vars["κzθ"])
    GTF[expname] = (vars["GTH"])
end

integrate(start, x) = cumsum([start, x...])[1:end-1] .* (100 * 2.628e+6)
function integration(t, x)
    int_x = 3.154e+7* cumul_integrate(tecco, x)
    int_x .-= int_x[1]
    return 100 * int_x
end
E, F = trend_matrices(tecco)
trend(x) = (F * x[:])[2]

sns.set_style("darkgrid", Dict("axes.facecolor" => ".95", "grid.color" => "0.5"))

fig, axs = plt.subplots(1, 3, figsize = (12, 6), sharex = false, sharey = true)
axs[3].set_title(L"\mathbf{A}_{\mathcal{S}}")
axs[2].set_title(L"\mathbf{A}_{\mathcal{T}}")
axs[1].set_title(L"\mathbf{A}_{\mathcal{B}}")

axs[1].set_ylabel("[cK]", fontweight = "bold")
[ax.set_xlabel("time", fontweight = "bold") for ax in axs]
fig.tight_layout()

lw = 2.5; α = 0.9
plot_labels_list = [
"Seasonal First-Guess Wind", 
"Seasonal Adjusted Wind"]

# plot_labels_list = ["Iteration 0 (Climatological Winds)", "Iteration 129 (Climatological Winds)"]
E, F = trend_matrices(tecco)
trend(x) = (F * x[:])[2]
trends_dict = Dict()
exp_colors["mean_tau_noadjust_redo_bf"] = exp_colors["iter0_bulkformula"] 
exp_colors["mean_tau_yesadjust_redo_bf"] = exp_colors["only_wind"] 

for (i, expname) in enumerate(exps)
    trends_dict[expname] = zeros(3)
    println(expname)

    axs[2].plot(tecco, integrate(0, Topadvection[expname]), label = plot_labels_list[i], 
    c = exp_colors[expname], linewidth = lw)
    trends_dict[expname][2] =  100 * trend(integrate(0, Topadvection[expname]))
    println("WU:", trends_dict[expname][2])
    axs[1].plot(tecco, integrate(0, Botadvection[expname]), label = plot_labels_list[i], c = exp_colors[expname], linewidth = lw)
    trends_dict[expname][1] =  100 * trend(integrate(0, Botadvection[expname]))
    println("WB:", trends_dict[expname][1])

    axs[3].plot(tecco, integrate(0, Hadvection[expname]), label = plot_labels_list[i], c = exp_colors[expname], linewidth = lw)
    trends_dict[expname][3] =  100 * trend(integrate(0, Hadvection[expname]))
    println("VS:", trends_dict[expname][3])

end
axs[1].tick_params(which="both", bottom=true, left = true)
axs[2].tick_params(which="both", bottom=true)
axs[3].tick_params(which="both", bottom=true)

fig_labs = uppercase.(["a", "b", "c", "d", "e"])
for (i, a) in enumerate(axs)
    a.annotate(fig_labs[i], (0.05, 0.05), fontsize = 17.5, 
    xycoords="axes fraction", fontweight = "bold")
end
[a.grid(alpha = 0.3) for a in axs]
fig.subplots_adjust(wspace = 0.05)
axs[2].legend(loc = "lower center", frameon = false, bbox_to_anchor = (0.5, -0.33), ncols = 3)
fig.savefig(plotsdir("wind_rewrite/1.3.HeatBudgetTimeSeries_external_decomp_faces_all.png"), bbox_inches = "tight", dpi = 1000)
fig