include("../../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall, Statistics
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

vars_names =   ["iter0_bulkformula", "only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]

wθ_resid_top = Dict(); wθ_eul = Dict(); wθ_bol = Dict()
wθ_resid_bot = Dict();
Δwθ_top = Dict()
wΔθ_top = Dict()
ΔwΔθ_top = Dict()

Δwθ_bot = Dict()
wΔθ_bot = Dict()
ΔwΔθ_bot = Dict()

integrate(start, x) = cumsum([start, x...])[1:end-1]
include(srcdir("plot_and_dir_config.jl"))



for (i, expname) in enumerate(["mean_tau_yesadjust_redo_bf", "mean_tau_noadjust_redo_bf"])
    fname = datadir("native/" * expname * region * "_THETA_budget_Wθ_eul_bol_ΔDecomp_climatological" * suffix * ".jld2")
    print(expname)
    vars = jldopen(fname)["dθ"]

    wθ_resid_top[expname] = 100 .* integrate(0, - vars["wθTop"]).* 2.628f+6
    Δwθ_top[expname] = 100 .* integrate(0, - vars["ΔW_θTop"]).* 2.628f+6
    wΔθ_top[expname] = 100 .* integrate(0, - vars["W_ΔθTop"]).* 2.628f+6
    ΔwΔθ_top[expname] = 100 .* integrate(0, - vars["ΔW_ΔθTop"]).* 2.628f+6

    wθ_resid_bot[expname] = 100 .* integrate(0, vars["wθBot"]).* 2.628f+6
    Δwθ_bot[expname] = 100 .* integrate(0, vars["ΔW_θBot"] ).* 2.628f+6
    wΔθ_bot[expname] = 100 .* integrate(0, vars["W_ΔθBot"] ).* 2.628f+6
    ΔwΔθ_bot[expname] = 100 .* integrate(0, vars["ΔW_ΔθBot"] ).* 2.628f+6
end

exp_colors["mean_tau_noadjust_redo_bf"] = exp_colors["iter0_bulkformula"] 
exp_colors["mean_tau_yesadjust_redo_bf"] = exp_colors["only_wind"] 

sns.set_style("darkgrid", Dict("axes.facecolor" => ".95", "grid.color" => "0.5"))

vars =   [ "only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]
E,F = trend_matrices(tecco)
get_trend(x) = (F * x)[2]
var_explained(y, x) = 100 * ( 1 - (var(y - x) / var(y))) 



wθ_resid_list = [wθ_resid_bot, wθ_resid_top]
Δwθ_list = [Δwθ_bot, Δwθ_top]
wΔθ_list = [wΔθ_bot, wΔθ_top]
ΔwΔθ_list = [ΔwΔθ_bot, ΔwΔθ_top]

lw = 2.5

fig, ax = plt.subplots(1, 2, figsize = (12, 6), sharey = true)
ax[1].set_title(L"\mathbf{A}^{\mathcal{B}}" * " Decomposition")
ax[2].set_title(L"\mathbf{A}^{\mathcal{T}}"  * " Decomposition")
var_name = [L"^\mathcal{B}", L"^\mathcal{U}"]
fig
colors_alt = Dict(
    "Δwθ" => "#c5b0d5",  
    "wΔθ" => "#ffbb78", 
    "ΔwΔθ" => "#98df8a", 
    "mean_tau_noadjust_redo_bf" => exp_colors["mean_tau_noadjust_redo_bf"], 
    "mean_tau_yesadjust_redo_bf" => exp_colors["mean_tau_yesadjust_redo_bf"]
)
ax[1].grid(alpha = 0.3)
ax[2].grid(alpha = 0.3)
for (i, axs) in enumerate(ax)
    wθ_resid = wθ_resid_list[i]
    Δwθ = Δwθ_list[i]
    wΔθ = wΔθ_list[i]
    ΔwΔθ = ΔwΔθ_list[i]
    axs.set_xlabel("time", fontweight = "bold")

    axs.plot(tecco, wθ_resid["mean_tau_noadjust_redo_bf"], label =  "Seasonal First-Guess Wind, " * L"\mathbf{A}_{\overline{0}}", 
    linewidth = 2.5, color = colors_alt["mean_tau_noadjust_redo_bf"], alpha = 0.8)

    axs.plot(tecco,  wθ_resid["mean_tau_yesadjust_redo_bf"], label =  "Seasonal Adjusted Wind, " * L"\mathbf{A}_{\overline{129}}", 
    linewidth = 2.5, color = colors_alt["mean_tau_yesadjust_redo_bf"], alpha = 0.8)

    axs.plot(tecco, Δwθ["mean_tau_yesadjust_redo_bf"] , label =  "Boundary Velocity Response,\n" * L"\mathbf{A}({\Delta w^{res}, \theta_\overline{0}})", 
    linewidth = 2.5, colors_alt["Δwθ"], alpha = 1.0)

    axs.plot(tecco, wΔθ["mean_tau_yesadjust_redo_bf"]  , label =  "Boundary Temperature Response,\n" * L"\mathbf{A}({ w^{res}_\overline{0}}, \Delta \theta)", 
    linewidth = lw,  colors_alt["wΔθ"], alpha = 0.95)

    axs.plot(tecco, ΔwΔθ["mean_tau_yesadjust_redo_bf"]  , label = "Boundary Response Interaction,\n" * L"\mathbf{A}({\Delta w^{res}, \Delta \theta})", 
    linewidth = lw,  colors_alt["ΔwΔθ"], alpha = 0.5, linestyle = "--")
    # axs.axhline(0, c = "k", linestyle = "--", alpha = 0.5)
end

handles, labels = ax[1].get_legend_handles_labels()

# First legend (2 items)
ax[1].legend(handles[1:2], labels[1:2], loc="upper center", bbox_to_anchor=(0.5, 0-.275), ncol=1, frameon=false)
ax[2].legend(handles[3:end], labels[3:end], loc="upper center", bbox_to_anchor=(0.6, 0-.25), ncol=1, frameon=false)

# ax[1].legend(loc = "lower left", frameon = false, bbox_to_anchor = (0.05, 0-.7), ncols = 2)


ax[1].set_ylabel("[cK]", fontweight = "bold")
fig_labs = uppercase.(["a", "b", "c", "d", "e", "f"])
for (i, a) in enumerate(ax)
    a.annotate(fig_labs[i], (0.05, 0.02), fontsize = 20, color = "black", 
    xycoords="axes fraction", fontweight = "bold")
end

ax[1].tick_params(which="both", bottom=true, left = true)
ax[2].tick_params(which="both", bottom=true)
fig.subplots_adjust(wspace = 0.05)

fig.savefig(plotsdir("wind_rewrite/2.1.HeatBudgetTimeSeries_vertical_heatbudget_mathlabels.png"), bbox_inches = "tight", dpi = 400)
fig