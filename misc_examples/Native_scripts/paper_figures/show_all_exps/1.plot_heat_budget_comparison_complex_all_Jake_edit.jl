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
gridspec = gs.GridSpec
#blue, red, #green, orangeDataFrames

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)
region = "NPAC"; 

uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= z[:].<= uplvl)
tecco = 1992+1/24:1/12:2018

include(srcdir("plot_and_dir_config.jl"))
sns.set_theme(context = "notebook", style = "ticks",
              palette = colors, rc = custom_params);

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
function integration(t, x)
    int_x = 3.154e+7* cumul_integrate(tecco, x)
    int_x .-= int_x[1]
    return 100 * int_x
end


fig, axs = plt.subplots(1, 2, figsize = (10, 6), sharey = true)
axs[1].set_title("Mid-depth North Pacific\nTemperature Anomaly\n" * L"\theta'")
axs[2].set_title("Heat Budget Term Contributions")
[ax.set_xlabel("time") for ax in axs]
# fig.tight_layout()


lw = 2.5; α = 0.8

exps =  ["iter0_bulkformula", "iter129_bulkformula", "only_buoyancy", 
"only_init", "only_kappa", "only_wind"]
plot_labels_list = ["Iteration 0", "Iteration 129", "Buoyancy Forcing", "Initial Condition", "Mixing Parameters", "Wind Stress"]
plot_labels["Difference"] = "Difference"
exp_colors["Difference"] = "k"
fig.tight_layout()
fig
E, F = trend_matrices(tecco)
trend(x) = (F * x[:])[2]
trends_dict = Dict()
for (i, expname) in enumerate(exps)
    trends_dict[expname] = zeros(5)
    println(expname)

    axs[1].plot(tecco, 100 .* (temps[expname] .- temps[expname][1]), label = plot_labels_list[i], c = exp_colors[expname], linewidth = lw)



    axs[2].plot(tecco, integration(0, int_advection[expname]), label = plot_labels_list[i], c = exp_colors[expname], linewidth = lw)
    trends_dict[expname][2] =  trend(integration(0, int_advection[expname]))

    axs[2].plot(tecco, integration(0, diffusion[expname]), label = plot_labels_list[i], c = exp_colors[expname])
    trends_dict[expname][3] =  trend(integration(0, diffusion[expname]))

    axs[2].plot(tecco, integration(0, GTF[expname]), label = plot_labels_list[i], c = exp_colors[expname], linewidth = lw)
    trends_dict[expname][4] =  trend(integration(0, GTF[expname]))

    axs[2].plot(tecco, integration(0, ext_advection[expname]), label = plot_labels_list[i], c = exp_colors[expname], linewidth = lw, zorder = 10)
    trends_dict[expname][1] =  trend(integration(0, ext_advection[expname]))
    
    trends_dict[expname][5] =  100 * trend(temps[expname])

end

fig_labs = uppercase.(["a", "b", "c", "d", "e"])
for (i, a) in enumerate(axs)
    a.annotate(fig_labs[i], (0.05, 0.05), fontsize = 17, 
    xycoords="axes fraction", fontweight = "bold")
end

axs[2].annotate("Advection\n" * L"\mathbf{A}^{}", (0.5, 0.15), fontsize = 12, 
xycoords="axes fraction", ha="center", fontweight = "normal")
axs[2].annotate("Diffusion\n"* L"\mathbf{F}_\kappa", (0.5, 0.85), fontsize = 12, 
xycoords="axes fraction", ha="center", fontweight = "normal")
axs[2].annotate("Residual\n" * L"\mathbf{R}", (0.8, 0.57), fontsize = 12, 
xycoords="axes fraction", ha="center", fontweight = "normal")
axs[2].annotate("Geothermal\nHeating\n" * L"\mathbf{F}_{geo}", (0.75, 0.73), fontsize = 12, 
xycoords="axes fraction", ha="center", fontweight = "normal")

[a.grid() for a in axs]
fig

trends_df = DataFrame(trends_dict)
trends_df = 1 .* round.(trends_df, digits = 3)
trends_df[!,:variable] = ["External Advection", "Internal Advection", "Diffusion", "Geothermal", "Temperature"]
trends_df

axs[1].legend(loc = "lower left", frameon = false, bbox_to_anchor = (0.4, 0-.3), ncols = 3)
# fig.subplots_adjust(wspace = 0.25)
fig
fig.savefig(plotsdir("native/paper_figures/1.HeatBudgetTimeSeries_complex_all_2panels.png"), bbox_inches = "tight", dpi = 400)
