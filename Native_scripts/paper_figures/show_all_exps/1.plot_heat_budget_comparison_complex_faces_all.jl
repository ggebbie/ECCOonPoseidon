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

#blue, red, #green, orangeDataFrames

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)
region = "NPAC"; 

uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= z[:].<= uplvl)
tecco = 1992+1/24:1/12:2018

include(srcdir("plot_and_dir_config.jl"))
sns.set_theme(context = "paper", style = "ticks",
              palette = colors, rc = custom_params);
exps =  ["iter0_bulkformula", "iter129_bulkformula", "only_buoyancy", 
              "only_init", "only_kappa", "only_wind"]
Hadvection = Dict(); diffusion = Dict(); GTF = Dict()
Vadvection = Dict()
Topadvection = Dict()
Botadvection = Dict()

temps = Dict()
for (i, expname) in enumerate(exps)
    fname = datadir("native/" * expname * region * "_THETA_budget_ref_with_Bolus" * suffix * ".jld2")
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

temps["Difference"] = temps["iter129_bulkformula"] .- temps["iter0_bulkformula"]

Hadvection["Difference"] = Hadvection["iter129_bulkformula"] .- Hadvection["iter0_bulkformula"]
Vadvection["Difference"] = Vadvection["iter129_bulkformula"] .- Vadvection["iter0_bulkformula"]

diffusion["Difference"] = diffusion["iter129_bulkformula"] .- diffusion["iter0_bulkformula"]
GTF["Difference"] = GTF["iter129_bulkformula"] .- GTF["iter0_bulkformula"]

fig, axs = plt.subplots(1, 3, figsize = (8, 5), sharex = false, sharey = true)
axs[3].set_title("Southern Boundary \n Advective Heat Flux \n" * L"\mathbf{A}_{\mathcal{S}}")
axs[2].set_title("Upper Boundary \n Advective Heat Flux  \n" * L"\mathbf{A}_{\mathcal{U}}")
axs[1].set_title("Bottom Boundary \n Advective Heat Flux  \n " * L"\mathbf{A}_{\mathcal{B}}")

axs[1].set_ylabel("[cK]", fontweight = "bold")
[ax.set_xlabel("time", fontweight = "bold") for ax in axs]
fig.tight_layout()

integrate(start, x) = cumsum([start, x...])[1:end-1] .* (100 * 2.628e+6)

lw = 1.5; α = 0.8

exps =  ["iter0_bulkformula", "iter129_bulkformula", "only_buoyancy", 
"only_init", "only_kappa", "only_wind"]
plot_labels["Difference"] = "Difference"
exp_colors["Difference"] = "k"
plot_labels_list = ["Iteration 0", "Iteration 129", "Buoyancy Forcing", "Initial Condition", "Mixing Parameters", "Wind Stress"]

E, F = trend_matrices(tecco)
trend(x) = (F * x[:])[2]
trends_dict = Dict()
for (i, expname) in enumerate(exps)
    trends_dict[expname] = zeros(3)
    println(expname)
    axs[2].plot(tecco, integrate(0, Topadvection[expname]), label = plot_labels_list[i], c = exp_colors[expname], linewidth = lw)
    trends_dict[expname][2] =  trend(integrate(0, Topadvection[expname]))
    axs[1].plot(tecco, integrate(0, Botadvection[expname]), label = plot_labels_list[i], c = exp_colors[expname], linewidth = lw)
    trends_dict[expname][1] =  trend(integrate(0, Botadvection[expname]))

    axs[3].plot(tecco, integrate(0, Hadvection[expname]), label = plot_labels_list[i], c = exp_colors[expname], linewidth = lw)
    trends_dict[expname][3] =  trend(integrate(0, Hadvection[expname]))


end
fig_labs = uppercase.(["a", "b", "c", "d", "e"])
for (i, a) in enumerate(axs)
    a.annotate(fig_labs[i], (0.05, 0.05), fontsize = 15, 
    xycoords="axes fraction", fontweight = "bold")
end
trend(temps["iter0_bulkformula"])
trends_df = DataFrame(trends_dict)
trends_df = 1 .* round.(trends_df, digits = 3)
# trends_df.Difference .= trends_df.iter129_bulkformula - trends_df.iter0_bulkformula 
trends_df[!,:variable] = ["Bottom Advection", "Top Advection", "Horizontal Advection"]
trends_df
[a.grid() for a in axs]
axs[2].legend(loc = "lower center", frameon = false, bbox_to_anchor = (0.5, -0.35), ncols = 3)
fig.tight_layout()
fig.subplots_adjust(wspace = 0.1)
fig
fig.savefig(plotsdir("native/paper_figures/2.HeatBudgetTimeSeries_external_decomp_faces_all.png"), bbox_inches = "tight", dpi = 1000)
