include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson,
    BenchmarkTools, LaTeXStrings, PyCall, DataFrames
import NaNMath as nm
import PyPlot as plt

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

exps =  ["iter129_bulkformula", "iter0_bulkformula"]
Hadvection = Dict(); diffusion = Dict(); GTF = Dict()
Vadvection = Dict()
temps = Dict()
for (i, expname) in enumerate(exps)
    fname = datadir("native/" * expname * region * "_THETA_budget_ref_with_Bolus" * suffix * ".jld2")
    vars = load(fname)["dθ"]
    temps[expname] = vars["θ"]
    uθ = vars["VθSouth"]
    ∇wθ = vars["wθBot"] .- vars["wθTop"]
    Hadvection[expname] = (uθ)
    Vadvection[expname] = (∇wθ)

    diffusion[expname] = (vars["κxyθ"] .+ vars["κzθ"])
    GTF[expname] = (vars["GTH"])
end

temps["Difference"] = temps["iter129_bulkformula"] .- temps["iter0_bulkformula"]

Hadvection["Difference"] = Hadvection["iter129_bulkformula"] .- Hadvection["iter0_bulkformula"]
Vadvection["Difference"] = Vadvection["iter129_bulkformula"] .- Vadvection["iter0_bulkformula"]

diffusion["Difference"] = diffusion["iter129_bulkformula"] .- diffusion["iter0_bulkformula"]
GTF["Difference"] = GTF["iter129_bulkformula"] .- GTF["iter0_bulkformula"]

fig, axs = plt.subplots(2, 2, figsize = (6, 6), sharex = false, sharey = true)
axs[1].set_title(L"\mathbf{A}_Z")
axs[2].set_title(L"\mathbf{A}_H")
axs[3].set_title(L"\mathbf{F}_\kappa")
axs[4].set_title(L"\mathbf{F}_{geo}")
axs[1].set_ylabel("[cK]", fontweight = "bold")
axs[2].set_ylabel("[cK]", fontweight = "bold")
[ax.set_xlabel("time", fontweight = "bold") for ax in axs]
fig.tight_layout()

integrate(start, x) = cumsum([start, x...])[1:end-1] .* (100 * 2.628e+6)

lw = 2.5; α = 0.8

exps =  ["iter0_bulkformula", "iter129_bulkformula", "Difference"]
plot_labels["Difference"] = "Effect of All Control Adjustments"
exp_colors["Difference"] = "k"

E, F = trend_matrices(tecco)
trend(x) = (F * x[:])[2]
trends_dict = Dict()
for (i, expname) in enumerate(exps)
    trends_dict[expname] = zeros(5)
    println(expname)
    axs[1].plot(tecco, integrate(0, Vadvection[expname]), label = plot_labels[expname], c = exp_colors[expname], linewidth = lw)
    trends_dict[expname][1] =  trend(integrate(0, Vadvection[expname]))
    println("Vertical")
    println(round(integrate(0, Vadvection[expname])[end], digits = 2))
    axs[2].plot(tecco, integrate(0, Hadvection[expname]), label = plot_labels[expname], c = exp_colors[expname], linewidth = lw)
    trends_dict[expname][2] =  trend(integrate(0, Hadvection[expname]))

    println("Horizontal")
    println(round(integrate(0, Hadvection[expname])[end], digits = 2))
    axs[3].plot(tecco, integrate(0, diffusion[expname]), label = plot_labels[expname], c = exp_colors[expname])
    trends_dict[expname][3] =  trend(integrate(0, diffusion[expname]))

    println("Diff")
    println(round(integrate(0, diffusion[expname])[end], digits = 2))
    axs[4].plot(tecco, integrate(0, GTF[expname]), label = plot_labels[expname], c = exp_colors[expname], linewidth = lw)
    trends_dict[expname][4] =  trend(integrate(0, GTF[expname]))

    println("Geo")
    println(round(integrate(0, GTF[expname])[end], digits = 2))
    trends_dict[expname][5] =  100 * trend(temps[expname])

end
trend(temps["iter0_bulkformula"])
trends_df = DataFrame(trends_dict)
trends_df = 1 .* round.(trends_df, digits = 3)
trends_df.Difference .= trends_df.iter129_bulkformula - trends_df.iter0_bulkformula 
trends_df
[a.grid() for a in axs]
axs[2].legend(loc = "lower left", frameon = false, bbox_to_anchor = (0.1, -0.4), ncols = 5)
fig.subplots_adjust(wspace = 0.1)
fig
fig.savefig(plotsdir("native/paper_figures/1.HeatBudgetTimeSeries_simple.png"), bbox_inches = "tight", dpi = 500)
