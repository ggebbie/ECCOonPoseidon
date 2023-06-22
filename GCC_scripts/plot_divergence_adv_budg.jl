#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings,
    PyCall, ColorSchemes
using .OHC_helper
import PyPlot as plt
@pyimport seaborn as sns;
@pyimport pandas as pd;

colors =  sns.color_palette("colorblind")[1:4]
labels = ["ECCO", "ECCO Forcing", "ECCO Init. and κ", "CTRL"]

sns.set_theme(context = "poster", style = "ticks", font_scale = 1.2,
              palette = sns.color_palette("colorblind"));

include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());

# abbreviations for each experiment for labels, etc.
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
marks = expsymbols()
nexp = length(shortnames) # number of experiments

tecco = 1992+1/24:1/12:2018
region = "NPAC30_full"; suffix = "2to3"
Trueθ_Dict = jldopen(datadir("HeatBudgetTrue" * "_" * region * "_" * suffix * ".jld2"))["Heat_Budgets"]
Approxθ_Dict = jldopen(datadir("HeatBudgetApprox" * "_" * region * "_" * suffix * ".jld2"))["Advective_Heat_Budgets"]
pygui(false)
fig,axs=plt.subplots(1,3, figsize = (25, 7.5))

@time for (i, expname) in enumerate(keys(shortnames))
    Heat_Budget = Trueθ_Dict[expname]
    geo = jldopen(datadir("Geothermal_"*region*"_2to3.jld2"))["geo"]
    other_terms = Heat_Budget["DiffR"] .+ Heat_Budget["DiffH"] .+ Heat_Budget["AdvR"]
    other_terms = other_terms[1:end-1] .+ geo
    meridional_advection = Heat_Budget["AdvH"][1:end-1]
    θ = Heat_Budget["θ"]
    θ_approx = Float32.(meridional_advection .+ other_terms .+ θ[1])
    axs[1].plot(tecco, θ .- θ[1], color = colors[i]); 
    axs[2].plot(tecco, meridional_advection, color = colors[i], label = labels[i])
    axs[3].plot(tecco, other_terms, color = colors[i]); 
end

axs[1].set_title("North Pacific Temperature \n "*L"\theta - \theta_0", y=1.07)
axs[2].set_title("Meridional Advection \n " * L"\int F_{adv}^Y dt", y=1.05)
axs[3].set_title("Other Terms \n " *L"\int F_{adv}^Z + F_{diff} + F_{geo} dt", y=1.05)
axs[1].set_ylabel("[" *L"^\circ" * " C]")
axs[2].set_ylabel("")
axs[3].set_ylabel("")

axs[2].legend()
sns.move_legend(axs[2], "lower center",  bbox_to_anchor=(.5, -0.34), ncol=4, frameon=true, borderaxespad=0.)
fig.subplots_adjust(wspace=0.3)
[ax.set_xticks([1995, 2005, 2015],[1995, 2005, 2015], rotation=30, ha="right") for ax in axs]
[ax.tick_params(bottom=true) for ax in axs]
for ax in axs
    ax.tick_params(bottom=true, left=true)
    ax.yaxis.set_ticks_position("both")
    ax.tick_params(axis = "x", direction="inout")
    ax.tick_params(axis = "y", direction="in")
end


fig
# fig.tight_layout()
fig.savefig(plotsdir() * "/AGU_Plots/HeatBudget_" * region * ".pdf",bbox_inches = "tight")
fig