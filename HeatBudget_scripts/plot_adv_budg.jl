#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings
using .OHC_helper
import PyPlot as plt
using PyCall
using ColorSchemes
@pyimport seaborn as sns;
@pyimport pandas as pd;

colors =  sns.color_palette("deep")[1:4]
sns.set_theme(context = "talk", font_scale = 1.0,
              palette = sns.color_palette("deep"));

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
region = "NPAC30"; suffix = "2to3"
Trueθ_Dict = jldopen(datadir("HeatBudgetTrue" * "_" * region * "_" * suffix * ".jld2"))["Advective_Heat_Budgets"]
Approxθ_Dict = jldopen(datadir("HeatBudgetApprox" * "_" * region * "_" * suffix * ".jld2"))["Advective_Heat_Budgets"]

pygui(false)
fig, axs = plt.subplots(1, 3, figsize =(15, 7.5), sharey = "col", sharex = "col");
vθ_sym = L"\int \frac{\partial}{\partial y}(v \theta) dt "
wθ_sym = L"\int \frac{\partial}{\partial z}(w \theta) dt"
vθwθ_sym = L"\int\frac{\partial}{\partial y}(v\theta)+\frac{\partial }{\partial z}(w\theta)dt"
axs[1].set_title("Lateral Advection \n" * vθ_sym)
axs[2].set_title("Vertical Advection \n" *wθ_sym)
axs[3].set_title("Sum \n" * vθ_sym * L"+" * wθ_sym)
labels_L = [L"\theta^{\Delta F, \Delta T}", L"\theta^{\Delta F}",L"\theta^{\Delta T}", L"\theta^{0}"]
expname = "iter129_bulkformula"; i = 1

Trueθ = Trueθ_Dict[expname]
Approxθ = Approxθ_Dict[expname]
axs[1].plot(tecco, Trueθ["AdvH"][1:end-1], color = colors[i])
axs[1].plot(tecco, Approxθ["AdvH"][1:end-1], linestyle = "--", color = colors[i], alpha = 0.9)
axs[1].plot(tecco, Approxθ["AdvHmean"][1:end-1], linestyle = "dotted", color = colors[i], linewidth=3, alpha = 0.8)
axs[1].set_ylabel( labels_L[i] * " " * L"[^\circ C]")

axs[2].plot(tecco, Trueθ["AdvR"][1:end-1], label = "truth", color = colors[i])
axs[2].plot(tecco, Approxθ["AdvR"][1:end-1], label = "FV", linestyle = "--", color = colors[i], alpha = 0.9)
axs[2].plot(tecco, Approxθ["AdvRmean"][1:end-1], label = "FV (constant θ)", linestyle = "dotted", color = colors[i], 
linewidth=3, alpha = 0.8)

axs[3].plot(tecco, (Trueθ["AdvH"] .+ Trueθ["AdvR"])[1:end-1], color = colors[i])
axs[3].plot(tecco, (Approxθ["AdvH"] .+ Approxθ["AdvR"])[1:end-1], linestyle = "--", color = colors[i], alpha = 0.9)
axs[3].plot(tecco, (Approxθ["AdvHmean"] .+ Approxθ["AdvRmean"])[1:end-1], linestyle = "dotted", 
color = colors[i], linewidth=3, alpha = 0.8)
pad = 5 # in points


[ax.set_xlabel("time") for ax in axs]
axs[2].legend(); 
sns.move_legend(axs[2], "lower center",  bbox_to_anchor=(.5, -0.3), ncol=3, frameon=true, borderaxespad=0.)

fig
# fig.tight_layout()
# fig.savefig(plotsdir() * "/OHC_Divergence/Reconstructions/" * "AdvHeadBudget" * "_" * expname * "_" * region * "_" * suffix * ".png",
# dpi = 1000)



t = derivative(Heat_Budget_approx["AdvH"], 2.628e+6) .\ (Heat_Budget_approx["Vtrsp"]  .* Heat_Budget_approx["θSouth"])

mean(t)
plot( )
