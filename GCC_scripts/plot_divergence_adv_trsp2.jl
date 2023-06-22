#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings, Statistics
using .OHC_helper
import PyPlot as plt
using PyCall
using ColorSchemes
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
E = ones(length(tecco), 2); E[:, 2] .= tecco
F = (E'*E)\E' # least squares estimator

fig,axs=plt.subplots(1,2, figsize = (19 , 8))
# fig
β = [0., 0.]
Vtrsp_mean = zeros(312)
linregress(y) = E * F * y
derivative(x, dt) = (x[2:end] .- x[1:end-1]) ./ dt
using LinearAlgebra
# for (i, expname) in enumerate(keys(shortnames))
i  = 1
expname = "iter129_bulkformula"
# print(i)
# print(expname)
Heat_Budget_true = Trueθ_Dict[expname]
Heat_Budget_approx = Approxθ_Dict[expname]

meridional_advection = Heat_Budget_true["AdvH"][1:end-1]
NPAC_Volume =  2.2557688f16
Vtrsp = Float32.(Heat_Budget_approx["Vtrsp"])
E_θ = (Vtrsp)

F_θ = (E_θ'*E_θ)\E_θ' # least squares estimator
Vθsouth_true = derivative(Heat_Budget_true["AdvH"], 2.628e+6) .* NPAC_Volume
Vtrsp_est = F_θ*(Vθsouth_true)

meridional_advection_approx = repeat([mean(E_θ) .* Vtrsp_est], 312) ./ NPAC_Volume
meridional_advection_approx_mean = cumsum(cat(0.0, meridional_advection_approx .* 2.628e+6, dims = 1))
Vθp = (Vθsouth_true .- (E_θ .* Vtrsp_est))./ NPAC_Volume 
thetap = Vθp .* NPAC_Volume ./ E_θ
meridional_advection_approx_prime = cumsum(cat(0.0, Vθp .* 2.628e+6, dims = 1))
if i == 1 
    cmaps = sns.dark_palette(colors[i], reverse = true)
    # r2 = cor(meridional_advection[1:end], meridional_advection_approx_mean[1:end-1])^2
    axs[1].plot(tecco, meridional_advection_approx_prime[1:end-1], color = cmaps[3], 
    linestyle = ":", label = L" \int  \frac{\partial}{\partial y} (\theta V ') dt" ); 
    axs[1].plot(tecco, meridional_advection_approx_mean[1:end-1], color = cmaps[2], 
    linestyle = "--", label = L" \int  \frac{\partial}{\partial y} ( \theta \bar{V}) dt" ); 
    axs[1].plot(tecco, meridional_advection, color = cmaps[1], label = "SUM"); 
    axs[1].set_ylabel("[" * L"^\circ " * "C]")
end
axs[2].plot(tecco, Vtrsp, color = colors[i], label =labels[i]); 
Vtrsp_mean.+=Vtrsp./4
println(expname)
println(meridional_advection_approx_mean[end-1] / meridional_advection[end])
# axs[1].hlines([0], minimum(tecco), maximum(tecco), color = "grey", alpha = 0.7); 
for ax in axs
    ax.tick_params(bottom=true, left=true)
    ax.yaxis.set_ticks_position("both")
    ax.tick_params(axis = "x", direction="inout")
    ax.tick_params(axis = "y", direction="in")
end
# in ECCO " * L"[F_{adv}^Y = \frac{\partial}{\partial y} (V \theta)]
axs[1].legend(ncol = 1); axs[1].set_title("Decomposition of ECCO \n Meridional Advection")
axs[2].set_ylabel("[Sv]")

# axs[1].set_title("Not sure what to title")
axs[2].set_title("Meridional Volume Transport \n")
axs[2].legend()
β .= F * Vtrsp_mean
axs[2].plot(tecco, E * β, color = "k"); 
# sns.move_legend(axs[1], "lower center",  bbox_to_anchor=(0., 0), ncol=2, frameon=true, borderaxespad=0.)
# sns.move_legend(axs[2], "lower center",  bbox_to_anchor=(.5, -0.3), ncol=2, frameon=true, borderaxespad=0.)
sns.move_legend(axs[2], "upper right",  bbox_to_anchor=(1.8, 1), ncol=1, frameon=true, borderaxespad=0.)

fig.subplots_adjust(wspace=0.5)
fig
# fig.savefig(plotsdir() * "/AGU_Plots/TransportS_" * region * "_" * suffix * ".pdf",bbox_inches = "tight")


