#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings, DSP, LinearAlgebra
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
expname = "iter0_bulkformula"
Approxθ_Dict = jldopen(datadir("HeatBudgetApprox" * "_" * region * "_" * suffix * ".jld2"))["Advective_Heat_Budgets"]

fig, axsm = plt.subplots(4, 2, figsize =(20, 15), sharex = "col", sharey = "col");

E = ones(length(tecco), 2); E[:, 1] .= tecco
get_trend(y) = E * pinv(E) * y
detrend(y) = y .- get_trend(y)
labels_L = [L"V^{\Delta F, \Delta T}", L"V^{\Delta F}",L"V^{\Delta T}", L"V^{0}"]

@time for (i, expname) in enumerate(keys(shortnames))
    axs = axsm[i, :]
    ApproxV = Float32.(Approxθ_Dict[expname]["Vtrsp"])
    detrend_approxV = ApproxV
    T = 1/12
    p = periodogram(detrend(ApproxV), fs = 1/T)
    p_ifreq  = inv.(freq(p))
    p_power = power(p)
    axs[1].plot(tecco, 1e-6 .* ApproxV , color = colors[i])
    axs[1].set_ylabel("Sv")
    β = round((pinv(E) * ApproxV)[1] .* 1e-6, sigdigits = 1)
    noise_std = round(1e-6 .* std(get_trend(ApproxV) .- ApproxV), sigdigits = 2)
    axs[1].plot(tecco, 1e-6 .* get_trend(ApproxV), color = "k", label = string(β) * " [Sv/yr]")
    axs[1].plot(tecco[1], NaN, label = "error = " * string(noise_std) * " [Sv]")    
    axs[1].set_ylabel(labels_L[i] * " [Sv]"); axs[1].legend()

    axs[2].set_xscale("log")
    axs[2].plot(p_ifreq, p_power, color = colors[i])
    axs[2].scatter(p_ifreq, p_power, color = colors[i], s = 5)
    axs[2].set_ylabel("power")


end

# axsm[4, 1].tick_params(bottom=true, left = true)
# axsm[4, 2].tick_params(bottom=true, left = true)
axsm[1, 1].set_title("Meridional Transport")
axsm[1, 2].set_title("Periodogram (trend removed)")
axsm[4, 1].set_xlabel("Time (years)")
axsm[4, 2].set_xlabel("Period (years)")

# axsm[4, 2].legend(); 
# sns.move_legend(axsm[4, 2], "lower center",  bbox_to_anchor=(.5, -0.75), ncol=3, frameon=true, borderaxespad=0.)
fig
# fig.tight_layout()
# fig.savefig(plotsdir() * "/OHC_Divergence/Reconstructions/" * "MeridionalTrsp" * region * "_" * suffix * ".png",
# dpi = 1000)

