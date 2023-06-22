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
Approxθ_Dict = jldopen(datadir("HeatBudgetApprox" * "_" * region * "_" * suffix * ".jld2"))["Advective_Heat_Budgets"]

fig, axs = plt.subplots(3, 2, figsize =(20, 15), sharex = "col");

E = ones(length(tecco), 2); E[:, 1] .= tecco
get_trend(y) = E * pinv(E) * y
detrend(y) = y .- get_trend(y)
expname = "nosfcadjust"; i = 3

# axs[1, 1].set_title("Meridional Transport")
labels_L = [L"W^{\Delta F, \Delta T}", L"W^{\Delta F}",L"W^{\Delta T}", L"W^{0}"]
vars = ["Vtrsp", "Wbot", "Wtop"]
titles = ["Meridional Transport", "Vertical Transport (3km)", "Vertical Transport (2km)"]

for (j, var) in enumerate(vars)
        ApproxW = Float32.(Approxθ_Dict[expname][var])
        axs[j, 1].plot(tecco, 1e-6 .* ApproxW , color = colors[i])
        axs[j, 1].set_ylabel("Sv")
        β = round((pinv(E) * ApproxW)[1] .* 1e-6, sigdigits = 1)
        noise_std = round(1e-6 .* std(get_trend(ApproxW) .- ApproxW), sigdigits = 2)
        axs[j, 1].plot(tecco, 1e-6 .* get_trend(ApproxW), color = "k", label = string(β) * " [Sv/yr]")
        axs[j, 1].plot(tecco[1], NaN, label = "error = " * string(noise_std) * " [Sv]")    
        axs[j, 1].set_ylabel(labels_L[i] * " [Sv]"); axs[j, 1].legend()
        
        T = 1/12
        p = periodogram(detrend(ApproxW), fs = 1/T)
        p_ifreq  = inv.(freq(p)); p_power = power(p)
        axs[j, 2].plot(p_ifreq, p_power, color = colors[i])
        axs[j, 2].scatter(p_ifreq, p_power, color = colors[i], s = 5)
        axs[j, 2].set_xscale("log")
        axs[j, 2].set_ylabel("power")
        axs[j, 1].set_title(titles[j])
        axs[j, 2].set_title("Periodogram (trend removed)")
end
labels_L = [L"V^{\Delta F, \Delta T}", L"V^{\Delta F}",L"V^{\Delta T}", L"V^{0}"]
axs[1, 1].set_ylabel(labels_L[i] * " [Sv]")

axs[3, 1].set_xlabel("Time (years)")
axs[3, 2].set_xlabel("Period (years)")
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/Reconstructions/" * "AllTrsps" * "_" * expname * "_" *region * "_" * suffix * ".png",
dpi = 1000)

