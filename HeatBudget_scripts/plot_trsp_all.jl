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
E = ones(length(tecco), 2); E[:, 2] .= tecco 
OLS(E, y, W) = inv(E'*inv(W)*E)*E'*inv(W)*y
W = ones(length(tecco));# W[1:12*5] .= 5; W[(12*5)+1:end-(12*5)] .= 50; W[end-(12*5)+1:end] .= 5
get_trend(E, y) = E * OLS(E, y .- mean(y), Diagonal(W)) .+ mean(y)
detrend(x) = x .- get_trend(E, x)
region = "NPAC30_full"; suffix = "2to3"
Approxθ_Dict = jldopen(datadir("HeatBudgetApprox" * "_" * region * "_" * suffix * ".jld2"))["Advective_Heat_Budgets"]
fig, axs = plt.subplots(4, 4, figsize =(20, 15), sharex = "col", sharey = "row");

# axs[1, 1].set_title("Meridional Transport")
labels_L = [L"U^{\Delta F, \Delta T}", L"U^{\Delta F}",L"U^{\Delta T}", L"U^{0}"]
vars = ["Vtrsp", "Wbot", "Wtop"]
titles = [L"V_{in}", L"W_{in}" *" (3km)", L"W_{out}"*"(2km)"]
mult = [1]
for (i, expname) in enumerate(keys(shortnames))
    axs[1, i].set_title(labels_L[i])
    sum = 0.0 .* Approxθ_Dict[expname][vars[1]]
    for (j, var) in enumerate(vars)
            (j < 3) ? mult[1] = 1 : mult[1] = -1
            print(j)
            ApproxVel = mult[1] * 1e-6.*Float32.(Approxθ_Dict[expname][var])
            sum.+=ApproxVel
            axs[j, i].plot(tecco, ApproxVel , color = colors[i])
            a, β = round.(OLS(E, ApproxVel .- mean(ApproxVel), I), sigdigits = 2); 
            noise_std = round(std(get_trend(E, ApproxVel) .- ApproxVel), sigdigits = 2)
            axs[j, i].plot(tecco,get_trend(E, ApproxVel), color = "k", label = string(β) * " [Sv/yr]")
            axs[j, i].plot(tecco[1], NaN, label = "error = " * string(noise_std) * " [Sv]")    
            axs[j, 1].set_ylabel(titles[j] * "\n [Sv]"); axs[j, i].legend()
            
            #plot frequencies
            # T = 1/12
            # p = periodogram(detrend(ApproxVel), fs = 1/T)
            # p_ifreq  = inv.(freq(p)); p_power = power(p)
            # axs[j, 2].plot(p_ifreq, p_power, color = colors[i])
            # axs[j, 2].scatter(p_ifreq, p_power, color = colors[i], s = 5)
            # axs[j, 2].set_xscale("log")
            # axs[j, 2].set_ylabel("power")
            # axs[j, 1].set_title(titles[j])
            # axs[j, 2].set_title("Periodogram (trend removed)")
    end
    axs[4, i].plot(tecco, sum , color = colors[i])

end
fig.subplots_adjust(wspace=0.1, hspace=0.0)
[axs[3, i].set_xlabel("Time (years)") for i=1:4]

# fig.tight_layout()
# fig.savefig(plotsdir() * "/OHC_Divergence/Reconstructions/" * "AllTrsps" * "_" *region * "_" * suffix * ".png",
# dpi = 1000)

