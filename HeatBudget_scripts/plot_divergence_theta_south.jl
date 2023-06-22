#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings, DSP 
using .OHC_helper
import PyPlot as plt
using PyCall
using ColorSchemes
@pyimport seaborn as sns;
@pyimport pandas as pd;

colors =  sns.color_palette("deep")[1:4]
sns.set_theme(context = "talk", font_scale = 1.0,
              palette = sns.color_palette("deep"));

# tyle("ticks")
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
region = "NPAC"; suffix = "2to3"
expname = "iter0_bulkformula"
Approxθ_Dict = jldopen(datadir("HeatBudgetApprox" * "_" * region * "_" * suffix * ".jld2"))["Advective_Heat_Budgets"]

fig, axsm = plt.subplots(4, 1, figsize =(10, 10), sharex = "col", sharey = "col");
axsm[4, 1].tick_params(bottom=true, left=true)

E = ones(length(tecco), 2); E[:, 1] .= tecco
get_trend(y) = E * pinv(E) * y
detrend(y) = y .- get_trend(y)
pygui(true)
@time for (i, expname) in enumerate(keys(shortnames))
    axs = axsm[i]
    Approxθ = Float32.(Approxθ_Dict[expname]["θSouth"])

    axs.plot(tecco, Approxθ , color = colors[i])
    axs.set_ylabel("Sv")
    β = round((pinv(E) * Approxθ)[1], sigdigits = 3)
    # axs.plot(tecco, get_trend(Approxθ), color = "k", label = string(β) * L"^\circ C/yr")
    axs.set_ylabel(L"^\circ C"); axs.legend()

    # pygui(true)
end


axsm[1].set_title("Average Temperature on south face")
# axsm[1, 2].set_title("Periodogram (trend removed)")
axsm[4].set_xlabel("Time (years)")
# axsm[4, 2].set_xlabel("Period (years)")

# axsm[4, 2].legend(); 
# sns.move_legend(axsm[4, 2], "lower center",  bbox_to_anchor=(.5, -0.75), ncol=3, frameon=true, borderaxespad=0.)

fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/Reconstructions/" * "Southθ" * region * "_" * suffix * ".png",
dpi = 1000)

