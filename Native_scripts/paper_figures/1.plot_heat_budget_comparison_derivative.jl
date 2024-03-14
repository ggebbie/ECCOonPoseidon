include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour, Statistics,
    MeshArrays, MITgcmTools, JLD2, DrWatson,
    BenchmarkTools, LaTeXStrings, PyCall, DataFrames
import NaNMath as nm
import PyPlot as plt

import NumericalIntegration
include(srcdir("config_exp.jl"))
0
cmo = pyimport("cmocean.cm");

#blue, red, #green, orangeDataFrames

# setup some data
x = collect(-π/2 : π/1000 : π/2)
y = sin.(x)
dx = x[2:end] .- x[1:end-1]

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

exps =  ["iter0_bulkformula", "iter129_bulkformula"]
Hadvection = Dict(); diffusion = Dict(); GTF = Dict()
Vadvection = Dict();
tot_advection = Dict()
temps = Dict()
for (i, expname) in enumerate(exps)
    fname = datadir("native/" * expname * region * "_THETA_budget_ref_with_Bolus_wextra" * suffix * ".jld2")
    vars = load(fname)["dθ"]
    temps[expname] = vars["θ"]
    uθ = vars["VθSouth"]
    ∇wθ = vars["wθBot"] .- vars["wθTop"]
    Hadvection[expname] = (uθ)
    Vadvection[expname] = (∇wθ)
    tot_advection[expname] = vars["wθ"] .+ vars["uvθ"]
    diffusion[expname] = (vars["κxyθ"] .+ vars["κzθ"])
    GTF[expname] = (vars["GTH"])
end
fig

diffusion["Difference"] = diffusion["iter129_bulkformula"] .- diffusion["iter0_bulkformula"]
GTF["Difference"] = GTF["iter129_bulkformula"] .- GTF["iter0_bulkformula"]

lw = 2.5; α = 0.8
crude_integrate(start, x) = cumsum([start, x...])[1:end-1] .* (100 * 2.628e+6)

exps =  ["iter0_bulkformula", "iter129_bulkformula", "only_init", "only_wind", "only_buoyancy", "only_kappa"]
E, F = trend_matrices(tecco)
trend(x) = (F * x[:])[2]
trends_dict = Dict()

fig, ax = plt.subplots()

for (i, expname) in enumerate(["iter0_bulkformula", "iter129_bulkformula"])
    println(expname)

    approx = Vadvection[expname] .+ Hadvection[expname] 
    actual = tot_advection[expname]
    resid = actual .- approx
    ax.plot(approx, label = "approx 1")
    ax.plot(actual, label = "actual")
    ax.plot(resid, label = "resid")

    # println((actual))
    println(100 * 3.154e+7 * mean(resid))
    println(100 * 3.154e+7 * mean(Vadvection[expname]))
    println(100 * 3.154e+7 * mean(Hadvection[expname]))

    println(100 * 3.154e+7 * mean(actual))
    println(mean(resid) / mean(actual))

end
covar
fig
t = tecco
for (i, expname) in enumerate(["only_kappa"])
    println(expname)

    approx = cumul_integrate(t, Vadvection[expname] .* 3.154e+7) .+ cumul_integrate(t, Hadvection[expname] .* 3.154e+7) .+ 
    cumul_integrate(t, diffusion[expname] .* 3.154e+7) .+ cumul_integrate(t, GTF[expname] .* 3.154e+7)
    approx .-= approx[1]; approx= 100 .* approx
    actual = 100 .* (temps[expname] .- temps[expname][1])
    resid = actual .- approx
    ax.plot(approx, label = "approx 2")

    # println(actual[2:5] .- approx[2:5])
    # println((resid))
    # println((actual))
    # println(std(resid))
    # println(std(actual))
    println(std(resid) / std(actual))

end
ax.legend()

fig

cumul_integrate(t, Vadvection["iter0_bulkformula"] )
cumul_integrate(t, Vadvection["iter0_bulkformula"]) 

fig
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
