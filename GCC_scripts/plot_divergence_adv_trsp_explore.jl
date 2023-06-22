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
sns.set_theme(context = "poster", font_scale = 1.2,
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
Trueθ2_Dict = jldopen(datadir("HeatBudgetTrue" * "_" * region * "_" * suffix * "_2.jld2"))["Heat_Budgets"]
Approxθ_Dict = jldopen(datadir("HeatBudgetApprox" * "_" * region * "_" * suffix * ".jld2"))["Advective_Heat_Budgets"]

expname = "iter129_bulkformula"

Heat_Budget_true2 = Trueθ2_Dict[expname]

Heat_Budget_true = Trueθ_Dict[expname]
Heat_Budget_approx = Approxθ_Dict[expname]
# Heat_Budget_approx2 = Approxθ_Dict2[expname]
using Plots

plot(Heat_Budget_true["AdvH"])
plot!(Heat_Budget_approx["Vtrsp"])

# meridional_advection = Heat_Budget_true["AdvH"][1:end-1]
# E_trsp = Heat_Budget_approx["Vtrsp"]
# F_trsp = inv(E_trsp'E_trsp) * E_trsp'
# x = F_trsp * meridional_advection
# Vtrsp = 1e-6.*Heat_Budget_approx["Vtrsp"] #convert to SV 
derivative(x, dt) = (x[2:end] .- x[1:end-1]) ./ dt

#estimate Θ
Vθsouth_true = derivative(Heat_Budget_true["AdvH"], 2.628e+6) .*  Heat_Budget_approx["Volume"]
Vθsouth_true2 = derivative(Heat_Budget_true2["DiffH"], 2.628e+6) .*  Heat_Budget_approx["Volume"]

using LinearAlgebra
E = Heat_Budget_approx["Vtrsp"]
F = (E'*E)\E' # least squares estimator
θ_est = F*Vθsouth_true
θ_est2 = F*Vθsouth_true2

#estimate Θ
Vθsouth_est = (Heat_Budget_approx["Vtrsp"] .* θ_est)

fig,axs=plt.subplots(1,1, figsize = (25, 7.5))
axs.plot(Vθsouth_true, label = "true")
axs.plot(Vθsouth_est, label = "estimate")
axs.plot(Vθsouth_est, label = "estimate")

axs.legend()
fig

#estimate Θ
Vθsouth_true = derivative(Heat_Budget_true["AdvH"], 2.628e+6) .*  Heat_Budget_approx["Volume"]
Vθsouth_est = (Heat_Budget_approx["Vtrsp"] .* θ_est2)

fig,axs=plt.subplots(1,1, figsize = (25, 7.5))
axs.plot(Vθsouth_true, label = "true")
axs.plot(Vθsouth_est, label = "estimate")
axs.legend()
fig


fig,axs=plt.subplots(1,1, figsize = (25, 7.5))
axs.plot(Heat_Budget_approx["θSouth"], label = "true")
axs.plot(θ_est, label = "estimate")
axs.legend()
axs.set_ylim(1.54, 1.56)
fig


fig,axs=plt.subplots(1,1, figsize = (25, 7.5))
axs.plot(Heat_Budget_approx["AdvHmean"], label = "true")
axs.plot(Heat_Budget_true["AdvH"], label = "estimate")
axs.legend()
fig

fig,axs=plt.subplots(1,1, figsize = (25, 7.5))
axs.plot(Heat_Budget_approx["Vtrsp"], label = "true")
axs.legend()
fig

fig,axs=plt.subplots(1,1, figsize = (25, 7.5))
axs.plot(Vθsouth_true .- Vθsouth_est, label = "true")
axs.legend()
fig

fig,axs=plt.subplots(1,1, figsize = (25, 7.5))
axs.plot(Heat_Budget_true["AdvH"] .+ Heat_Budget_true["AdvR"] 
.+ Heat_Budget_true["DiffH"] .+ Heat_Budget_true["DiffR"] .+ Heat_Budget_true["θ"][1], label = "estimate")
axs.plot(Heat_Budget_true["θ"], label = "estimate")
geo_therm = Heat_Budget_true["AdvH"][1:end-1] .+ Heat_Budget_true["AdvR"][1:end-1] .+ Heat_Budget_true["DiffH"][1:end-1] .+ Heat_Budget_true["DiffR"][1:end-1]  .- Heat_Budget_true["θ"] .+ Heat_Budget_true["θ"][1]
axs.legend()
fig

fig,axs=plt.subplots(1,1, figsize = (25, 7.5))
geo_therm = Heat_Budget_true["AdvH"][1:end-1] .+ Heat_Budget_true["AdvR"][1:end-1] .+ Heat_Budget_true["DiffH"][1:end-1] .+ Heat_Budget_true["DiffR"][1:end-1]  .- Heat_Budget_true["θ"] .+ Heat_Budget_true["θ"][1]
axs.plot(geo_therm[2:end] .- geo_therm[1:end-1], label = "estimate")

axs.legend()
fig

Vθ = 1e-6.*derivative(Heat_Budget_true2["DiffH"], 2.628e+6)  .*  Heat_Budget_approx["Volume"]
plot(6 .* Vθ ./ 20)
