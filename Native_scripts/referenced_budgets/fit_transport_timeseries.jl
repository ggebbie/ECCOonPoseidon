#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, PyCall, LsqFit
using .OHC_helper
import PyPlot as plt
@pyimport seaborn as sns;
@pyimport pandas as pd;

colors =  sns.color_palette("colorblind")[1:4]
labels = ["ECCO", "ECCO Forcing", "ECCO Init. and κ", "CTRL"]

sns.set_theme(context = "poster", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));

include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());
tecco = 1992+1/24:1/12:2018

ocean_mask = OHC_helper.wet_pts(Γ)
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; region, extent = "not")
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"


fig,axs=plt.subplots(1,1, sharey = true, figsize = (10, 7.5))
alpha = [1, 1, 1, 1]
colors =  sns.color_palette("colorblind")[1:5]
colors = colors[[1, 3, 5, 4]]

exps = ["nosfcadjust", "noinitadjust", "iter0_bulkformula", "iter129_bulkformula"]
exps = ["seasonalclimatology", "seasonalclimatology_iter0"]

# labels = ["Iteration 129", "Initial 129", "Initial 0", "Iteration 0"]

@time for (i, expname) in enumerate(exps)

    fname = datadir("native/" * expname * region * "_THETA_budget_ref_" * suffix * ".jld2")
    θ = load(fname)["dθ"]["θ"]
    p0 = [0.1, 0.01, 1.5]
    t = (collect(tecco) .- tecco[1]);
    model(t, p) = (p[1] .* exp.(-t .* p[2])) .+ p[3] 
    fit = curve_fit(model, t, θ, p0; lower = [0.01, 0.005, 0.00], upper = [1.0, 1, Inf], autodiff=:forwarddiff)
    fit_estimate = model(0:26, fit.param)
    axs.plot(t, θ, alpha  = alpha[i], color = colors[i], label = expname)
    axs.plot(0:26, fit_estimate, alpha  = alpha[i], color = "k"); 
    println(fit.converged); println(round.(fit.param, digits = 5));
    println(stderror(fit))
end
axs.legend()
fig


fig,axs=plt.subplots(1,1, sharey = true, figsize = (10, 7.5))
alpha = [1, 1, 1, 1]
colors =  sns.color_palette("colorblind")[1:5]
colors = colors[[1, 3, 5, 4]]

exps = ["nosfcadjust", "noinitadjust", "iter0_bulkformula", "iter129_bulkformula"]
exps = ["seasonalclimatology", "seasonalclimatology_iter0"]

# labels = ["Iteration 129", "Initial 129", "Initial 0", "Iteration 0"]

@time for (i, expname) in enumerate(exps)

    fname = datadir("native/" * expname * region * "_THETA_budget_ref_" * suffix * ".jld2")
    θ = load(fname)["dθ"]["θ"]
    p0 = [0.1, 0.01, 1.5]
    t = (collect(tecco) .- tecco[1])
    model(t, p) = (p[1] .* exp.(-t .* p[2])) .+ p[3] 
    fit = curve_fit(model, t, θ, p0; lower = [0.01, 0.005, 0.00], upper = [1.0, 1, minimum(θ)], autodiff=:forwarddiff)
    t_est = 0:10000; fit_estimate = model(t_est, fit.param); 
    axs.plot(t, θ, alpha  = alpha[i], color = colors[i], label = expname)
    axs.plot(t_est, fit_estimate, alpha  = alpha[i], color = "k"); 
    # println(fit.converged); println(round.(fit.param, digits = 4));
    # println(stderror(fit))
    # println(confidence_interval(fit, 0.05))
end
axs.legend()
fig
