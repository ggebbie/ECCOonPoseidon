#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, PyCall, 
    Distributions

using .OHC_helper
import PyPlot as plt
@pyimport seaborn as sns;
@pyimport pandas as pd;

include(srcdir("config_exp.jl"))

sns.set_theme(context = "talk", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));

runpath,diagpath = listexperiments(exprootdir());

tecco = 1992+1/24:1/12:2018; tecco = collect(tecco)
uplvl = -2e3; botlvl = -3e3
lvls = findall( botlvl .<= z[:].<= uplvl)
nz = 50; nt = length(tecco)

fnames(expname) = datadir(region * "_" * expname * "_30N_THETA_levels" * ".jld2")

fig,axs=plt.subplots(1,1, sharey = true, figsize = (10, 7.5))
θ = load(fnames("iter129_bulkformula"))["θ"][lvls, :]
ΔV = load(fnames("iter129_bulkformula"))["ΔV"][lvls]
volume_weight_column(x, ΔV) =  sum(x .* ΔV, dims =1) ./ sum(ΔV)

θz = volume_weight_column(θ, ΔV)[:]
E,F = trend_matrices(tecco)
a, b = F*θz
fit_line(x, a, b) = b .* (x .- mean(x)) .+ a
y0b = fit_line(tecco, a, b)
σ = std(θz .- y0b)
axs.plot(tecco, θz, label = "ECCO Output");
axs.plot(tecco, y0b, label = "Trend = " * string(round(1000 * b, digits = 2)) * " m°C per year");
axs.fill_between(tecco, θz .- (2 .*σ), θz .+ (2 .*σ), alpha=0.3, color="k", label = "2σ estimate")
axs.legend(frameon=false)
axs.set_xlabel("time"); axs.set_ylabel("°C")
ax.set_title("Line P02 Mid-Depth Temperature")
fig
fig.savefig(plotsdir("native/P02_MidDepth_Temperatures.png"))



n_iterations = Int(1e7)
trends = zeros(n_iterations)
for k = 1:n_iterations
    t_rand = zeros(Int, 3)
    for (i, year) in enumerate([1994, 2004, 2013])
        random_day(year) = rand(tecco[findall( year .<= collect(tecco).<= year+1)])
        t_rand[i] = findmin(abs.(collect(tecco) .- random_day(year)))[2]
    end
    tecco_sample = tecco[t_rand]
    E,F = trend_matrices(tecco_sample)

    θ_sample = θz[t_rand] .+ rand(Normal(0, σ), 3)

    a, b = F*θ_sample
    trends[k] = 1000 * b
end

fig, ax = plt.subplots()
ax.set_title("Distribution of Mid-Depth P02 Temperature Trends")
sns.histplot(trends, ax = ax, stat = "probability")
tμ = mean(trends)
tσ = std(trends)
ax.axvline(tμ + (2 * tσ), c = "k", linestyle = "--")
ax.axvline(tμ - (2 * tσ), c = "k", linestyle = "--", label = "2σ")
pos_trend_est = sum(trends .> 0.0) / sum(trends .> -Inf)
pos_trend_est = round(pos_trend_est, digits = 2)
ax.scatter(0, nothing, label = "P(β > 0) = " * string(pos_trend_est), alpha = 0.0) 
ax.legend(frameon=false)
ax.set_xlabel("mK per year")
fig
fig.savefig(plotsdir("native/P02_MidDepth_Temp_Trend_Dist.png"))

# fig.savefig(plotsdir("native/generals/temptimeseries" * region * "_all.png"),
#  dpi = 900, bbox_inches = "tight")
# fig