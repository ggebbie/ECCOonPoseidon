#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
using .OHC_helper
import NaNMath as nm
import PyPlot as plt
include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");
@pyimport seaborn as sns;
@pyimport pandas as pd;

colors =  sns.color_palette("colorblind")[1:4]
labels = ["ECCO", "ECCO Forcing", "ECCO Init. and κ", "CTRL"]

sns.set_theme(context = "poster", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));
#using . is faster 
(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());

# abbreviations for each experiment for labels, etc.
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
marks = expsymbols()
nexp = length(shortnames) # number of experiments
 
ocean_mask = OHC_helper.wet_pts(Γ)
region = "PAC56"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));

tecco = 1992+1/24:1/12:2018
nz = 50; tstart = 12*3; tstop = 312; nt = tstop - tstart + 1;
E,F = trend_matrices(tecco[tstart:tstop])
#pre-allocate
# GTF = OHC_helper.get_geothermalheating(γ, Γ)
# GTF_zonal = ma_zonal_avg(GTF, cell_volumes)
ΔV = ma_zonal_sum(cell_volumes)
ΔV[ΔV .== 0] .= NaN

sum_fluxes(d) = sum(d[varn] for varn in ["κxyθ", "uvθ", "wθ", "κzθ"])
expname = "iter129_bulkformula"
fname = datadir("THETA_budget_zonal_" * region * "_" * expname *".jld2")
dθ = load(fname)["dθ"]

dθ_plot = Dict()
dθ_plot["θ Trend"] = deepcopy(dθ["dθ"] )
dθ_plot["∇(⟨uθ⟩)"] = dθ["uvθ"] .+ dθ["wθ"]
dθ_plot["∇²(κ⟨θ⟩)"] = dθ["κxyθ"] .+ dθ["κzθ"]
dθ_plot["d⟨θ⟩/dt"] = sum_fluxes(dθ)
dθ_plot["Residual"] = dθ_plot["θ Trend"] .- dθ_plot["d⟨θ⟩/dt"]

nm.extrema(dθ_plot["Residual"][30:end, :].* (100 * 3.154e+7) )
X = OHC_helper.ma_zonal_sum(ϕ .* area) ./ OHC_helper.ma_zonal_sum(area)
Y = -z[:]
bounds = 0.30
levels = -bounds:0.01:bounds

fig, axes = plt.subplots(1, 3, figsize=(30,8), sharey = true)
CF = Any[]
for (i, key) in enumerate(["d⟨θ⟩/dt", "∇(⟨uθ⟩)", "∇²(κ⟨θ⟩)"])
    ax = axes[i]
    data = dθ_plot[key] .* (100 * 3.154e+7)
    println(nm.extrema(data[30:end, :]))
    cf = ax.contourf(X, Y,  data, cmap = cmo.balance, 
                     vmin = -bounds, vmax = bounds, 
                     levels = levels, extend = "both")   
    push!(CF, cf)
    ax.set_title(key)
end
axes[1].set_ylim(1000, 5000); axes[1].invert_yaxis(); 
axes[1].set_ylabel("Depth [m]")
[ax.set_xticks(-40:10:60) for ax in axes]
[ax.set_xlim(-39, 59) for ax in axes]
[ax.set_xlabel("Latitude") for ax in axes]
fig.colorbar(CF[1], ax = axes, orientation = "horizontal", fraction = 0.07, 
label = L"[°C / century ]",pad=0.2)
fig.savefig(plotsdir("HeatBudgetZonal_" * region * "_" * expname * ".png"), bbox_inches = "tight")
fig

fig, axes = plt.subplots(1, 3, figsize=(30,8), sharey = true)
CF = Any[]
for (i, key) in enumerate(["θ Trend", "d⟨θ⟩/dt", "Residual"])
    ax = axes[i]
    data = dθ_plot[key] .* (100 * 3.154e+7)
    println(nm.extrema(data[30:end, :]))
    cf = ax.contourf(X, Y,  data, cmap = cmo.balance, 
                     vmin = -bounds, vmax = bounds, 
                     levels = levels, extend = "both")   
    push!(CF, cf)
    ax.set_title(key)
end
axes[1].set_ylim(1000, 5000); axes[1].invert_yaxis(); 
axes[1].set_ylabel("Depth [m]")
[ax.set_xticks(-40:10:60) for ax in axes]
[ax.set_xlim(-39, 59) for ax in axes]
[ax.set_xlabel("Latitude") for ax in axes]
fig.colorbar(CF[1], ax = axes, orientation = "horizontal", fraction = 0.07, 
label = L"[°C / century ]",pad=0.2)
fig.savefig(plotsdir("HeatBudgetZonalDifference_" * region * "_" * expname * ".png"), bbox_inches = "tight")
fig