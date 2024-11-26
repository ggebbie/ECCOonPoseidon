#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
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
fname = datadir("native/THETA_budget_zonal_" * region *"_1995_2017.jld2")
dθ_iter129 = load(fname)["dθs"]["iter129_bulkformula"]
dθ_iter0 = load(fname)["dθs"]["iter0_bulkformula"]

dθ_plot = Dict()
dθ_plot["Iteration 129"] = dθ_iter129["dθ"] 
dθ_plot["Iteration 0"] = dθ_iter0["dθ"] 
dθ_plot["Iteration 129 minus Iteration 0"] = dθ_plot["Iteration 129"] .- dθ_plot["Iteration 0"]

X = OHC_helper.ma_zonal_sum(ϕ .* area) ./ OHC_helper.ma_zonal_sum(area)
Y = -z[:]
bounds = 0.30
levels = -bounds:0.025:bounds

#plot the differences between temperature trends 
fig, axes = plt.subplots(1, 3, figsize=(30,8), sharey = false)
CF = Any[]
for (i, key) in enumerate(["Iteration 129", "Iteration 0", "Iteration 129 minus Iteration 0"])
    ax = axes[i]
    data = dθ_plot[key] .* (100 * 3.154e+7)
    println(nm.maximum(data))
    cf = ax.contourf(X, Y,  data, cmap = cmo.balance, 
                     vmin = -bounds, vmax = bounds, 
                     levels = levels, extend = "both")   
    push!(CF, cf)
    ax.set_title(key)
end
[ax.set_yticks(1000:1000:5000) for ax in axes]
[ax.set_ylim(1000, 5000) for ax in axes]; 
[ax.invert_yaxis() for ax in axes]
[ax.set_ylabel("Depth [m]") for ax in [axes[1]]]
[ax.set_xticks(-40:10:60) for ax in axes]
[ax.set_xlim(-39, 59) for ax in axes]
[ax.set_xlabel("Latitude") for ax in axes]
fig.colorbar(CF[1], ax = axes, ticks = -bounds:0.1:bounds, orientation = "horizontal", fraction = 0.07, 
label = L"[°C / century ]",pad=0.2)
fig.suptitle("Pacific Ocean Temperature Trends in ECCO [t = 1995-2017] (colors)", y=1)
fig.savefig(plotsdir("native/TempTrendsZonalDiffs_" * region * ".png"), bbox_inches = "tight")

#plot the differences between temperature trends now with sigma2 surfaces overlaid
ρfname = datadir("native/native_sigma2_zonalavg_" * region *"_1995_2017.jld2")
ρ = Dict()
ρ["Iteration 129"] = load(ρfname)["ρθSavg"]["iter129_bulkformula"]["σ2"]
ρ["Iteration 0"] = load(ρfname)["ρθSavg"]["iter0_bulkformula"]["σ2"]
ρ["Iteration 129 minus Iteration 0"] = 0.0 * ρ["Iteration 0"]
fig, axes = plt.subplots(1, 3, figsize=(30,15), sharey = false)
CF = Any[]
for (i, key) in enumerate(["Iteration 129", "Iteration 0", "Iteration 129 minus Iteration 0"])
    ax = axes[i]
    data = dθ_plot[key] .* (100 * 3.154e+7)
    println(nm.maximum(data))
    cf = ax.contourf(X, Y,  data, cmap = cmo.balance, 
                     vmin = -bounds, vmax = bounds, 
                     levels = levels, extend = "both")   
    push!(CF, cf)
    ax.set_title(key)
    CS = ax.contour(X, Y, ρ[key], colors="black", levels = sigma2grid(), linewidths = 1);
    ax.clabel(CS, fontsize=15, inline=true)

end
[ax.set_yticks(1000:1000:5000) for ax in axes]
[ax.set_ylim(1000, 5000) for ax in axes]; 
[ax.invert_yaxis() for ax in axes]
[ax.set_ylabel("Depth [m]") for ax in [axes[1]]]
[ax.set_xticks(-40:10:60) for ax in axes]
[ax.set_xlim(-39, 59) for ax in axes]
[ax.set_xlabel("Latitude") for ax in axes]
fig.colorbar(CF[1], ax = axes, ticks = -bounds:0.1:bounds, orientation = "horizontal", fraction = 0.04, 
label = L"[°C / century ]",pad=0.2)
fig.suptitle("Pacific Ocean Temperature Trends in ECCO [t = 1995-2017] (colors)\n Time-mean σ₂ isopycnals(black contours)", y=1)
fig
fig.savefig(plotsdir("native/TempTrendsZonalDiffswSigma2_" * region * ".png"), bbox_inches = "tight")

#plot the differences between temperature trends WITH isopycnals plotted over
