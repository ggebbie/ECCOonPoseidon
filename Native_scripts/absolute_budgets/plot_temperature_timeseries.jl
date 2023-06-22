#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings,
    PyCall, ColorSchemes
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

ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
marks = expsymbols()
nexp = length(shortnames) # number of experiments

tecco = 1992+1/24:1/12:2018

ocean_mask = OHC_helper.wet_pts(Γ)
region = "NPAC"; 
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));
uplvl = -2e3; botlvl = -3e3
lvls = findall( botlvl .<= z[:].<= uplvl)

# GTF = OHC_helper.get_geothermalheating(γ, Γ) .* PAC_msk
nz = 50

ΔV = zeros(Float32, nz)
GTF_z=  zeros(Float32, nz)
[ΔV[k] = Float32(sum(cell_volumes[:, k])) for k=lvls]

# [GTF_z[k] = sum(GTF[:, k]) for k=lvls]
GTF_z = repeat(GTF_z, 1, 312)
V = sum(ΔV)

integrate(x) = cumsum(hcat([0.0], x .* 2.628e+6), dims=2)
volume_weight_column(x) =  sum(x .* ΔV, dims =1) ./ V 
sum_fluxes(d) = sum(d[varn][:] for varn in ["κxyθ", "uvθ", "wθ", "κzθ",  "GTF"])
anomaly(x) = x .- mean(x)

fname = datadir("native/" * region * "_THETA_levels" * ".jld2")
θ = load(fname)["θ"]
θ_budg_int = Dict()

for var in keys(θ)
    tmp = volume_weight_column(θ[var])
    θ_budg_int[var] = tmp[:]
end

sns.set_theme(context = "talk", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));

fig,axs=plt.subplots(1,1, sharey = true, figsize = (7.5, 7.5))
labels = ["Iteration 129", "Initial 129", "Initial 0", "Iteration 0"]
alpha = [0.5, 1, 1, 0.5]
colors =  sns.color_palette("colorblind")[1:5]
colors = colors[[1, 3, 5, 4]]
@time for (i, expname) in enumerate(keys(θ_budg_int))
    println(expname)
    θz = vec(θ_budg_int[expname])

    axs.plot(tecco, 100 .* (θz .- θz[1]), label = labels[i], alpha  = alpha[i], color = colors[i]); 
end
axs.legend(frameon=false)
axs.set_xlabel("time")
axs.set_ylabel(["cK"])
fig

fig.savefig(plotsdir("native/generals/temptimeseries" * region * "_all.png"),
 dpi = 900, bbox_inches = "tight")
# fig