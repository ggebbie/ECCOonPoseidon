#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools,
    PyCall
using .OHC_helper
import PyPlot as plt
@pyimport seaborn as sns;
@pyimport pandas as pd;

include(srcdir("config_exp.jl"))


sns.set_theme(context = "talk", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());

tecco = 1992+1/24:1/12:2018

ocean_mask = OHC_helper.wet_pts(Γ)
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = OHC_helper.get_cell_volumes(area, cell_depths);
uplvl = -2e3; botlvl = -3e3
lvls = findall( botlvl .<= z[:].<= uplvl)

nz = 50

ΔV = zeros(Float32, nz)
[ΔV[k] = Float32(sum(cell_volumes[:, k])) for k=lvls]

V = sum(ΔV)

integrate(x) = cumsum(hcat([0.0], x .* 2.628e+6), dims=2)
volume_weight_column(x) =  sum(x .* ΔV, dims =1) ./ V 
sum_fluxes(d) = sum(d[varn][:] for varn in ["κxyθ", "uvθ", "wθ", "κzθ",  "GTF"])
anomaly(x) = x .- mean(x)

fname(expname) = datadir(region * "_" * expname * "_THETA_levels" * ".jld2")

fig,axs=plt.subplots(1,1, sharey = true, figsize = (10, 7.5))
# labels = ["Iteration 129", "Initial 129", "Initial 0", "Iteration 0"]
alpha = [0.5, 1, 1, 0.5]

exps = ["climatological_tau", "nosfcadjust"]
load(fname(exps[1]))["θ"]
@time for (i, expname) in enumerate(exps)
    println(expname)
    θ = load(fname(expname))["θ"]
    θz = volume_weight_column(θ)[:]

    axs.plot(tecco, θz, label = expname, alpha  = alpha[i]); 
end
axs.legend(frameon=false)
axs.set_xlabel("time")
axs.set_ylabel("cK")
fig

# fig.savefig(plotsdir("native/generals/temptimeseries" * region * "_all.png"),
#  dpi = 900, bbox_inches = "tight")
# fig