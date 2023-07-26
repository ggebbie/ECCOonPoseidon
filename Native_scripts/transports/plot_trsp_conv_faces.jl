#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
using .OHC_helper
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");
@pyimport seaborn as sns;

colors =  sns.color_palette("deep")[1:4]
sns.set_theme(context = "talk", font_scale = 1.0,
              palette = sns.color_palette("deep"));
sns.set_style("white")
(ϕ,λ) = latlonC(γ)
area = readarea(γ)
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
tecco = 1992+1/24:1/12:2018

ocean_mask = OHC_helper.wet_pts(Γ)
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; region, extent = "not", include_bering = false)
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= z[:].<= uplvl)

fname = datadir("native/" * region * "_TRSP_" * ".jld2")
transports_dict = load(fname)["transports_dict"]
fig, axes = plt.subplots(1, 3, figsize = (10, 8), sharey = true)
alphas = [0.5, 1]
nt = 312
labels = ["Initial 0", "Initial 129"]

for (i, var) in enumerate(keys(transports_dict))
    println(var)
    withoutBering = 1
    Vout = withoutBering .* transports_dict[var].Vout
    axes[1].plot(tecco[1:nt], 1e-6 .* transports_dict[var].Vin, label = var, c = colors[i], alpha = alphas[i]); 
    axes[2].plot(tecco[1:nt], 1e-6 .* transports_dict[var].Wtop, label = var, c = colors[i], alpha = alphas[i]); 
    axes[3].plot(tecco[1:nt], 1e-6 .* transports_dict[var].Wbot, label = var, c = colors[i], alpha = alphas[i]); 

    axes[1].set_title( L"V_{in}"); axes[2].set_title( L"W_{out}"); ; axes[2].set_title( L"W_{in}"); 
    fig.suptitle(region * " Volume Fluxes [z = 2 - 3 km]")
end
[a.legend() for a in axes]
fig