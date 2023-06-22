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
sns.set_theme(context = "notebook", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));
colors =  sns.color_palette("deep")[1:4]

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
area_mask = area .* PAC_msk

for expname in ["iter129_bulkformula", "iter0_bulkformula", "seasonalclimatology"]
    fname = datadir("native/" * region * "_TRSP_hovmoller.jld2")
    W_dict = load(fname)["W_dict"][expname]

    fname = datadir("native/" * region * "_areamean_τcurl.jld2")
    curlτ_dict = load(fname)["curlτ_dict"][expname]

    fname = datadir("native/" * region * "_NTRSP_hovmoller.jld2")
    N_dict = load(fname)["N_dict"][expname]

    fig, axes = plt.subplots(2, 1, figsize = (10, 8), sharey = false, sharex = true)
    levels = -10:2:10
    Wresid = W_dict["We"]
    cf = axes[1].contourf(tecco, z, 1e-6 .* (Wresid), vmin = -10, vmax = 10, cmap = cmo.balance, levels = levels, extend = "both")
    axes[2].plot(tecco, curlτ_dict, linewidth = 2, color = "k", alpha = 0.8)
    axes[2].axhline(mean(curlτ_dict), linewidth = 1, color = "k", alpha = 0.5)
    axes[1].set_title("Net Vertical Transport in North Pacific")
    axes[2].set_title("Area-Averaged North Pacific Wind Stress Curl")
    fig.tight_layout()
    fig.colorbar(cf, ax = axes, orientation = "horizontal", fraction = 0.04, label = "[Sv]")
    axes[1].set_ylabel("Depth [m]")
    axes[2].set_ylabel("N per m³")
    fig.savefig(plotsdir("native/WTRSP_τcurl_Hovmoller" * expname * ".png"), dpi = 400)

    fig, axes = plt.subplots(1, figsize = (10, 4), sharey = false, sharex = false)
    levels = -2:0.2:2
    Nresid = N_dict["Nr"]
    cf = axes.contourf(tecco, z, 1e-6 .* (Nresid), vmin = -2, vmax = 2, cmap = cmo.balance, levels = levels, extend = "both")
    axes.set_title("Net Meridional Transport in North Pacific")
    fig.tight_layout()
    fig.colorbar(cf, ax = axes, orientation = "horizontal", fraction = 0.04, label = "[Sv]")
    axes.set_ylabel("Depth [m]")
    fig.savefig(plotsdir("native/NTRSP_Hovmoller" * expname * ".png"), dpi = 400)

    fig, axes = plt.subplots(1, figsize = (4, 8), sharey = false, sharex = false)
    axes.plot(1e-6 .* mean(Nresid, dims = 2), z, label = "Meridional Transport")
    axes.plot(1e-6 .* mean(Wresid, dims = 2), z, label = "Vertical Transport")
    axes.set_title("Time-mean transport in North Pacific")
    axes.legend()
    axes.set_ylabel("Depth [m]")
    axes.set_xlabel("Sv")
    axes.set_ylim(-5000, -500)
    fig.savefig(plotsdir("native/TRSP_timemean" * expname * ".png"), dpi = 400)
end