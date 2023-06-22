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

fname = datadir("native/" * region * "_WTRSP_hovmoller.jld2")
W_dict = load(fname)["W_dict"]
W_dict["Difference"] = Dict(); W_dict["Difference"]["Wr"] = W_dict["iter129_bulkformula"]["Wr"] .- W_dict["seasonalclimatology"]["Wr"]
fig, axes = plt.subplots(1, 3, figsize = (20, 3), sharey = true)
levels = -10:2:10
fig.suptitle("Net Residual Vertical Transport in North Pacific")
CF = []
for (i, expname) in enumerate(["iter129_bulkformula", "seasonalclimatology", "Difference"])
    Wresid = W_dict[expname]["Wr"]
    cf = axes[i].contourf(tecco, z, 1e-6 .* (Wresid), vmin = -10, vmax = 10, cmap = cmo.balance, levels = levels, extend = "both")
    axes[i].set_title(expname)
    push!(CF, cf)
end
fig
fig.colorbar(CF[1], ax = axes, orientation = "horizontal", fraction = 0.02, label = "[Sv]")
fig

fname = datadir("native/" * region * "_NTRSP_hovmoller.jld2")
N_dict = load(fname)["N_dict"]
N_dict["Difference"] = Dict(); N_dict["Difference"]["Nr"] = N_dict["iter129_bulkformula"]["Nr"] .- N_dict["seasonalclimatology"]["Nr"]
fig, axes = plt.subplots(1, 3, figsize = (20, 3), sharey = true)
levels = -1:0.2:1
fig.suptitle("Net Residual Horizontal Transport in North Pacific")
CF = []
for (i, expname) in enumerate(["iter129_bulkformula", "seasonalclimatology", "Difference"])
    Nresid = N_dict[expname]["Nr"]
    cf = axes[i].contourf(tecco, z, 1e-6 .* (Nresid), vmin = -1, vmax = 1, cmap = cmo.balance, levels = levels, extend = "both")
    axes[i].set_title(expname)
    push!(CF, cf)
end
fig.colorbar(CF[1], ax = axes, orientation = "horizontal", fraction = 0.1, label = "[Sv]")
fig

fig, axes = plt.subplots(1, 3, figsize = (20, 3), sharey = true)
levels = -10:2:10
fig.suptitle("Net Residual Horizontal Transport in North Pacific")

fig, axes = plt.subplots(1, 3, figsize = (20, 3), sharey = true)
fig
for (i, expname) in enumerate(["iter129_bulkformula", "seasonalclimatology"])
    Nresid = N_dict[expname]["Nr"]
    println(1e-6 .* maximum(Nresid))
    axes[i].plot(tecco, 1e-6 .* sum(Nresid[30:end], dims = 1)[:])
    axes[i].set_title(expname)
    # axes[i].set_ylim((-0.25, 0.25))
end
fig

fig, axes = plt.subplots(1, 3, figsize = (10, 5), sharey = true)
levels = -10:2:10
fig.suptitle("Net Residual Vertical Transport in North Pacific")
CF = []
W_dict["Difference"] = Dict(); W_dict["Difference"]["Wr"] = W_dict["iter129_bulkformula"]["Wr"] .- W_dict["seasonalclimatology"]["Wr"]

for (i, expname) in enumerate(["iter129_bulkformula", "iter0_bulkformula", "Difference"])
    Wresid = W_dict[expname]["Wr"]
    Wresid_climatology = hcat([mean(Wresid[:, i:12:end], dims = 2)[:] for i in 1:12]...)
    cf = axes[i].contourf(1:12, z, 1e-6 .* (Wresid_climatology), vmin = -10, vmax = 10, cmap = cmo.balance, levels = levels, extend = "both")
    axes[i].set_title(expname)
    axes[i].set_xlabel("time")
    push!(CF, cf)
end
fig.tight_layout()
fig
fig.colorbar(CF[1], ax = axes, orientation = "horizontal", fraction = 0.05, label = "[Sv]")
fig