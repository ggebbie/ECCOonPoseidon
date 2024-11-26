#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings
import PyPlot as plt
 
import NaNMath as nm
using .OHC_helper
using PyCall

@pyimport seaborn as sns;
@pyimport pandas as pd;
colors =  sns.color_palette("colorblind")[1:4]
labels = ["ECCO", "ECCO Forcing", "ECCO w/out Forcing Adjust.", "CTRL"]
sns.set_theme(context = "poster", style = "ticks", font_scale = 1.2,
              palette = sns.color_palette("colorblind"));
pygui(false)
cm = pyimport("cmocean.cm");colorway = cm.balance;
include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
tecco = 1992+1/24:1/12:2018

runpath,diagpath = listexperiments(exprootdir());
diagpath["seasonalclimatology_iter0"]
diagpath["seasonalclimatology_iter0_multarr0"] = "/batou/ECCOv4r4/exps/seasonalclimatology_iter0/run_run_multarr0/diags/"
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
marks = expsymbols()
#load in latitude mask 
ocean_mask = wet_pts(Γ)
region = "PAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not", include_bering = false)

#create volume mask
area = readarea(γ)
cell_depths = get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
H = smush(cell_depths, γ); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H
cell_volumes = OHC_helper.get_cell_volumes(area, cell_depths);

#load in temperatures
fig, ax = plt.subplots(1, 1, figsize=(9,12))
ax.vlines(0, 0, 5, color = "grey", alpha = 0.9)

jldopen(datadir("OPT-0015_GH19_PAC.jld2"), "r"; compress = true) do f
    dt_GH19 = 20 / 100 #2015 - 1990 years * 100 to get units percentury
    ax.plot(100 * f["ΔT_GH19"], f["depth_GH19"].*1e-3, color = "k", 
    label = "OPT-15", linestyle = "--")
end

fname = datadir("native/" * region * "_THETA_depth_trends_" * ".jld2")
β = load(fname)["dθ"]

#need to save these fields for faster plotting! They will not change anymore 
@time for (i, expname) in enumerate(["iter129_bulkformula"])#enumerate(keys(shortnames))

    ax.plot(1e2.*β[expname], -z.*1e-3, color = colors[i], label = expname)
    # ax.scatter(1e2.*β, -z.*1e-3, color = colors[i], label = labels[i])
end

ax.set_xlabel(L"^\circ" * "C per century")
ax.set_ylabel("Depth [km]")

ax.tick_params(bottom=true, left=true)
ax.set_xlim(-0.15, 0.15)
ax.invert_yaxis()
ytx = [1, 2,3, 4]
xtx = [-0.1, 0.0, 0.1]
ax.set_yticks(ytx,ytx)
ax.set_xticks(xtx,xtx)
ax.set_ylim(4, 0.5)
ax.yaxis.set_ticks_position("both")
ax.xaxis.set_ticks_position("both")
ax.tick_params(direction="in")
ax.legend(loc = "lower right")
ax.legend(fontsize=24)
ax.set_title("Pacific Ocean Temperature Trends")
sns.move_legend(ax, "lower center",  bbox_to_anchor=(.5, -0.5), ncol=2, frameon=true, borderaxespad=0.)
fig
# fig.savefig(plotsdir() * "/DepthTrendsPACcrop_" * region * ".png",bbox_inches = "tight")

