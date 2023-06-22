#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings,
Plots
gr()
import NaNMath as nm
using .OHC_helper
import PyPlot as plt
using PyCall

@pyimport seaborn as sns;
@pyimport pandas as pd;
colors =  sns.color_palette("colorblind")[1:4]
labels = ["ECCO", "ECCO Forcing", "ECCO Init. Cond.", "CTRL"]
sns.set_theme(context = "poster", style = "ticks", font_scale = 1.2,
              palette = sns.color_palette("colorblind"));
pygui(false)
cm = pyimport("cmocean.cm");colorway = cm.balance;
include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
tecco = 1992+1/24:1/12:2018
runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
marks = expsymbols()
#load in latitude mask 
ocean_mask = wet_pts(Γ)
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
suffix = "2to3"
uplvl = -2e3; botlvl = -3e3
lvls = findall( botlvl .<= z[:].<= uplvl)

#create volume mask
area = readarea(γ)
cell_depths = get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
H = smush(cell_depths); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H
cell_volumes = get_cell_volumes(area, cell_depths);


function volume_average_by_depth(ds::MeshArray, ΔV::MeshArray, γ)
    nz = size(ds, 2)
    vol_avg = zeros(Float32, nz)
    V = zeros(Float32, nz)

    for k=1:nz
        V[k] = Float32(sum(ΔV[:, k]))
    end

    for ff=1:5, k=1:nz
        vol_avg[k] += sum(ds[ff, k] .* ΔV[ff, k]) / V[k]
    end

    return vol_avg
end

#load trend trend_matrices, F is LS estimator
E,F = trend_matrices(tecco)

#load in temperatures
fig, ax = plt.subplots(1, 1, figsize=(9,12))
ax.vlines(0, 0, 5, color = "grey", alpha = 0.9)

# jldopen(datadir("OPT-0015_GH19.jld2"), "r"; compress = true) do f
#     dt_GH19 = 1.2 #1.2 century is the sampling period of ΔT
#     ax.plot(f["ΔT_GH19"] ./ dt_GH19, f["depth_GH19"].*1e-3, color = "k", 
#     label = "GH19", linestyle = "--")
# end

@time for (i, expname) in enumerate(keys(shortnames))#enumerate(keys(shortnames))
    filelist = searchdir(diagpath[expname], "state_2d_set1") # first filter for state_3d_set1
    datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
    datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    β = zeros(50)
    nt = length(datafilelist_S)
    # F =  F
    tt = 1
    fnameS = datafilelist_S[tt]
    fnameθ = datafilelist_θ[tt]
    sθ = extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_H)
    sθavg = volume_average_by_depth(sθ, cell_volumes, γ)
    ax.plot(sθavg, -z.*1e-3, color = colors[i], label = labels[i])
        # ax.scatter(1e2.*β, -z.*1e-3, color = colors[i], label = labels[i])
    
end

ax.set_xlabel(L"^\circ" * "C per century")
ax.set_ylabel("Depth [km]")

ax.tick_params(bottom=true, left=true)
# ax.set_ylim(1, 4)
ax.set_xlim(0.8, 2)

ax.invert_yaxis()
# ytx = [2,3]
# xtx = [-0.1, 0.0, 0.1]
# ax.set_yticks(ytx,ytx)
# ax.set_xticks(xtx,xtx)
ax.yaxis.set_ticks_position("both")
ax.xaxis.set_ticks_position("both")
ax.tick_params(direction="in")
ax.legend(loc = "lower right")
ax.legend(fontsize=24)
ax.set_title("North Pacific Ocean Temperature Trends")
# sns.move_legend(ax, "lower center",  bbox_to_anchor=(.5, -0.5), ncol=2, frameon=true, borderaxespad=0.)
fig
fig.savefig(plotsdir() * "/AGU_Plots/ThetaAvgPAC_" * region * ".pdf",bbox_inches = "tight")

