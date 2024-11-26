include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings
using .OHC_helper
using ColorSchemes
import NaNMath as nm
import PyPlot as plt
using PyCall

include(srcdir("config_exp.jl"))
@pyimport cmocean as cmocean
@pyimport cmocean.cm as cm
@pyimport seaborn as sns;
@pyimport pandas as pd;

subtract_dict(D, subkey) = Dict(key => value .- D[subkey] for (key, value) in D)
mask_dict(D, lvls) = Dict(key => value[lvls, :] for (key, value) in D)

colors =  sns.color_palette("deep")[1:4]
labels = ["ECCO", "ECCO Forcing", "ECCO Init. Cond. and κ", "CTRL"]
sns.set_theme(context = "poster", style = "white", font_scale = 1.0,
              palette = sns.color_palette("deep"));

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
marks = expsymbols()
 
ocean_mask = wet_pts(Γ)
region = "PAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = OHC_helper.get_cell_volumes(area, cell_depths);

tecco = collect(1992+1/24:1/12:2018)

#load in velocities 
κunadjasted_f = "/batou/eccodrive/files/Version4/Release4/input_init/total_diffkr_r009bit11.bin"
κunadjasted = read_bin(κunadjasted_f, MeshArray(γ,Float32,50))
xx_κadjasted_f = "/batou/eccodrive/files/Version4/Release4/input_init/xx_diffkr.0000000129.data"
xx_κadjasted = γ.read(xx_κadjasted_f,MeshArray(γ,Float32,50))


Y = abs.(z)
X = OHC_helper.ma_zonal_sum(ϕ .* area) ./ OHC_helper.ma_zonal_sum(area)
X_mask = X[X .> -40]
zlev = findall(-4500 .< z[:] .< -2000)

κadjasted_zonal = OHC_helper.ma_zonal_avg(xx_κadjasted, cell_volumes)
# κadjasted_zonal = log.(abs.(κadjasted_zonal))

lightcmap = cmocean.tools.lighten(cm.curl_r, 0.8)

fig, axes = plt.subplots(1,1, figsize = (20, 12.5))
bounds = round(nm.maximum(abs.(κadjasted_zonal[zlev, X .> -40])), digits = 3)
levels = round.(LinRange(-bounds, bounds, 11), digits = 3)

CS = axes.contourf(X_mask, Y, κadjasted_zonal[:, X .> -40], 
cmap=lightcmap, vmin = -bounds, vmax = bounds, 
levels = levels, extend = "both");

axes.set_xlim(-40, 60)
fig.colorbar(CS, orientation = "horizontal", fraction = 0.05, 
label = L"[m^2 / s]")
axes.invert_yaxis()
axes.set_xticks(-40:10:60)
axes.set_xlim(-39, 60)
fig
fig.savefig(plotsdir("DiffusionAdjustmentinPAC.png"), dpi = 1000)

# fig.savefig(plotsdir("ECCOWaterMassDifference_" * titles[i] * ".png"), dpi = 1000)

κunadjasted_f = "/batou/eccodrive/files/Version4/Release4/input_init/total_diffkr_r009bit11.bin"
κunadjasted = read_bin(κunadjasted_f, MeshArray(γ,Float32,50))
κadjastedtotal_zonal = OHC_helper.ma_zonal_avg(κunadjasted, cell_volumes)
κadjastedtotal_zonal = 1e6 .* κadjastedtotal_zonal

fig, axes = plt.subplots(1,1, figsize = (20, 12.5))
bounds = round.(nm.extrema(κadjastedtotal_zonal[zlev, X .> -40]), digits = 1)

levels = round.(LinRange(bounds[1], bounds[2], 5), digits = 1) 

# good_color = sns.color_palette("Spectral", as_cmap=true)
# good_color = sns.color_palette("vlag", as_cmap=true)
lightcmap = cmocean.tools.lighten(cm.matter, 0.9)
CS = axes.contourf(X_mask, Y, κadjastedtotal_zonal[:, X .> -40], 
cmap=lightcmap, vmin = bounds[1] - 1, vmax = bounds[2] + 1, 
levels = levels, extend = "both");
# CS = axes.contour(X_mask, Y, κadjastedtotal_zonal[:, X .> -40], colors="black", levels = levels);
# axes.clabel(CS, fontsize=25, inline=true)
axes.set_xlim(-40, 60)
# axes.set_title(titles[i])
axes.invert_yaxis()
axes.set_xticks(-40:10:60)
axes.set_xlim(-39, 60)
fig.colorbar(CS, orientation = "horizontal", fraction = 0.05, 
label = L"[mm^2 / s]")
fig
fig.savefig(plotsdir("DiffusioninPAC.png"), dpi = 1000)

total = κunadjasted