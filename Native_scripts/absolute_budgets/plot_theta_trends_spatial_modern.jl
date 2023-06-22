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
const pyslice=pybuiltin(:slice)
@pyimport seaborn as sns;
@pyimport pandas as pd;
colors =  sns.color_palette("colorblind")[1:4]
labels = ["ECCO", "ECCO Forcing", "ECCO Init. and κ", "CTRL"]
sns.set_theme(context = "poster", style = "ticks", font_scale = 1.2,
              palette = sns.color_palette("colorblind"));
pygui(false)
cm = pyimport("cmocean.cm");colorway = cm.balance;
include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
tecco = 1992+1/24:1/12:2018
runpath,diagpath = listexperiments(exprootdir());


region = "PAC56"; ocean_mask = wet_pts(Γ)

PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
suffix = "2to3"
ocn_reg = LLCcropC(ocean_mask,γ) 

#load in temperatures
expname = "nosfcadjust"
fname = expname * "_θ_trends_2to3km.jld2"
βθ_ECCO = load(datadir(fname))
fname = "OPT-0015_θ_trends_2to3km.jld2"
βθ_OPT0015 = load(datadir(fname))

bounds = 0.15
levels = -bounds:0.01:bounds

proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

fig, ax = plt.subplots(1, 1, figsize=(15,10), subplot_kw=Dict("projection"=> proj0))

reg_mask = LLCcropC(PAC_msk,γ) 
reg_mask[reg_mask .== 0] .=NaN

data = deepcopy(βθ_ECCO)
data["β"][ocn_reg .== 0] .= NaN
data["β"] = data["β"] .* 100 

cf = ax.contourf(data["λ"], data["ϕ"],  data["β"], transform=projPC, 
levels = levels, cmap = cm.balance, vmin = -bounds, vmax = bounds, extend = "both")   
ax.coastlines(resolution="110m")
ax.set_extent((-180, 180, -70, 56),crs=projPC)

gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
gl.xlabels_top = false
gl.ylabels_right = false
ax.set_title("ECCO without Forcing Adjustments")
fig.savefig(plotsdir(expname *"_θ_trends_2to3km.png"), dpi = 1000)

cbar = fig.colorbar(cf, orientation="horizontal",  ticks=-bounds:0.05:bounds,
extend = "both", label = L"^\circ" *"C per century")
fig
fig.savefig(plotsdir(expname *"_θ_trends_2to3km_wcbar.png"), dpi = 1000)

fig, ax = plt.subplots(1, 1, figsize=(15,10), subplot_kw=Dict("projection"=> proj0))
data = deepcopy(βθ_OPT0015)
data["β"] = data["β"] .* 100 

ax.contourf(data["λ"], data["ϕ"],  data["β"], transform=projPC, 
levels = levels, cmap = cm.balance, vmin = -bounds, vmax = bounds, extend = "both")        
ax.coastlines(resolution="110m")
# ax.add_feature(ECCOonPoseidon.cartopy.feature.LAND.with_scale("110m"), facecolor="grey")
# feature = pyimport("cartopy.feature")
# ax.add_feature(feature.LAND, color="lightgray")
ax.set_title("OPT-15")
ax.set_extent((-180, 180, -70, 56),crs=projPC)
gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
                    
gl.xlabels_top = false
gl.ylabels_right = false
fig
# fig.savefig(plotsdir(expname * "OPT0015_θ_trends_2to3km.png"), dpi = 1000)

