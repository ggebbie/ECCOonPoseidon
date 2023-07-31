#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics,
NCDatasets, Printf, 
DataFrames, LaTeXStrings
import NaNMath as nm
using .OHC_helper
import PyPlot as plt
using PyCall
const pyslice=pybuiltin(:slice)
@pyimport seaborn as sns;
@pyimport pandas as pd;
colors =  sns.color_palette("colorblind")[1:4]
labels = ["ECCO", "ECCO Forcing", "ECCO Init. and κ", "CTRL"]
sns.set_theme(context = "talk", style = "ticks", font_scale = 1.2,
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

fname = "iter129_bulkformula_θ_trends_0to700.jld2"
βθ_ECCO = load(datadir(fname))
# fname = "OPT-0015_θ_trends_2to3km.jld2"
# βθ_OPT0015 = load(datadir(fname))


proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

fig, ax = plt.subplots(1, 1, figsize=(15,10), subplot_kw=Dict("projection"=> proj0))
data = deepcopy(βθ_ECCO); data["β"][ocn_reg .== 0] .= NaN


cf = ax.pcolormesh(data["λ"], data["ϕ"],  data["β"], transform=projPC,
cmap = cm.balance, vmin = -0.1, vmax = 0.1)   
ax.coastlines(resolution="110m")

gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
gl.xlabels_top = false
gl.ylabels_right = false

cbar = fig.colorbar(cf, orientation="horizontal",
extend = "both", label = L"^\circ" *"C per year", fraction = 0.04)
ax.set_title("Upper Ocean Temperature Trends in ECCO \n (2003 - 2012)")

fig
fig.savefig(plotsdir("UpperOceanTrendsECCO.png"), dpi = 400)
