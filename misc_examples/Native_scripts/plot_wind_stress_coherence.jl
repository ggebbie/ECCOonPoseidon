#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
Printf, LinearAlgebra, LaTeXStrings, PyCall
import PyPlot as plt
import NaNMath as nm
@pyimport matplotlib.animation as anim

cm = pyimport("cmocean.cm");
@pyimport seaborn as sns;

sns.set_theme(context = "poster", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));
              
include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

region = "NPAC"; 
ocean_mask = OHC_helper.wet_pts(Γ)
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
sns.set_theme(context = "poster", style = "white", font_scale = 1.0,
              palette = sns.color_palette("deep"));

#computed on standard pressures, 
#may have a different value if using pressures diagnosed from ECCO
#https://ecco-v4-python-tutorial.readthedocs.io/VectorCalculus_ECCO_barotropicVorticity.html
runpath,diagpath = listexperiments(exprootdir());
τx_dict = Dict(); τy_dict = Dict()
#read in the first time step of S and θ
expname = "iter129_bulkformula"

filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
τdatafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
nt = length(τdatafilelist)
tecco= 1992+1/24:1/12:2018 # ecco years

for tt = 1:nt
    println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
    Tname = τdatafilelist[tt]
    
    @time EXF = γ.read(diagpath[expname]*τdatafilelist[tt],MeshArray(γ,Float32,15))
    τx = EXF[:, 14]; τy = EXF[:, 15]; 

    τE, τN = rotate_uv(τx, τy, Γ)
    τx_reg = LLCcropC(τE,γ); τy_reg = LLCcropC(τN,γ); 
    
    τxs .+= τx_reg ./ nt
    τys .+= τy_reg ./ nt
end
UVtoUEVN
proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

fig, ax = plt.subplots(2, figsize=(15,20), 
subplot_kw=Dict("projection"=> proj0), sharey = true)
vars = [τxs, τys]; var_labels = ["τ_x", "τ_y"]
bounds = maximum([nm.maximum(abs.(var .* reg_PAC)) for var in vars])

levels = -bounds:0.1:bounds
for (i, var) in enumerate(vars)
    cf = ax[i].contourf(reg_λ, reg_ϕ,  var, vmin = -bounds, vmax = bounds,
    transform=projPC, cmap = cm.curl, extend = "both")

    gl = ax[i].gridlines(crs=projPC, draw_labels=true, linewidth=2, 
                    color="gray", alpha=0, linestyle="--")
    gl.top_labels = false; gl.right_labels = false

    # ax.coastlines()
    ax[i].set_title(var_labels[i])
    ax[i].set_extent((-326, -60, -60, 70),crs=projPC)
end
# # fig.suptitle("Correlation of dθ/dt with atmospheric variables")
# # fig.tight_layout()
# fig.colorbar(cf, ax=ax, orientation = "horizontal", 
# fraction=0.03, pad = 0.05, extend = "both", label = L"Sv")
fig