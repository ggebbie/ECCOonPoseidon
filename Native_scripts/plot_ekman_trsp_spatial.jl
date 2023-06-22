#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
Printf, LinearAlgebra, LaTeXStrings, PyCall, GibbsSeaWater
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

region = "PAC56"; 
ocean_mask = OHC_helper.wet_pts(Γ)
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
sns.set_theme(context = "poster", style = "white", font_scale = 1.0,
              palette = sns.color_palette("deep"));

#computed on standard pressures, 
#may have a different value if using pressures diagnosed from ECCO
#https://ecco-v4-python-tutorial.readthedocs.io/VectorCalculus_ECCO_barotropicVorticity.html
runpath,diagpath = listexperiments(exprootdir());

#read in the first time step of S and θ
expname = "iter129_bulkformula"

filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
datafilelist_τ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

nt = length(datafilelist_τ); nz = length(z)
tecco= 1992+1/24:1/12:2018 # ecco years

reg_λ = LLCcropC(λ,γ); reg_ϕ = LLCcropC(ϕ,γ)

reg_mask = LLCcropC(ocean_mask, γ)

mask_array!(x, mask) = (x[mask .== 0] .= NaN)
similar_zeros(x) = zeros(size(x)); 
τxs = similar_zeros(reg_ϕ); 
τys = similar_zeros(reg_ϕ); 
σ0  =similar_zeros(reg_ϕ); 

reg_PAC = LLCcropC(PAC_msk, γ)
dxG = Γ.DXG .* Γ.hFacS[:, 1]; dyG  = Γ.DYG .* Γ.hFacW[:, 1]; 

for tt = 1:nt
    println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
    τname = datafilelist_τ[tt]
    θname = datafilelist_θ[tt]

    @time EXF = γ.read(diagpath[expname]*τname,MeshArray(γ,Float32,15))
    @time θSz = γ.read(diagpath[expname]*θname,MeshArray(γ,Float32,2*nz))
    θz = LLCcropC(θSz[:, 1],γ); Sz = LLCcropC(θSz[:, nz+1],γ)

    τx = EXF[:, 14] .* dxG; τy = EXF[:, 15] .* dyG; 
    τE, τN = rotate_uv(τx, τy, Γ)
    τx_reg = LLCcropC(τE,γ); τy_reg = LLCcropC(τN,γ); 
    
    τxs .+= τx_reg ./ nt
    τys .+= τy_reg ./ nt

    SA = gsw_sa_from_sp.(Sz,0,0,30)
    CT = gsw_ct_from_pt.(SA,θz)
    σ0 .+= gsw_rho.(SA,CT,0) ./ nt #reference at the surface

end

#plot spatial 
proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
fig, ax = plt.subplots( figsize=(15,20), 
subplot_kw=Dict("projection"=> proj0), sharey = true)
bounds = nm.maximum(abs.(σ0 .* reg_mask))
cf = ax.contourf(reg_λ, reg_ϕ,  σ0,
transform=projPC, cmap = cm.curl, extend = "both")

gl = ax.gridlines(crs=projPC, draw_labels=true, linewidth=2, 
                color="gray", alpha=0, linestyle="--")
gl.top_labels = false; gl.right_labels = false
ax.coastlines()
ax.set_title("Meridional Ekman Transport")
ax.set_extent((-326, -60, -60, 70),crs=projPC)
fig.colorbar(cf, ax=ax, orientation = "horizontal", 
fraction=0.03, pad = 0.05, extend = "both", label = L"Sv")
fig
