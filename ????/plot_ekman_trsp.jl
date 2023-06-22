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

include(srcdir("config_exp.jl"))
include("./IsopycnalHelpers.jl")

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

region = "PAC56"; 
ocean_mask = OHC_helper.wet_pts(Γ)
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")

#computed on standard pressures, 
#may have a different value if using pressures diagnosed from ECCO
#https://ecco-v4-python-tutorial.readthedocs.io/VectorCalculus_ECCO_barotropicVorticity.html
runpath,diagpath = listexperiments(exprootdir());

#read in the first time step of S and θ
expname = "iter129_bulkformula"

# # Get list of files for salinity on sigma1
# filelist = searchdir(runpath[expname]*"sigma1/",TSroot) # first filter for state_3d_set1
# datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
τdatafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"

nt = length(τdatafilelist)
tecco= 1992+1/24:1/12:2018 # ecco years

reg_λ = LLCcropC(λ,γ); reg_ϕ = LLCcropC(ϕ,γ)

reg_mask = LLCcropC(ocean_mask, γ)
τxs = zeros(size(reg_ϕ)); 
τys = zeros(size(τxs));
# curlτs = zeros(size(τxs))τxs[reg_mask .== 0] .= NaN
reg_PAC = LLCcropC(PAC_msk, γ)
dxG = Γ.DXG .* Γ.hFacS[:, 1]; dyG  = Γ.DYG .* Γ.hFacW[:, 1]; 

for tt = 1:nt
    println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
    # Tname = datafilelist[tt]
    
    @time EXF = γ.read(diagpath[expname]*τdatafilelist[tt],MeshArray(γ,Float32,15))
    τx = EXF[:, 14] .* dxG; τy = EXF[:, 15] .* dyG; 
    # curlτ = curl(τx,τy,Γ)
    τE, τN = rotate_uv(τx, τy, Γ)

    τx_reg = LLCcropC(τE,γ); τy_reg = LLCcropC(τN,γ); 
    # curlτ_reg = LLCcropC(curlτ,γ);
    
    τxs .+= τx_reg ./ nt
    τys .+= τy_reg ./ nt
    # curlτs .+= curlτ_reg ./ nt
end

ρ = 1045
Ω = 2π/86400
f = 2Ω .* sind.(reg_ϕ)
V_ek = -τxs ./ (f .* ρ)
V_ek[abs.(reg_ϕ) .< 5] .= 0.0
V_ek_PAC = 1e-6 .*  sum(V_ek.* reg_PAC, dims = 1)[:]

fig, ax = plt.subplots( figsize=(15,20))
ax.plot(reg_ϕ[1, :][:], V_ek_PAC)
ax.set_xlim(0, 65)
ax.set_ylim(-10, 10)
fig
# vars = [τxs, τys, V_ek]; var_labels = ["τ_x", "τ_y", "Meridional Ekman Transport"]
# bounds = maximum([nm.maximum(abs.(var .* reg_PAC)) for var in vars])

bounds = nm.maximum(abs.(V_ek .* reg_mask))

# cf = ax.pcolormesh(reg_λ, reg_ϕ,  V_ek, vmin = -bounds, vmax = bounds, 
#     shading="nearest", transform=projPC, rasterized = true, cmap = cm.curl)
levels = -1:0.1:1
cf = ax.contourf(reg_λ, reg_ϕ,  V_ek, vmin = -1, vmax = 1, levels = levels,
transform=projPC, cmap = cm.curl, extend = "both")

gl = ax.gridlines(crs=projPC, draw_labels=true, linewidth=2, 
                color="gray", alpha=0, linestyle="--")
gl.top_labels = false; gl.right_labels = false

# ax.coastlines()
ax.set_title("Meridional Ekman Transport")
ax.set_extent((-326, -60, -60, 70),crs=projPC)

# fig.suptitle("Correlation of dθ/dt with atmospheric variables")
# fig.tight_layout()
fig.colorbar(cf, ax=ax, orientation = "horizontal", 
fraction=0.03, pad = 0.05, extend = "both", label = L"Sv")
fig