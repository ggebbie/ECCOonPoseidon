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

sig1grid = sigma1grid()
nσ = length(sig1grid)

#get root for salinity on sigma1
#THETA for iter129
TSroot = "THETA_on_sigma1" 

# Get list of files for salinity on sigma1
filelist = searchdir(runpath[expname]*"sigma1/",TSroot) # first filter for state_3d_set1
datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
nt = length(datafilelist)

tecco= 1992+1/24:1/12:2018 # ecco years

reg_mask = LLCcropC(PAC_msk,γ)
reg_mask[reg_mask .== 0] .= NaN

reg_λ = LLCcropC(λ,γ); reg_ϕ = LLCcropC(ϕ,γ)
junk,ilonilat = findmin((reg_ϕ.-0).^2 .+ (reg_λ.-0).^2)        

Us = zeros(size(reg_mask)..., nt)
bs = zeros(nt)

for tt = 1:nt
    println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
    Tname = datafilelist[tt]
    @time θσ1 = γ.read(runpath[expname]*"sigma1/"*Tname,MeshArray(γ,Float32,nσ))
    θ_reg = LLCcropC(θσ1[:, 66],γ)
    Us[:, :, tt] .= θ_reg
    bs[tt] = θ_reg[ilonilat]
end

Up = Us .- mean(Us, dims = 3)
bp = bs .- mean(bs)
σU = std(Up, dims = 3)[:, :, 1]; σb = std(bp)
σUσb = σU .* σb

c = zeros(size(reg_mask)...)
[c[:, :] .+= (Up[:, :, tt] .* bp[tt] ./ nt) ./ σUσb for tt in 1:nt]

fig, ax = plt.subplots(figsize=(15,10), subplot_kw=Dict("projection"=> proj0))
bounds = nm.maximum(abs.(c)) 

proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

ax.coastlines()
cf = ax.pcolormesh(reg_λ, reg_ϕ,  c, 
                   vmin = -bounds, vmax = bounds, shading="nearest", 
                   transform=projPC, rasterized = true, cmap = cm.curl)

fig.colorbar(cf, ax=ax, orientation = "horizontal", 
fraction=0.04, pad = 0.1, extend = "both", label = "ρ")

ax.scatter(reg_λ[ilonilat], reg_ϕ[ilonilat], s=100, marker = "*",c = "blue")

gl = ax.gridlines(crs=projPC, draw_labels=true, linewidth=2, 
                  color="gray", alpha=0, linestyle="--")
gl.top_labels = false; gl.right_labels = false

ax.set_title("Correlation")

fig.tight_layout()
fig