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
filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
τdatafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"

nt = length(datafilelist)
tecco= 1992+1/24:1/12:2018 # ecco years

reg_λ = LLCcropC(λ,γ); reg_ϕ = LLCcropC(ϕ,γ)
junk,ilonilat = findmin((reg_ϕ.-0).^2 .+ (reg_λ.-0).^2)        

τxs = zeros(size(reg_ϕ)..., nt)
τys = zeros(size(τxs))
curlτs = zeros(size(τxs))

bs = zeros(nt)

for tt = 1:nt
    println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
    Tname = datafilelist[tt]
    @time θσ1 = γ.read(runpath[expname]*"sigma1/"*Tname,MeshArray(γ,Float32,nσ))
    
    @time EXF = γ.read(diagpath[expname]*τdatafilelist[tt],MeshArray(γ,Float32,15))
    τx = EXF[:, 14]; τy = EXF[:, 15]; curlτ = curl(τx,τy,Γ)
    τE, τN = rotate_uv(τx, τy, Γ)

    τx_reg = LLCcropC(τE,γ); τy_reg = LLCcropC(τN,γ); 
    curlτ_reg = LLCcropC(curlτ,γ);
    
    θ_reg  = LLCcropC(θσ1[:, 66],γ)

    τxs[:, :, tt] .= τx_reg
    τys[:, :, tt] .= τy_reg
    curlτs[:, :, tt] .= curlτ_reg

    bs[tt] = θ_reg[ilonilat]
end

function spatial_correlation(bs, Us)
    nt = length(bs)
    Up = Us .- mean(Us, dims = 3)
    bp = bs .- mean(bs)
    σU = std(Up, dims = 3)[:, :, 1]; σb = std(bp)
    σUσb = σU .* σb
    
    c = zeros(size(Us[:, :, 1])...)
    [c[:, :] .+= (Up[:, :, tt] .* bp[tt] ./ nt) ./ σUσb for tt in 1:nt]
    return c
end


proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

fig, axes = plt.subplots(ncols= 3, figsize=(15,20), 
subplot_kw=Dict("projection"=> proj0), sharey = true)
vars = [τxs, τys, curlτs ]; var_labels = ["τ_x", "τ_y", "curl τ"]
bounds = maximum([nm.maximum(abs.(var)) for var in vars])
cfs = Any[]

dtbs =  (bs[2:end] .- bs[1:end-1]) ./ 2.628e+6 #time derivative

for (i, ax) in enumerate(axes)
    # linear interpolate 
    var = (vars[i][:, :, 1:end-1] .+ vars[i][:, :, 2:end]) ./ 2
    cs = spatial_correlation(dtbs, var)
    cf = ax.pcolormesh(reg_λ, reg_ϕ,  cs, 
                   vmin = -bounds, vmax = bounds, shading="nearest", 
                   transform=projPC, rasterized = true, cmap = cm.curl)
    push!(cfs, cf)
    ax.scatter(reg_λ[ilonilat], reg_ϕ[ilonilat], s=100, marker = "*",c = "blue")
    gl = ax.gridlines(crs=projPC, draw_labels=true, linewidth=2, 
                  color="gray", alpha=0, linestyle="--")
    gl.top_labels = false; gl.right_labels = false
    if i > 1
        gl.left_labels = false
    end
    # ax.coastlines()
    ax.set_title(var_labels[i])
    ax.set_extent((-326, -60, -60, 70),crs=projPC)

end
# fig.suptitle("Correlation of dθ/dt with atmospheric variables")
# fig.tight_layout()
fig.colorbar(cfs[1], ax=axes, orientation = "horizontal", 
fraction=0.01, pad = 0.05, extend = "both", label = L"\rho")
fig