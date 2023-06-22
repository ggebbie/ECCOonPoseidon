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
#reshape λ for plotting 
λ_offset = deepcopy(λ)
for ff ∈ [3,4] 
    λ_offset[ff][λ_offset[ff] .<= 0 ] = λ_offset[ff][λ_offset[ff] .<= 0 ] .+ 360
end


area = readarea(γ)

#computed on standard pressures, 
#may have a different value if using pressures diagnosed from ECCO
#https://ecco-v4-python-tutorial.readthedocs.io/VectorCalculus_ECCO_barotropicVorticity.html
runpath,diagpath = listexperiments(exprootdir());

#read in the first time step of S and θ
expname = "iter129_bulkformula"

filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
τdatafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"

nt = length(datafilelist)
tecco= 1992+1/24:1/12:2018 # ecco years

tt = 200
@time EXF = γ.read(diagpath[expname]*τdatafilelist[tt],MeshArray(γ,Float32,15))
τx = EXF[:, 14]; τy = EXF[:, 15]; 
τE, τN = UVtoUEVN(τx, τy, Γ)
τE2, τN2 = rotate_uv(τx, τy, Γ)

proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

fig, axes = plt.subplots(nrows= 1, ncols = 3, figsize=(15,10), 
                        subplot_kw=Dict("projection"=> proj0))

vars = [τx, τE, τE2]; 
var_labels = ["Unrotated", "Gael's Rotation", "Jake's Rotation"]

nanmax(ma) = maximum([nm.maximum(ma[ff]) for ff = 1:5])
bounds = maximum([nanmax(abs.(var)) for var in vars])

for (i, ax) in enumerate(axes)
    cff = Any[]
    for ff = 1:5
        cf = ax.pcolormesh(λ_offset[ff], ϕ[ff],  vars[i][ff], 
                    vmin = -bounds, vmax = bounds, shading="nearest", 
                    transform=projPC, rasterized = true, cmap = cm.curl)
        push!(cff, cf)
    end

    push!(cfs, cff[1])
    gl = ax.gridlines(crs=projPC, draw_labels=true, linewidth=2, 
                  color="gray", alpha=0, linestyle="--")
    gl.top_labels = false; gl.right_labels = false
    # ax.coastlines()
    ax.set_title(var_labels[i])
end
fig.colorbar(cf, ax = axes[:], orientation = "horizontal", 
fraction=0.04, pad = 0.1, label = "ρ")
# fig.tight_layout()
fig1 = fig; fig1

#################### Doing Vertical Velocity
fig, axes = plt.subplots(nrows= 1, ncols = 3, figsize=(15,10), 
                        subplot_kw=Dict("projection"=> proj0))

vars = [τy, τN, τN2]; 
var_labels = ["Unrotated", "Gael's Rotation", "Jake's Rotation"]

bounds = maximum([nanmax(abs.(var)) for var in vars])


for (i, ax) in enumerate(axes)
    cff = Any[]
    for ff = 1:5
        cf = ax.pcolormesh(λ_offset[ff], ϕ[ff],  vars[i][ff], 
                    vmin = -bounds, vmax = bounds, shading="nearest", 
                    transform=projPC, rasterized = true, cmap = cm.curl)
        push!(cff, cf)
    end

    push!(cfs, cff[1])
    gl = ax.gridlines(crs=projPC, draw_labels=true, linewidth=2, 
                  color="gray", alpha=0, linestyle="--")
    gl.top_labels = false; gl.right_labels = false
    # ax.coastlines()
    ax.set_title(var_labels[i])
end
fig.colorbar(cf, ax = axes[:], orientation = "horizontal", 
fraction=0.04, pad = 0.1, extend = "both", label = "ρ")
# fig.tight_layout()
fig2 = fig; fig2
