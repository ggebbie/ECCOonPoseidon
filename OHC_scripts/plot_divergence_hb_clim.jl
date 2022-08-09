using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools,
PyPlot, JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf, RollingFunctions
using PyPlot   # important!
using PyCall
using DataFrames, GLM
import Plots: Animation, frame, gif, cgrad
import Plots: contourf as jcontourf
@pyimport imageio
@pyimport matplotlib.animation as animation

θ_flat = Dict()
for (key,values) in shortnames
    expname = key; println(key)
    @time sθ = load(datadir(filedir*fname1*"_"*expname * suffix * ".jld2"), "var_exp")
    θ_flat[expname] = Vector{MeshArrays.gcmarray{Float64, 1, Matrix{Float64}}}(undef, length(tecco))
    println("constructing time series")
    @time for tt in 1:length(tecco)
        θ_flat[expname][tt] = OHC_helper.smush(sθ[tt] .* crop_vols) ./ OHC_helper.smush(crop_vols)
    end
    @time GC.gc(true) #garbage collecting 
end
for ff ∈ [3,4] 
    λ[ff][λ[ff] .<= 0 ] = λ[ff][λ[ff] .<= 0 ] .+ 360
end

proj = ECCOonPoseidon.cartopy.crs.PlateCarree()
proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
filenames = [plotsdir() * "/OHC_Divergence/temp"*string(tt)*".png"
             for tt in 1:length(tecco)]

nanmaxi(a) = maximum(filter(!isnan,a)) 
nanmini(a) = minimum(filter(!isnan,a)) 

function nanmin(ma::MeshArrays.gcmarray)
    temp = [Inf]
    temp_ma = MeshArrays.mask(ma, Inf)
    for ff in 1:5
        (nanmini(temp_ma[ff][:]) < temp[1]) && (temp .= nanmini(temp_ma[ff][:]))
    end
    return temp[1]
end

function nanmax(ma::MeshArrays.gcmarray)
    temp = [-Inf]
    temp_ma = MeshArrays.mask(ma, -Inf)
    for ff in 1:5
        (nanmaxi(temp_ma[ff][:]) > temp[1]) && (temp .= nanmaxi(temp_ma[ff][:]))
    end
    return temp[1]
end

tt = 1
cf = Vector{Any}(undef ,1)
fig, ax = plt.subplots(1, 1, figsize=(7,5), subplot_kw=Dict("projection"=> proj0))
ax.set_title("Average Temperature in ECCO, " * region)
# θ_flat["iter129_bulkformula"][tt][findall(θ_flat["iter129_bulkformula"][tt] .== 0)] .= NaN
b1 = nanmin(θ_flat["iter129_bulkformula"][tt])
b2 = nanmax(θ_flat["iter129_bulkformula"][tt])
ax.set_global()
ax.coastlines()
ax.set_extent((110, -70, -0, 60))
ax.gridlines(crs=proj, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
for ff in 1:5
        cf[1] = ax.pcolormesh(λ[ff], ϕ[ff],  θ_flat["iter129_bulkformula"][tt][ff],
        vmin = b1, vmax = b2, shading="nearest", transform=projPC, rasterized = true, cmap = colorway)            
end

cbar = fig.colorbar(cf[1], label = L"^\circ"*"C",
orientation = "horizontal")
tight_layout()

fig.savefig(plotsdir() * "/OHC_Divergence/" * "θClim_" * region * suffix * ".png",
dpi = 1500)

close("all")