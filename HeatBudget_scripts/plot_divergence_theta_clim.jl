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

import NaNMath as nm
nanmin(ma::MeshArrays.gcmarray) = minimum(nm.minimum.(ma.f))
nanmax(ma::MeshArrays.gcmarray) = maximum(nm.maximum.(ma.f))

θ_flat = Dict()
λ360 = deepcopy(λ)
expname = key; println(key)
θ_flat[expname] = Vector{MeshArrays.gcmarray{Float64, 1, Matrix{Float64}}}(undef, length(tecco))
println("constructing time series")
@time for tt in 1:length(tecco)
    fdir1 = filedir * expname * "/" * fname1 * "_" * suffix * "_"*"$tt" *"_.jld2"
    sθ = load_object(datadir(fdir1))
    θ_flat[expname][tt] = OHC_helper.smush(sθ .* crop_vols) ./ OHC_helper.smush(crop_vols)
end
@time GC.gc(true) #garbage collecting 

for ff ∈ [3,4] 
    λ360[ff][λ[ff] .<= 0 ] = λ360[ff][λ[ff] .<= 0 ] .+ 360
end


proj = ECCOonPoseidon.cartopy.crs.PlateCarree()
proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

filenames = [plotsdir() * "/OHC_Divergence/temp"*string(tt)*".png"
             for tt in 1:length(tecco)]

tt = 312
fig, ax = plt.subplots(2, 2, figsize=(12,5), subplot_kw=Dict("projection"=> proj0))
cf = Vector{Any}(undef ,1)
fig.suptitle("Δθ in ECCO")

for (i, expname) in enumerate(collect(keys(θ_flat)))
    ax[i].set_title(region)
    # θ_flat["iter129_bulkformula"][tt][findall(θ_flat["iter129_bulkformula"][tt] .== 0)] .= NaN
    var = θ_flat[expname]
    baseline = θ_flat["iter0_bulkformula"][tt] 
    b1 = nanmin(θ_flat[expname][tt] .- baseline)
    b2 = nanmax(θ_flat[expname][tt] .- baseline)
    ax[i].set_global()
    ax[i].coastlines()
    ax[i].set_extent((110, -70, -0, 60))
    ax[i].gridlines(crs=proj, draw_labels=true,
                        linewidth=2, color="gray", alpha=0, linestyle="--")
    for ff in 1:5
            cf[1] = ax[i].pcolormesh(λ[ff], ϕ[ff],  θ_flat[expname][tt][ff],
            vmin = b1, vmax = b2, shading="nearest", transform=projPC, rasterized = true, cmap = colorway)            
    end

end
cbar = fig.colorbar(cf[1], label = L"^\circ"*"C",
orientation = "horizontal")
fig.tight_layout()
fig.savefig(plotsdir() * "/OHC_Divergence/Clims/" * "θClimAnom_$tt" * region * suffix *".png",
dpi = 1500)
close("all")