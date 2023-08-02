# this script will compile and save a basinwide average of OHC 
# for a particular basin (basin_name) 

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise 
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools,
    PyPlot, PyCall, JLD2, DrWatson, Statistics, JLD2
using .OHC_helper
include(srcdir("config_exp.jl"))
ccrs = pyimport("cartopy.crs")

do_constant_density=true 
do_total=true 
output2file = true
ocean_mask = Γ.hFacC[:,1] 
ocean_mask[findall(ocean_mask.>0.0)].=1.0
ocean_mask[findall(ocean_mask.<=0.0)].=NaN

(ϕ,λ) = latlonC(γ) #each tile has dims=(λ, ϕ)
area = readarea(γ)
for ff ∈ [3,4] 
    λ[ff][λ[ff] .<= 0 ] = λ[ff][λ[ff] .<= 0 ] .+ 360
end
runpath,diagpath = listexperiments(exprootdir())
shortnames = expnames()
marks = expsymbols()
nexp = length(shortnames) # number of experiments
fileroot = "state_3d_set1"
filedir = "ECCO_vars/"
filename = "THETAs"

i = 0 
#animation tutorial found here https://genkuroki.github.io/documents/Jupyter/20170624%20Examples%20of%20animations%20in%20Julia%20by%20PyPlot%20and%20matplotlib.animation.html
for (keys,values) in shortnames
    expname = keys
    println(expname)
    lvl = 1
    time = 300
    @load datadir(filedir * filename*"_"*expname*".jld2") var_exp
    
    fig, ax = plt.subplots(1, 1, figsize=(5,5), subplot_kw=Dict("projection"=> ccrs.PlateCarree(central_longitude = -150)) )
    ax.set_extent((110, -70, -55, 60))
    # ax.coastlines()
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=true,
                         linewidth=2, color="gray", alpha=0, linestyle="--")
    var_sel = var_exp[time][:, lvl].*ocean_mask
    cf = Vector{Any}(undef ,1)

    # cf[1] = ax.contourf(λ[2], ϕ[2], var_sel[2]/[2],transform=ccrs.PlateCarree())=
    min_val = min_mesharray(var_sel)
    max_val = max_mesharray(var_sel)
    for ff in 1:length(var_sel)
            #pcolormesh can handle nans, contourf cannot
            cf[1] = ax.pcolormesh(λ[ff], ϕ[ff], var_sel[ff],shading="nearest", 
            transform=ccrs.PlateCarree(), rasterized = true, vmin = min_val, vmax = max_val)
    end
    cbar = fig.colorbar(cf[1])
    fig.tight_layout()
    folder = "/home/ameza/ECCOonPoseidon/plots/OHC Climatologies/"
    outputfile = plotsdir(folder * "climatology_" * expname * ".pdf")
    println("saving plots at.. " * outputfile)
    fig.savefig(outputfile)
    close("all")
end

