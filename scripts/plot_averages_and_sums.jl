
include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyPlot, PyCall
using .OHC_helper

include(srcdir("config_exp.jl"))
include(srcdir("MeshArraysPlots.jl"))

@pyimport cmocean.cm as cmo

(ϕ,λ) = latlonC(γ)

runpath,diagpath = listexperiments(exprootdir());

#testing vertical sum function
cell_depths = get_cell_thickness(ocean_mask, ΔzF, Γ.hFacC); 
Depths = vertical_sum(cell_depths)

proj0 = crs.PlateCarree(central_longitude=-150)
projPC = crs.PlateCarree()

fig, ax = plt.subplots(1, figsize=(15,20), 
subplot_kw=Dict("projection"=> proj0), sharey = true)
difference = Depths .- Γ.Depth #should be the same
CB = pcolormesh_ma(ax, λ, ϕ, difference, cmo.balance)
fig.colorbar(CB, ax=ax, orientation = "horizontal", fraction = 0.02)
println("Maximal Error: ", maximum(abs.(difference)), " meters")
ax.coastlines()
fig
