#include("MeshArray_helper.jl")
include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise 
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools,
    PyPlot, JLD2, DrWatson, Statistics, JLD2
using .OHC_helper

import CairoMakie as Mkie
import GeoMakie
include(srcdir("config_exp.jl"))

do_constant_density=true 
do_total=true 
# basin_name="Pacific"
basin_name = "Atlantic"
output2file = true

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir())
# abbreviations for each experiment for labels, etc.
shortnames = expnames()
marks = expsymbols()
nexp = length(shortnames) # number of experiments
fileroot = "state_3d_set1"

# non-zero Γ.hFacC[:,1] indicate wet points
ocean_mask = Γ.hFacC[:,1] 
ocean_mask[findall(ocean_mask.>0.0)].=1.0
ocean_mask[findall(ocean_mask.<=0.0)].=0.0

basinID=findall(basin_list.==basin_name)[1]
basin_mask=similar(basins)
cell_depths = similar(Γ.hFacC) .* 0.0 

for ff in 1:length(area)
    #creating a mask for the selected basin (:basin_name)
    basin_mask[ff] .= ocean_mask[ff].*(basins[ff].==basinID) 
end

cell_depths = get_basin_depths(basin_mask, Δz, Γ.hFacC)
basin_volume = get_basin_volumes(area, cell_depths)

θ_native = Dict{String,Array{Any,2}}()
tecco = 1992+1/24:1/12:2018
monthsperyear = 12

for (keys,values) in shortnames
    expname = keys
    println(expname)
    nz = length(Δz)
    filelist = searchdir(diagpath[expname],fileroot) # 1st filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # 2nd filter for "data"
    if do_constant_density
        @time θ = calc_θ_bar(diagpath, expname,datafilelist,γ,basin_volume, 
                    nz, level_volumes)
    else
        nothing
    end
    θ_native[expname] = θ  
end

outputfile, output_tree = OHC_outputpath("full_θ_comp_"*basin_name*"_.pdf", 
                plotsdir(), do_constant_density)

θ_datadir = joinpath(datadir(), "θ_data")
output_tree= joinpath(θ_datadir, output_tree[length(plotsdir()) + 2:end])
isdir(output_tree) ? nothing : mkpath(output_tree)

@save joinpath(output_tree, "full_θ_comp_"*basin_name*"_native_scale_.jld2") θ_native
isdir(plotsdir()) ? nothing : mkdir(plotsdir())

level_volumes = zeros(length(z))
for ff in eachindex(basin_volume)
    level_volumes[ff[2]] += sum(basin_volume[ff])
end

@save joinpath(output_tree, basin_name*"LevelVolumes.jld2") level_volumes
