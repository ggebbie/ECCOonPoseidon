# this script will compile and save a basinwide average of OHC 
# for a particular basin (basin_name) 

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise 
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools,
    PyPlot, JLD2, DrWatson, Statistics, JLD2
import CairoMakie as Mkie
import GeoMakie
using .OHC_helper
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

for ff in 1:length(area)
    above_SO = (ϕ[ff] .> -56.0) #removes southern ocean 
    basin_mask[ff] .= ocean_mask[ff].*(basins[ff].==basinID) .* above_SO
end

cell_depths = get_basin_depths(basin_mask, Δz, Γ.hFacC)
basin_volume = get_basin_volumes(area, cell_depths)
plot_patch(Γ, basin_mask, basin_name)

OHC_native = Dict{String,Array{Any,2}}() # don't forget trailing parentheses
tecco = 1992+1/24:1/12:2018
monthsperyear = 12

temp_grid = similar(basin_volume)
for (keys,values) in shortnames
    expname = keys
    println(expname)
    nz = length(Δz)
    filelist = searchdir(diagpath[expname],fileroot) # 1st filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # 2nd filter for "data"
    if do_constant_density
        @time OHC = calc_OHC(diagpath, expname,datafilelist,γ,basin_volume, 
                    temp_grid, nz)
    else
        nothing
    end
    OHC_native[expname] = OHC  #OHC returned in joules 
end

filedir = "OHC_data/"
filename = "full_OHC_"*basin_name*"_native_scale_.jld2"
@save datadir(filedir * filename) OHC_native

filename = basin_name*"LevelVolumes.jld2" 
level_volumes = zeros(length(z))
for ff in eachindex(basin_volume)
    level_volumes[ff[2]] += sum(basin_volume[ff])
end
@save datadir(filedir * filename) level_volumes

isdir(plotsdir()) ? nothing : mkdir(plotsdir())