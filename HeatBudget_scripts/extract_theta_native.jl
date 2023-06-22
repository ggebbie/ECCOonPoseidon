#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings,
Plots
gr()
import NaNMath as nm
using .OHC_helper
import PyPlot as plt
using PyCall

@pyimport seaborn as sns;
@pyimport pandas as pd;
colors =  sns.color_palette("deep")[1:4]
sns.set_theme(context = "talk", font_scale = 1.0,
              palette = sns.color_palette("deep"));#sns.set_context("talk")
pygui(false)
cm = pyimport("cmocean.cm");colorway = cm.balance;
include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
area = readarea(γ)
runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
labels_L = [L"\theta^{\Delta F, \Delta T}", L"\theta^{\Delta F}",L"\theta^{\Delta T}", L"\theta^{0}"]
tecco = 1992+1/24:1/12:2018

ocean_mask = OHC_helper.wet_pts(Γ)
region = "PAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
                            region, include_bering = false)
msk = PAC_msk;
(msk == ocean_mask) && (region = "GLOB")
suffix = "sfctobot"
uplvl = Inf; botlvl = -Inf
lvls = findall( botlvl .<= z[:].<= uplvl)

cell_depths = get_cell_depths(msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths);
smush_depths = smush(cell_depths); smush_depths[findall(smush_depths .== 0)] = Inf
inv_depths = 1 ./ smush_depths

function filter_θ(γ::gcmgrid, diagpath::Dict{String, String}, expname::String, 
    θz::Dict, θ_depths::Dict, θ_zonal::Dict, θsurf::Dict, 
    lvls::Vector{Int64}, 
    mask ::MeshArray,
    inv_depths::MeshArrays.gcmarray{T, 1, Matrix{T}},
    cell_volumes::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T <: Real

    vol_weight(x::MeshArray) = Float32(sum(x .* mask) / tot_vol)

    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
    datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    nt = length(datafilelist_θ)
    nlev = length(lvls)

    θsurf[expname] = zeros(Float32, nt)
    θz[expname] = zeros(Float32, nt)
    θ_depths[expname] = zeros(Float32, nlev, nt)
    θ_zonal[expname] = zeros(nlev, 270, nt)

    Threads.@threads for tt = 1:nt
        fnameS = datafilelist_S[tt]
        fnameθ = datafilelist_θ[tt]
        sθ = extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_depths)
        sθ_mask = sθ[:, lvls] .* mask #apply mask and crop 

        θsurf[expname][tt] = Float32.(volume_mean(sθ[:, 1] .* mask; weights = cell_volumes[:, 1]))
        θz[expname][tt] = Float32.(volume_mean(sθ_mask; weights = cell_volumes[:, lvls]))
        θ_depths[expname][:, tt] .= OHC_helper.ma_horiz_avg(sθ_mask, cell_volumes[:, lvls])
        θ_zonal[expname][:, :, tt] .= OHC_helper.ma_zonal_avg(sθ_mask, cell_volumes[:, lvls])
        GC.safepoint()
    end
    θ_zonal[expname][θ_zonal[expname] .== 0 ] .= NaN
end

θz, θ_depths, θ_zonal = Dict(), Dict(), Dict()
θsurf = Dict()
@time for expname in keys(shortnames)
    println(expname)
    filter_θ(γ, diagpath, expname, 
             θz, θ_depths, θ_zonal, θsurf, 
             lvls, msk, inv_depths, cell_volumes)
end

# jldsave(datadir("θzonal"*region*"_" * suffix * ".jld2"); θ_zonal)

#=Tries to reconstruct iter129 from iter0, nosfc & no init 
using MLR
=#
@time include("plot_divergence_theta_LinReg.jl")

#constrained parameter search on convex domain
@time include("plot_divergence_theta_LinRegConvex.jl")


#=Plots temperature climatologies, could also plot changes
=#

@time include("plot_divergence_theta_timeseries.jl")

#=ceates a depth-latitude plot of temperature changes
with respect to iter0
=#
# @time include("plot_divergence_theta_verticalTS.jl")


# #creates a depth-latitude plot of temperature changes

# @time include("plot_divergence_theta_zonal.jl") 
@time include("plot_divergence_theta_zonal_mean.jl") 
@time include("plot_divergence_theta_zonal_derivative.jl") 