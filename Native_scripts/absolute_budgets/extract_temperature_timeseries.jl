#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, FLoops
using .OHC_helper
import PyPlot as plt 

include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());

# abbreviations for each experiment for labels, etc.
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
nexp = length(shortnames) # number of experiments
 
ocean_mask = OHC_helper.wet_pts(Γ)
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));

H = OHC_helper.smush(cell_depths, γ); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H
tecco = 1992+1/24:1/12:2018
nz = 50

ΔV = zeros(Float32, nz)
[ΔV[k] = Float32(sum(cell_volumes[:, k])) for k=1:nz]
lvls = findall( -3000 .<= z[:].<= -2000)

masked_volume = cell_volumes[:, lvls]
sum_masked_volume = Float32(sum(masked_volume))

function heat_flux_profile(ds::MeshArray, ΔV)
    nz = size(ds, 2)
    vol_avg = zeros(Float32, nz)

    for ff=1:5, k=1:nz
        vol_avg[k] += Float32(sum(ds[ff, k])) ./ ΔV[k]
    end
    return vol_avg
end

function filter_heat_budget(savename, diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid)

    filelist = searchdir(diagpath[expname],"state_3d_set1") 
    datafilelist_θ  = filter(x -> occursin("data",x),filelist)
    filelist = searchdir(diagpath[expname], "state_2d_set1") # first filter for state_3d_set1
    datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
 
    nt = length(datafilelist_θ); nz = 50
    vars = zeros(Float32, nz, nt)

    @time for tt = 1:nt
        println(tt)
        fnameθ = datafilelist_θ[tt]
        fnameS = datafilelist_S[tt]

        sθ = OHC_helper.extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_H)
        vars[:, tt]  = heat_flux_profile(sθ .* cell_volumes, ΔV); 
    end

    return vars
end

θ = Dict()
for expname in ["iter129_bulkformula", "iter0_bulkformula", "seasonalclimatology", "seasonalclimatology_iter0"]
    θ[expname] = filter_heat_budget(savename, diagpath, expname, γ)
end

svename = datadir("native/" * region * "_THETA_levels" * ".jld2")
jldsave(svename, θ = θ)

k = 40; z[k]
vm = maximum(abs.(θ["iter129_bulkformula"] - mean(θ["iter129_bulkformula"], dims = 2)))
fig, ax = plt.subplots(2)
for (i, expname) in enumerate(["iter129_bulkformula", "seasonalclimatology"])
    ax[i].contourf(tecco, z, θ[expname] .- mean(θ[expname], dims = 2), 
    label = expname)
end
ax.legend()
fig