include("../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings
import NaNMath as nm
using PyCall

include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

include(srcdir("plot_and_dir_config.jl"))
 
ocean_mask = wet_pts(Γ)
region = "PAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "false", include_bering = true)

tecco = 1992+1/24:1/12:2018

diagpath["mean_tau_adjusts"] = vastdiagdir("seasonalclimatology", "run_only_clim_tau_adjusts")
fileroot = "trsp_3d_set1"
filelist = searchdir(diagpath["mean_tau_adjusts"],fileroot) # first filter for state_3d_set1
datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
print("Available months for analysis: ",  100 * length(datafilelist) / 312, " %")

vars =  ["mean_tau_adjusts"]

@time for expname in vars
    println(expname)
    @time Ψ_exp_timeseries, ϕ_avg = extract_meridional_Ψ(expname,diagpath, Γ, γ, PAC_msk)
    jldsave(datadir("Ψ_EulBol_timeseries_"*region*"_" * expname *".jld2"); 
    Ψ_exp_timeseries, ϕ_avg)
end