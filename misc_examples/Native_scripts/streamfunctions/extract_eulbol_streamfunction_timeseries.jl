include("../../src/intro.jl")

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

vars =  ["only_init", "only_kappa", "only_sfc", "iter129_bulkformula",  "iter0_bulkformula"]
vars =  ["only_buoyancy", "only_wind"]

@time for expname in vars
    println(expname)
    @time Ψ_exp_timeseries, ϕ_avg = extract_meridional_Ψ(expname,diagpath, Γ, γ, PAC_msk)
    jldsave(datadir("Ψ_EulBol_timeseries_"*region*"_" * expname *".jld2"); 
    Ψ_exp_timeseries, ϕ_avg)
end