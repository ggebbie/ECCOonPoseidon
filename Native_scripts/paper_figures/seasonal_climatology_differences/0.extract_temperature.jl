include("../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools
import PyPlot as plt

include(srcdir("plot_and_dir_config.jl"))
@pyimport matplotlib.patches as patches

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)
region = "NPAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)

cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)
ΔV = lateral_sum(cell_volumes)

lvls = findall( -3000 .<= -z[:].<= -2000)
mid_depths(x) = vec(sum(Float32, x[lvls, :] .* ΔV[lvls], dims = 1) / sum(Float32, ΔV[lvls]))

tecco = 1992+1/24:1/12:2018; nz = 50
get_datafiles(expname, key) = filter(x -> occursin("data",x),searchdir(diagpath[expname],key) )

function get_temperature(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid, cell_volumes)

    datafilelist_θ  = get_datafiles(expname, "state_3d_set1")
    ΔV = lateral_sum(cell_volumes)
    nt = length(datafilelist_θ); nz = 50
    println(nt, " months available")
    θ_avg = zeros(Float32, nz, nt)
    ma_template = MeshArray(γ,Float32,50)
    @time for tt = 1:nt
        println(tt)
        fnameθ = datafilelist_θ[tt]
        @time θ = γ.read(diagpath[expname]*fnameθ,ma_template)

        θ_avg[:, tt] .= lateral_sum(θ .* cell_volumes); 
        θ_avg[:, tt] .= θ_avg[:, tt] ./ ΔV
    end

    return θ_avg
end

diagpath["mean_tau_adjusts"] = vastdiagdir("seasonalclimatology", "run_only_clim_tau_adjusts")
diagpath["mean_tau_noadjusts"] = vastdiagdir("seasonalclimatology_iter0", "run")

adjust_exps = Dict()
adjust_exps["iter0_bulkformula"] = get_temperature(diagpath, "iter0_bulkformula", γ, cell_volumes)
adjust_exps["only_wind"] = get_temperature(diagpath, "only_wind", γ, cell_volumes)
adjust_exps["mean_tau_adjusts"] = get_temperature(diagpath, "mean_tau_adjusts", γ, cell_volumes)
adjust_exps["mean_tau_noadjusts"] = get_temperature(diagpath, "mean_tau_noadjusts", γ, cell_volumes)

jldsave(datadir(region * "_temperature_sens_exps_seasonal_differences.jld2"), adjust_exps= adjust_exps)