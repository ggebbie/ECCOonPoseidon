using Pkg
Pkg.activate(".")
include("../../../src/intro.jl")

pwd()
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
diagpath["mean_tau_noadjusts"] = vastdiagdir("seasonalclimatology", "run_only_clim_iter0_tau2")
diagpath["mean_tau_noadjust_redo"] = vastdiagdir("seasonalclimatology", "run_only_clim_no_tau_adjusts")
diagpath["mean_tau_noadjust_redo_bf"] = vastdiagdir("seasonalclimatology", "run_only_clim_iter0_tau_bf")
# diagpath["mean_tau_yesadjust_redo_bf_old"] = vastdiagdir("seasonalclimatology", "run_only_clim_iter129_tau_bf")
diagpath["mean_tau_yesadjust_redo_bf"] = vastdiagdir("seasonalclimatology", "run_only_clim_iter129_tau_bf2")


adjust_exps = Dict()
adjust_exps["mean_tau_noadjust_redo_bf"] = get_temperature(diagpath, "mean_tau_noadjust_redo_bf", γ, cell_volumes)
adjust_exps["mean_tau_yesadjust_redo_bf"] = get_temperature(diagpath, "mean_tau_yesadjust_redo_bf", γ, cell_volumes)
adjust_exps["mean_tau_adjusts"] = get_temperature(diagpath, "mean_tau_adjusts", γ, cell_volumes)
adjust_exps["iter0_bulkformula"] = get_temperature(diagpath, "iter0_bulkformula", γ, cell_volumes)

nt = size(adjust_exps["mean_tau_yesadjust_redo_bf"], 2)
fig, ax = plt.subplots(figsize = (10, 5))
ax.plot(1:nt, mid_depths(adjust_exps["mean_tau_noadjust_redo_bf"])[1:nt], label = "noadjusts")
ax.plot(1:nt, mid_depths(adjust_exps["mean_tau_adjusts"])[1:nt], label = "ff_129")
ax.plot(1:nt, mid_depths(adjust_exps["iter0_bulkformula"])[1:nt], label = "iter0")
ax.plot(1:nt, mid_depths(adjust_exps["mean_tau_yesadjust_redo_bf"])[1:nt], label = "bf_129")
# ax.scatter(1:nt, mid_depths(adjust_exps["mean_tau_noadjust_redo_bf"])[1:nt], label = "noadjustsredobf", c = "k", s = 5)
ax.set_xlim(0, nt+10)
ax.legend()
fig
# nt = size(adjust_exps["mean_tau_noadjust_redo_bf"], 2)
# adjust_exps["mean_tau_noadjusts"][38:42, 1:nt] .- adjust_exps["mean_sfc_noadjusts"][38:42, :]
# adjust_exps["mean_sfc_noadjusts"][38:42, :] .- adjust_exps["mean_tau_noadjust_redo"][38:42, 1:nt] 
# adjust_exps["mean_tau_noadjusts"][38:42, 1:nt] .- adjust_exps["mean_tau_noadjust_redo"][38:42, 1:nt]

jldsave(datadir(region * "_temperature_sens_seasonal.jld2"), adjust_exps= adjust_exps)

z[38:43]