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

adjust_exps_NPAC = jldopen(datadir("NPAC" * "_temperature_sens_exps_seasonal_differences.jld2"))["adjust_exps"]
adjust_exps_NPAC_anomaly = Dict()
adjust_exps_NPAC_anomaly = Dict(key => (adjust_exps_NPAC[key] .- adjust_exps_NPAC[key][:, 1]) for key in keys(adjust_exps_NPAC))

adjust_exps_NPAC_test = jldopen(datadir("NPAC" * "_temperature_sens_seasonal_test.jld2"))["adjust_exps"]
adjust_exps_NPAC_test_anomaly = Dict()
adjust_exps_NPAC_test_anomaly = Dict(key => (adjust_exps_NPAC_test[key] .- adjust_exps_NPAC_test[key][:, 1]) for key in keys(adjust_exps_NPAC_test))

fig, ax = plt.subplots(figsize = (10, 10))
ax.plot(mid_depths(adjust_exps_NPAC["iter0_bulkformula"]), label = "CTRL")
ax.plot(mid_depths(adjust_exps_NPAC["only_wind"]), label = "only wind")
ax.plot(mid_depths(adjust_exps_NPAC_test["mean_tau_noadjusts"]), label = "CTRL Seasonal Tau")
ax.plot(mid_depths(adjust_exps_NPAC_test["mean_tau_adjusts"]), label = "129 Seasonal Tau")
ax.plot(mid_depths(adjust_exps_NPAC_test["mean_tau_noadjust_redo"]), label = "CTRL Seasonal Tau (redo1)")
ax.plot(mid_depths(adjust_exps_NPAC_test["mean_tau_noadjust_redo_bf"]), label = "CTRL Seasonal Tau (bf)")
ax.plot(mid_depths(adjust_exps_NPAC_test["mean_sfc_noadjusts"]), label = "CTRL Seasonal SFC")


ax.legend()
fig

E, F = trend_matrices(tecco)
F * mid_depths(adjust_exps_NPAC_test_anomaly["mean_tau_noadjust_redo"])
F * mid_depths(adjust_exps_NPAC_test_anomaly["mean_tau_noadjusts"])
F * mid_depths(adjust_exps_NPAC_test_anomaly["mean_sfc_noadjusts"])

F * mid_depths(adjust_exps_NPAC_anomaly["iter0_bulkformula"])
F * mid_depths(adjust_exps_NPAC_anomaly["only_wind"])

nt = length(mid_depths(adjust_exps_NPAC_test["mean_tau_noadjust_redo"]))

fig, ax = plt.subplots(figsize = (10, 10))
ax.plot(1:312, mid_depths(adjust_exps_NPAC_test["mean_tau_noadjusts_bf"])[1:312], label = "CTRL Seasonal Tau BF")
ax.plot(1:nt, mid_depths(adjust_exps_NPAC_test["mean_tau_yesadjusts_bf"]), label = "FULL Seasonal Tau BF")
ax.plot(1:312, mid_depths(adjust_exps_NPAC_test["mean_tau_adjusts"])[1:312], label = "FULL Seasonal Tau")
ax.plot(1:312, mid_depths(adjust_exps_NPAC_test["mean_tau_noadjusts"])[1:312], label = "CTRL Seasonal Tau")
ax.plot(1:312, mid_depths(adjust_exps_NPAC_test["mean_sfc_noadjusts"])[1:312], label = "CTRL Seasonal SFC")

ax.plot(1:312, mid_depths(adjust_exps_NPAC["iter0_bulkformula"])[1:312], label = "iter0")
ax.legend()
fig


fig, ax = plt.subplots(figsize = (10, 10))
ax.plot(tecco, mid_depths(adjust_exps_NPAC_anomaly["iter0_bulkformula"]), label = "CTRL")
ax.plot(tecco, mid_depths(adjust_exps_NPAC_test_anomaly["only_init"]), label = "ONLY INIT")
ax.plot(tecco[1:312], mid_depths(adjust_exps_NPAC_test_anomaly["mean_tau_noadjusts"])[1:312], label = "CTRL Seasonal Tau")
ax.plot(tecco[1:nt], mid_depths(adjust_exps_NPAC_test_anomaly["mean_tau_noadjusts_bf"])[1:nt], label = "CTRL Seasonal Tau BF")
ax.plot(tecco[1:312], mid_depths(adjust_exps_NPAC_test_anomaly["mean_tau_adjusts"])[1:312], label = "FULL Seasonal Tau")

ax.legend()
fig