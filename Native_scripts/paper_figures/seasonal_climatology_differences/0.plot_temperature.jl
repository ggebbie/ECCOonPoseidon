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
adjust_exps_NPAC_anomaly = Dict(key => (adjust_exps_NPAC[key] .- adjust_exps_NPAC[key][:, 1]) for key in keys(adjust_exps))


fig, ax = plt.subplots(figsize = (10, 10))
ax.plot(mid_depths(adjust_exps_NPAC["iter0_bulkformula"]), label = "CTRL")
ax.plot(mid_depths(adjust_exps_NPAC["iter129_bulkformula"]), label = "FULL")
ax.plot(mid_depths(adjust_exps_NPAC["noinitadjust"]), label = "NO INIT")
ax.plot(mid_depths(adjust_exps_NPAC["only_init"]), label = "ONLY INIT")
ax.legend()
fig

fig, ax = plt.subplots(figsize = (10, 10))
ax.plot(tecco, mid_depths(adjust_exps_NPAC_anomaly["iter0_bulkformula"]), label = "CTRL")
ax.plot(tecco, mid_depths(adjust_exps_NPAC_anomaly["iter129_bulkformula"]), label = "FULL")
ax.plot(tecco, mid_depths(adjust_exps_NPAC_anomaly["noinitadjust"]), label = "NO INIT")
ax.plot(tecco, mid_depths(adjust_exps_NPAC_anomaly["only_init"]), label = "ONLY INIT")
ax.legend()
fig

fig, ax = plt.subplots(figsize = (10, 10))
ax.plot(tecco, 100 .* mid_depths(adjust_exps_NPAC_anomaly["only_init"] .- adjust_exps_NPAC_anomaly["iter0_bulkformula"]), 
        label = "Linearizing about Iteration 0")
ax.plot(tecco, 100 .* mid_depths(adjust_exps_NPAC_anomaly["iter129_bulkformula"] .- adjust_exps_NPAC_anomaly["noinitadjust"]), 
        label = "Linearizing about Iteration 129")
ax.legend()
ax.set_ylabel("[cK]")
ax.set_xlabel("time")

x1 = 100 .* mid_depths(adjust_exps_NPAC_anomaly["only_init"] .- adjust_exps_NPAC_anomaly["iter0_bulkformula"])
x2 = 100 .* mid_depths(adjust_exps_NPAC_anomaly["iter129_bulkformula"] .- adjust_exps_NPAC_anomaly["noinitadjust"])
maximum(abs.(x1 .- x2))

ax.set_title("Effect of Initial Condition Control Adjustments \n on North Pacific Temperature Anomaly")
fig