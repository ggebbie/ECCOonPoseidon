include("../../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools
import DataFrames as DF
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

tecco = 1992+1/24:1/12:2018; nz = 50
get_datafiles(expname, key) = filter(x -> occursin("data",x),searchdir(diagpath[expname],key) )

adjust_exps_PAC = jldopen(datadir("θ_abs_DIV_PAC_all_times.jld2"))["adjust_exps"]
adjust_exps_NPAC = jldopen(datadir("θ_abs_DIV_NPAC_all_times.jld2"))["adjust_exps"]
adjust_exps_GLOB = jldopen(datadir("θ_ABS_DIV_GLOBAL_all_times.jld2"))["adjust_exps"]


fig, ax = plt.subplots(1, 3, figsize = (22.5, 7.5))
ax[1].plot(tecco, 1e4 * adjust_exps_NPAC["only_init"], label = "Different Initial Condition \n Iteration 0 Forcing")
ax[1].plot(tecco, 1e4 * adjust_exps_NPAC["noinitadjust"], label = "Different Initial Condition \n Iteration 129 Forcing")
ax[1].set_ylabel("Mean Squared Difference [cK]²")
ax[1].set_xlabel("time")
ax[1].set_title("N. PAC Mean Squared \n Temperature Difference\n z = 2000 - 3000")

ax[2].plot(tecco, 1e4 * adjust_exps_PAC["only_init"], label = "Different Initial Condition \n Iteration 0 Forcing")
ax[2].plot(tecco, 1e4 * adjust_exps_PAC["noinitadjust"], label = "Different Initial Condition \n Iteration 129 Forcing")
ax[2].set_ylabel("Mean Squared Difference [cK]²")
ax[2].set_xlabel("time")
ax[2].set_title("PAC Mean Squared \n Temperature Difference\n z = 2000 - 3000")
ax[2].legend(frameon = false, loc = (-0.3, -0.3), ncols = 2)

ax[3].plot(tecco, 1e4 * adjust_exps_GLOB["only_init"], label = "Different Initial Condition \n Iteration 0 Forcing")
ax[3].plot(tecco, 1e4 * adjust_exps_GLOB["noinitadjust"], label = "Different Initial Condition \n Iteration 129 Forcing")
ax[3].set_ylabel("Mean Squared Difference [cK]²")
ax[3].set_xlabel("time")
ax[3].set_title("Global Mean Squared \n Temperature Difference\n z = 2000 - 3000")
fig.savefig(plotsdir("native/paper_figures/diffs_initial_cond.png"), bbox_inches = "tight", dpi = 400)
fig


fig, ax = plt.subplots(1, 3, figsize = (22.5, 7.5))
ax[1].plot(tecco, 1e4 * adjust_exps_NPAC["only_sfc"], label = "Iteration 0 Initial Condition \n Different Forcings")
ax[1].plot(tecco, 1e4 * adjust_exps_NPAC["nosfcadjust"], label = "Iteration 129 Initial Condition \n Different Forcings")
ax[1].set_ylabel("Mean Squared Difference [cK]²")
ax[1].set_xlabel("time")
ax[1].set_title("N. PAC Mean Squared \n Temperature Difference\n z = 2000 - 3000")

ax[2].plot(tecco, 1e4 * adjust_exps_PAC["only_sfc"], label = "Iteration 0 Initial Condition \n Different Forcings")
ax[2].plot(tecco, 1e4 * adjust_exps_PAC["nosfcadjust"], label = "Iteration 129 Initial Condition \n Different Forcings")
ax[2].set_ylabel("Mean Squared Difference [cK]²")
ax[2].set_xlabel("time")
ax[2].set_title("PAC Mean Squared \n Temperature Difference\n z = 2000 - 3000")
ax[2].legend(frameon = false, loc = (-0.3, -0.3), ncols = 2)

ax[3].plot(tecco, 1e4 * adjust_exps_GLOB["only_sfc"], label = "Iteration 0 Initial Condition \n Different Forcings")
ax[3].plot(tecco, 1e4 * adjust_exps_GLOB["nosfcadjust"], label = "Iteration 129 Initial Condition \n Different Forcings")
ax[3].set_ylabel("Mean Squared Difference [cK]²")
ax[3].set_xlabel("time")
ax[3].set_title("Global Mean Squared \n Temperature Difference\n z = 2000 - 3000")
fig.savefig(plotsdir("native/paper_figures/diffs_forcing.png"), bbox_inches = "tight", dpi = 400)

fig


fig, ax = plt.subplots(1, 3, figsize = (22.5, 7.5))
ax[1].plot(tecco, 1e4 * adjust_exps_NPAC["only_kappa"], label = "Iteration 0 Initial Condition and Forcing \n Different Mixing Coefficient")
# ax[1].plot(tecco, 1e4 * adjust_exps_NPAC["nosfcadjust"], label = "Iteration 129 Initial Condition \n Different Forcings")
ax[1].set_ylabel("Mean Squared Difference [cK]²")
ax[1].set_xlabel("time")
ax[1].set_title("N. PAC Mean Squared \n Temperature Difference\n z = 2000 - 3000")

ax[2].plot(tecco, 1e4 * adjust_exps_PAC["only_kappa"], label = "Iteration 0 Initial Condition and Forcing \n Different Mixing Coefficient")
# ax[2].plot(tecco, 1e4 * adjust_exps_PAC["nosfcadjust"], label = "Iteration 129 Initial Condition \n Different Forcings")
ax[2].set_ylabel("Mean Squared Difference [cK]²")
ax[2].set_xlabel("time")
ax[2].set_title("PAC Mean Squared \n Temperature Difference\n z = 2000 - 3000")
ax[2].legend(frameon = false, loc = (-0.3, -0.3), ncols = 2)

ax[3].plot(tecco, 1e4 * adjust_exps_GLOB["only_kappa"], label = "Iteration 0 Initial Condition and Forcing \n Different Mixing Coefficient")
# ax[3].plot(tecco, 1e4 * adjust_exps_GLOB["nosfcadjust"], label = "Iteration 129 Initial Condition \n Different Forcings")
ax[3].set_ylabel("Mean Squared Difference [cK]²")
ax[3].set_xlabel("time")
ax[3].set_title("Global Mean Squared \n Temperature Difference\n z = 2000 - 3000")
fig.savefig(plotsdir("native/paper_figures/diffs_mixing.png"), bbox_inches = "tight", dpi = 400)

fig
