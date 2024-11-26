include("../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour, PyCall, 
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings
import NaNMath as nm
# import PyPlot as plt

include(srcdir("config_exp.jl"))
@pyimport seaborn as sns;

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

include(srcdir("plot_and_dir_config.jl"))

@pyimport cmocean.cm as cmos


tecco = collect(1992+1/24:1/12:2018)
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= -z[:].<= uplvl)

region = "NPAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)
ctrl_vol = sum(cell_volumes[:, lvls])

z_ref = findall( -2700 .<= -z[:].<= -2300)[1]

vars =  ["iter0_bulkformula", "iter129_bulkformula","only_init", "only_kappa", "only_wind", "only_buoyancy"]
ΨEul = Dict(); ΨBol = Dict()

wθ_resid = Dict(); wθ_eul = Dict(); wθ_bol = Dict()

for (i, expname) in enumerate(vars)
    read_file = jldopen(datadir("Ψ_EulBol_timeseries_PAC_" * expname *".jld2"))
    Ψ_EulBol = read_file["Ψ_exp_timeseries"]
    read_file = jldopen(datadir("Ψ_Eul_timeseries_PAC_" * expname *".jld2"))
    Ψ_Eul = read_file["Ψ_exp_timeseries"]
    Ψ_Bol = Ψ_EulBol .- Ψ_Eul

    ϕ_avg = read_file["ϕ_avg"]
    ϕ_ref = findall( 23 .<= ϕ_avg .<= 40)

    WEul[expname] = -Ψ_Eul[z_ref, ϕ_ref[1], :][:]
    WBol[expname] = -Ψ_Bol[z_ref, ϕ_ref[1], :][:]

    fname = datadir("native/" * expname * "NPAC_THETA_budget_Wθ_eul_bol_" * suffix * ".jld2")

    vars = jldopen(fname)["dθ"]

    wθ_resid[expname] =  vars["wθBot"] .- vars["wθTop"]
    wθ_bol[expname] = vars["wθBot_bol"] .- vars["wθTop_bol"]

    wθ_eul[expname] = wθ_resid[expname] .- wθ_bol[expname]

end

effect_exp =  ["only_wind", "only_kappa", "only_init"]

[WEul[expt] .-= WEul["iter0_bulkformula"] for expt in effect_exp]
[WBol[expt] .-= WBol["iter0_bulkformula"] for expt in effect_exp]

[wθ_eul[expt] .-= wθ_eul["iter0_bulkformula"] for expt in effect_exp]
[wθ_bol[expt] .-= wθ_bol["iter0_bulkformula"] for expt in effect_exp]

fig, ax = plt.subplots(1, 3, figsize = (20, 5), sharey = true)
vars =  ["only_wind",  "only_init", "only_kappa"]

for (i, expname) in enumerate(vars)

    ax[1].scatter(1e-6 .* (ΨEul[expname] .+ ΨBol[expname]), 
    wθ_eul[expname] .+ wθ_bol[expname], c = exp_colors[expname], 
    label = plot_labels_effects[expname], linewidth = 2.5)

    ax[2].scatter(1e-6 .* ΨEul[expname], wθ_eul[expname], 
    c = exp_colors[expname], label = plot_labels_effects[expname], linewidth = 2.5)
    
    ax[3].scatter(1e-6 .* ΨBol[expname], wθ_bol[expname], 
    c = exp_colors[expname], label = plot_labels_effects[expname], linewidth = 2.5)

end

ax[1].set_ylabel("Vertical \n Convergence Anomaly Rate \n [cK per year]")
ax[1].set_xlabel("Streamfunction Anomaly [Sv]")
ax[1].set_title("Residual")
ax[1].legend(frameon = false)

# ax[2].set_ylabel("Vertical \n Convergence Anomaly Rate \n [cK per year]")
ax[2].set_xlabel("Streamfunction Anomaly [Sv]")
ax[2].set_title("Eulerian")
ax[2].legend(frameon = false)

# ax[3].set_ylabel("Vertical \n Convergence Anomaly Rate \n [cK per year]")
ax[3].set_xlabel("Streamfunction Anomaly [Sv]")
ax[3].set_title("Bolus")
ax[3].legend(frameon = false)

# ax[1].set_title(" Reconstruction")
# fig.savefig(plotsdir("native/sensitivity_exps/ΨEul_Effect.png"))
fig



fig, ax = plt.subplots(1, 2, figsize = (20, 7.5), sharey = false)
vars =  ["only_wind",  "only_init", "only_kappa"]
mult = 100 * 86400 * 365 * 10 
for (i, expname) in enumerate(vars)

    ax[1].plot(tecco, mult .* (wθ_eul[expname] .+ wθ_bol[expname]), c = exp_colors[expname], 
    label = plot_labels_effects[expname], linewidth = 2.5)

    ax[2].plot(tecco, 1e-6 .* ΨEul[expname], c = exp_colors[expname], 
    label = plot_labels_effects[expname], linewidth = 2.5)
end
[a.set_ylim(-2.5, 2.5) for a in ax]
ax[1].set_ylabel("[cK per decade]")
ax[2].set_ylabel("[Sv]")

ax[1].set_title("Vertical Heat Advection Convergence")
ax[2].set_title("Vertical Transport [z = 2500]")
ax[2].invert_yaxis()
[a.legend(frameon = false) for a in ax]
fig
fig.savefig(plotsdir("native/sensitivity_exps/Comp_adve_vel_flipped.png"), bbox_inches = "tight")


fig, ax = plt.subplots(1, 2, figsize = (20, 7.5), sharey = false)
vars =  ["only_wind",  "only_init", "only_kappa"]
mult = 100 * 86400 * 365 * 10 
for (i, expname) in enumerate(vars)

    ax[1].plot(tecco, mult .* (wθ_eul[expname] .+ wθ_bol[expname]), c = exp_colors[expname], 
    label = plot_labels_effects[expname], linewidth = 3)

    ax[2].plot(tecco, 1e-6 .* ΨEul[expname], c = exp_colors[expname], 
    label = plot_labels_effects[expname], linewidth = 3)
end
ax[1].set_ylabel("[cK per decade]")
ax[2].set_ylabel("[Sv]")

ax[1].set_title("Vertical Heat Advection Convergence")
ax[2].set_title("Vertical Transport [z = 2500]")

[a.set_ylim(-2.5, 2.5) for a in ax]
[a.legend(frameon = false) for a in ax]

fig.savefig(plotsdir("native/sensitivity_exps/Comp_adve_vel.png"), bbox_inches = "tight")
fig