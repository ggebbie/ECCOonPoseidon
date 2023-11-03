include("../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour, PyCall, 
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))
@pyimport seaborn as sns;

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

include(srcdir("plot_and_dir_config.jl"))

@pyimport cmocean.cm as cmo


tecco = collect(1992+1/24:1/12:2018)
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= -z[:].<= uplvl)

region = "NPAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)
ctrl_vol = sum(cell_volumes[:, lvls])

region = "PAC"; 

cp = Vector{Any}(missing, 1)  
# alabels = ["Iteration 129", "Iteration 0", "Initial 129", "Initial 0"]
clevels = ["-9.0", "-7.5", "-1.5"]

ϕ_avg = jldopen(datadir("Ψ_Eul_timeseries_"*region*"_" * "only_init" *".jld2"))["ϕ_avg"]
ϕ_ref = findall( 23 .<= ϕ_avg .<= 40)
z_ref = findall( -2700 .<= -z[:].<= -2300)[1]

vars =  ["iter0_bulkformula", "iter129_bulkformula","only_init", "only_kappa", "only_wind", "only_buoyancy"]
ΨEul = Dict(); ΨBol = Dict()


for (i, expname) in enumerate(vars)
    read_file = jldopen(datadir("Ψ_EulBol_timeseries_"*region*"_" * expname *".jld2"))
    Ψ_EulBol = read_file["Ψ_exp_timeseries"]

    read_file = jldopen(datadir("Ψ_Eul_timeseries_"*region*"_" * expname *".jld2"))
    Ψ_Eul = read_file["Ψ_exp_timeseries"]

    # Ψs[expname] = Ψ_exp[z_ref, ϕ_ref, :][:]
    Ψ_Bol = Ψ_EulBol .- Ψ_Eul
    ΨEul[expname] = -Ψ_Eul[:, ϕ_ref[1], :]
    ΨBol[expname] = -Ψ_Bol[:, ϕ_ref[1], :]

end


effect_exp =  ["only_init", "only_kappa", "only_wind"]

[ΨEul[expt] .-= ΨEul["iter0_bulkformula"] for expt in effect_exp]
[ΨBol[expt] .-= ΨBol["iter0_bulkformula"] for expt in effect_exp]

fig, ax = plt.subplots(1, 2, figsize = (15, 5), sharex = true, sharey = true)
vars =  ["only_wind"]

for (i, expname) in enumerate(vars)
    Ψ_exp = 1e-6.* ΨEul[expname]
    vmax = maximum(abs.(Ψ_exp))
    ax[1].pcolormesh(tecco, z, Ψ_exp, cmap = cmo.balance, label = plot_labels_effects[expname], linewidth = 2.5, vmin = -vmax, vmax = vmax)
    Ψ_exp = 1e-6.* ΨBol[expname]
    ax[2].pcolormesh(tecco, z, Ψ_exp, cmap = cmo.balance, label = plot_labels_effects[expname], linewidth = 2.5)

    # ax.set_title(plot_labels[expname] * " Effect", fontweight = "bold")
end
ax[1].invert_yaxis()
ax[1].set_ylabel("[Sv]")
ax[1].legend(frameon = false)
ax[1].set_title("Eulerian Streamfunction Anomaly")
ax[2].set_title("Bolus Streamfunction Anomaly")

# fig.savefig(plotsdir("native/sensitivity_exps/ΨEul_Effect.png"))
fig


fig, ax = plt.subplots(1, 2, figsize = (15, 5), sharex = true, sharey = true)
vars =  ["only_wind"]

for (i, expname) in enumerate(vars)
    Ψ_exp = 1e-6.* ΨEul[expname]
    vmax = maximum(abs.(Ψ_exp))
    ax[1].plot(Ψ_exp[:, 37], z, label = plot_labels_effects[expname], linewidth = 2.5)
    Ψ_exp = 1e-6.* ΨBol[expname]
    ax[2].plot(Ψ_exp[:, 37], z, label = plot_labels_effects[expname], linewidth = 2.5)

    # ax.set_title(plot_labels[expname] * " Effect", fontweight = "bold")
end
ax[1].invert_yaxis()
ax[1].set_ylabel("[Sv]")
ax[1].legend(frameon = false)
ax[1].set_title("Eulerian Streamfunction Anomaly")
ax[2].set_title("Bolus Streamfunction Anomaly")

# fig.savefig(plotsdir("native/sensitivity_exps/ΨEul_Effect.png"))
fig
# tecco[140]