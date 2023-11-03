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

@pyimport cmocean.cm as cmos


region = "PAC"; 
tecco = collect(1992+1/24:1/12:2018)

tecco = 1992+1/24:1/12:2018
nz = 50; tstart = 1; tstop = 312; nt = tstop - tstart + 1;


uplvl = -1e3; botlvl = -Inf
lvls = findall( botlvl .<= z[:].<= uplvl)
Y_mask = reverse(z[lvls])

cp = Vector{Any}(missing, 1)  
# alabels = ["Iteration 129", "Iteration 0", "Initial 129", "Initial 0"]
clevels = ["-9.0", "-7.5", "-1.5"]


vars =  ["only_init", "only_kappa", "only_wind", "iter129_bulkformula",  "iter0_bulkformula"]
Ψ_means = Dict()
for (i, expname) in enumerate(vars)
    read_file = jldopen(datadir("Ψ_Eul_timeseries_"*region*"_" * expname *".jld2"))
    Ψ_expEul = read_file["Ψ_exp_timeseries"]

    read_file = jldopen(datadir("Ψ_Eul_timeseries_"*region*"_" * expname *".jld2"))
    Ψ_expEulBol = read_file["Ψ_exp_timeseries"]
    Ψ_exp = Ψ_expEulBol .- Ψ_expEul
    Ψ_mean = mean(Ψ_exp, dims = 3)[:, :, 1]
    Ψ_means[expname] = Ψ_mean
end

Ψ_means["only_init"] .= Ψ_means["only_init"] .- Ψ_means["iter0_bulkformula"]
Ψ_means["only_kappa"] .= Ψ_means["only_kappa"] .-Ψ_means["iter0_bulkformula"]
Ψ_means["only_wind"] .= Ψ_means["only_wind"] .-Ψ_means["iter0_bulkformula"]
# Ψ_means["SUM"] = Ψ_means["only_init"] .+ Ψ_means["only_kappa"] .+ Ψ_means["only_sfc"] .- 2*Ψ_means["iter0_bulkformula"]

ϕ_avg = jldopen(datadir("Ψ_Eul_timeseries_"*region*"_" * "only_init" *".jld2"))["ϕ_avg"]

fig,axs=plt.subplots(
    1, 3, figsize = (21, 7.5))
vars =  ["only_init", "only_kappa", "only_wind"]

for (i, expname) in enumerate(vars)
    Ψ_exp = 1e-6.* Ψ_means[expname]

    ax = axs[i]

    ax.set_facecolor("black")
    Ψ_bounds = 10.5
    levels = -3:0.5:3
    ax.contourf(ϕ_avg, z,  Ψ_exp, cmap=cmos.delta,levels = levels, 
    vmin = -Ψ_bounds, vmax = Ψ_bounds, extend = "both")
    cs2 = ax.contour(ϕ_avg, z[30:46], Ψ_exp[30:46, :], colors="k",levels = levels)
    ax.contour(ϕ_avg, z[1:30], Ψ_exp[1:30, :], colors="k",levels = levels)
    ax.contour(ϕ_avg, z[46:end], Ψ_exp[46:end, :], colors="k",levels = levels)
    labels = ax.clabel(cs2, fontsize=20, inline=true, fmt = "%.1f", 
    inline_spacing = 12, rightside_up = true, use_clabeltext = true)
    ax.invert_yaxis()
    ax.set_xticks(-40:20:60)
    ax.set_xlim(-34, 60)
    ax.set_title(plot_labels[expname] * " Effect", fontweight = "bold")
    ax.set_ylabel("Depth [m]")
end
fig
fig.savefig(plotsdir("native/sensitivity_exps/5.ΨBol_Effect.png"))
