
include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    LaTeXStrings, PyCall, RollingFunctions
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))

@pyimport seaborn as sns;

include(srcdir("plot_and_dir_config.jl"))

colors =  sns.color_palette("deep")[1:4]
sns.set_theme(context = "poster", font_scale = 1.0,
              palette = sns.color_palette("deep"));
sns.set_style("white")

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
tecco = 1992+1/24:1/12:2018

ocean_mask = wet_pts(Γ)
region = "PAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "false", include_bering = true)

ϕ_avg = zonal_average(ϕ, area .* PAC_msk)
ϕ_avg = ϕ_avg[isfinite.(ϕ_avg)]

region = "NPAC"; 
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall(botlvl .< -z[:] .< uplvl)

NPAC_boundidx = Base.findmin(abs.(ϕ_avg .- 23))[2]

obtain_V(Ψ, NPAC_boundidx) = vcat(Ψ[2:end, NPAC_boundidx, :] .- Ψ[1:end-1, NPAC_boundidx, :], -Ψ[end, NPAC_boundidx, :]')

function obtain_W(Ψ, NPAC_boundidx)
    Ws = Ψ[:, NPAC_boundidx, :]
    return -Ws
end

open_streamfunction(expname, type) = Float32.(jldopen(datadir("Ψ_" *type * "_timeseries_PAC_" * expname *".jld2"))["Ψ_exp_timeseries"])
get_transports(Ψ, NPAC_boundidx) = (obtain_V(Ψ, NPAC_boundidx), obtain_W(Ψ, NPAC_boundidx))

mean_pre2005(x) = Float32.(mean(x[:, 1:156], dims = 2)[:, 1])
mean_post2005(x) = Float32.(mean(x[:, 156:end], dims = 2)[:, 1])

round_mean(x) = round(mean(x), digits = 2)
include(srcdir("plot_and_dir_config.jl"))

vars =  ["iter0_bulkformula", "only_init", "only_kappa", "only_sfc", "iter129_bulkformula"]
var_labels = [plot_labels[expname] for expname in vars]
var_labels[1] = "CTRL"; var_labels[end] = "FULL"

fig, ax = plt.subplots(2, 2, figsize = (15, 10), sharey = true)
for (i, expname) in enumerate(vars)
    println(expname)

    ΨEulBol = open_streamfunction(expname, "EulBol")
    V, W = 1e-6 .* get_transports(ΨEulBol, NPAC_boundidx)

    ax[1, 1].plot(mean_pre2005(V), z, color = exp_colors[expname], alpha = 0.3); ax[1, 1].set_title("Merdional Transport \n 1993 - 2005")
    ax[1, 2].plot(mean_pre2005(W), z, color = exp_colors[expname], alpha = 0.3); ax[1, 2].set_title("Vertical Transport \n 1993 - 2005")

    ax[1, 1].plot(mean_post2005(V), z, color = exp_colors[expname]); ax[2, 1].set_title("Merdional Transport \n 2005 - 2017")
    ax[1, 2].plot(mean_post2005(W), z, color = exp_colors[expname], label = expname); ax[2, 2].set_title("Vertical Transport \n 2005 - 2017")

end
[a.set_xlim(-1, 0) for a in ax[:, 1]]
[a.set_xlim(0, 3) for a in ax[:, 2]]

ax[1].set_ylim(1700, 3200)
ax[1].invert_yaxis()
ax[4].legend()

fig



vars =  ["only_init", "only_kappa", "only_sfc", "iter129_bulkformula", "iter0_bulkformula"]
var_labels = [plot_labels[expname] for expname in vars]
# var_labels[end-1] = "CTRL"; var_labels[end] = "FULL"

fig, ax = plt.subplots(1, 2, figsize = (15, 7.5), sharey = true)
for (i, expname) in enumerate(vars)
    println(expname)
    ΨEulBol = open_streamfunction(expname, "Eul")
    V, W = get_transports(ΨEulBol, NPAC_boundidx)

    ax[1].plot(1e-6 .* mean(V, dims = 2)[:], cumsum(Γ.DRF), color = exp_colors[expname], label = var_labels[i]); 
    ax[1].set_title("Merdional Transport [Sv]")
    ax[2].plot(1e-6 .* mean(W, dims = 2)[:], cumsum(Γ.DRF), color = exp_colors[expname], label = var_labels[i]); 
    ax[2].set_title("Vertical Transport [Sv]")
end

ax[1].set_xlim(-0.6, 0); ax[2].set_xlim(0, 3); ax[1].set_ylabel("Depth [m]")
ax[1].set_yticks([1500, 2000, 2500, 3000])
ax[1].set_ylim(1499,  3250)
[a.legend(frameon = false) for a in ax]
ax[1].invert_yaxis()

fig.savefig(plotsdir("native/sensitivity_exps/Transport_Profile_Eulerian.png"))

fig

vars =  ["only_init", "only_kappa", "only_sfc", "iter129_bulkformula", "iter0_bulkformula"]
var_labels = [plot_labels[expname] for expname in vars]
# var_labels[end-1] = "CTRL"; var_labels[end] = "FULL"

fig, ax = plt.subplots(1, 2, figsize = (15, 7.5), sharey = true)
for (i, expname) in enumerate(vars)
    println(expname)
    ΨEulBol = open_streamfunction(expname, "EulBol")

    ΨEul = open_streamfunction(expname, "Eul")
    ΨBol = ΨEulBol .- ΨEul
    V, W = get_transports(ΨBol, NPAC_boundidx)

    ax[1].plot(1e-6 .* mean(V, dims = 2)[:], cumsum(Γ.DRF), color = exp_colors[expname], label = var_labels[i]); 
    ax[1].set_title("Merdional Transport [Sv]")
    ax[2].plot(1e-6 .* mean(W, dims = 2)[:], cumsum(Γ.DRF), color = exp_colors[expname], label = var_labels[i]); 
    ax[2].set_title("Vertical Transport [Sv]")
end

ax[1].set_xlim(-1, 1); ax[2].set_xlim(-1, 1); ax[1].set_ylabel("Depth [m]")
ax[1].set_yticks([1500, 2000, 2500, 3000])
ax[1].set_ylim(1499,  3250)
[a.legend(frameon = false) for a in ax]
ax[1].invert_yaxis()

fig.savefig(plotsdir("native/sensitivity_exps/Transport_Profile_Bolus.png"))

fig