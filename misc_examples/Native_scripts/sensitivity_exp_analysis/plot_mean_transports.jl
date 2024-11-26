
include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    LaTeXStrings, PyCall, RollingFunctions
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))

# cmo = pyimport("cmocean.cm");
@pyimport seaborn as sns;

include(srcdir("plot_and_dir_config.jl"))

colors =  sns.color_palette("deep")[1:4]
sns.set_theme(context = "talk", font_scale = 1.0,
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

obtain_Vin(Ψ, NPAC_boundidx) = Ψ[lvls[end] + 1, NPAC_boundidx, :] .- Ψ[lvls[1], NPAC_boundidx, :]
function obtain_Wtop(Ψ, NPAC_boundidx)
    Win = Ψ[lvls[1], NPAC_boundidx, :]
    return -Win
end

function obtain_Wbot(Ψ, NPAC_boundidx)
    Wout = Ψ[lvls[end] + 1, NPAC_boundidx, :]
    Wout = mean(Float32, Wout, dims = 1)[1, :]
    return -Wout
end

open_streamfunction(expname, type) = jldopen(datadir("Ψ_" *type * "_timeseries_PAC_" * expname *".jld2"))["Ψ_exp_timeseries"]
get_face_transports(Ψ, NPAC_boundidx) = (obtain_Vin(Ψ, NPAC_boundidx), obtain_Wtop(Ψ, NPAC_boundidx), obtain_Wbot(Ψ, NPAC_boundidx))

round_mean(x) = round(mean(x), digits = 2)
include(srcdir("plot_and_dir_config.jl"))

vars =  ["iter0_bulkformula", "only_init", "only_kappa", "only_sfc", "iter129_bulkformula"]
var_labels = [plot_labels[expname] for expname in vars]
var_labels[1] = "CTRL"; var_labels[end] = "FULL"

fig, ax = plt.subplots(1, 3, figsize = (15, 7.5), sharey = true)
ax[1].bar(1, NaN, color = "k", label = "Eulerian"); 
ax[1].bar(1, NaN, color = "k", alpha = 0.5, hatch="///", label = "Bolus"); 
# fig.suptitle("Bolus Transports")
for (i, expname) in enumerate(vars)
    println(expname)
    ΨEulBol = open_streamfunction(expname, "EulBol")
    ΨEul = open_streamfunction(expname, "Eul")
    ΨBol = ΨEulBol .- ΨEul

    Vin, Wout, Win = 1.f-6 .* mean.((get_face_transports(ΨEul, NPAC_boundidx) ))
    println(Vin)

    ax[1].bar(i, Vin, width=0.95, color = exp_colors[expname]); ax[1].set_title(L"V_{in}")
    ax[2].bar(i, Win, width=0.95, color = exp_colors[expname]); ax[2].set_title(L"W_{in}")
    ax[3].bar(i, Wout, width=0.95, color = exp_colors[expname]); ax[3].set_title(L"W_{out}")
    println(Vin + Win-Wout)

    Vin_bol, Wout_bol, Win_bol = 1.f-6 .* mean.((get_face_transports(ΨBol, NPAC_boundidx) ))
    println(Vin_bol- Wout_bol+ Win_bol)

    ax[1].bar(i, Vin_bol, width=0.95, bottom = Vin, color = exp_colors[expname], alpha = 0.5, hatch="////"); 
    ax[2].bar(i, Win_bol, width=0.95, bottom = Win, color = exp_colors[expname], alpha = 0.5, hatch="////"); 
    ax[3].bar(i, Wout_bol, width=0.95, bottom = Wout, color = exp_colors[expname], alpha = 0.5, hatch="////"); 

end
ax[1].set_ylim(-3, 3); ax[1].set_ylabel("Sv")
[a.set_xticks(1:length(vars), var_labels) for a in ax]
[a.set_xticklabels(var_labels, rotation = 45, va="center", position=(0,-0.05)) for a in ax]

ax[1].legend(ncols = 1, frameon = false)
fig
fig.savefig(plotsdir("native/sensitivity_exps/MeanEulBolTransports.png"))
