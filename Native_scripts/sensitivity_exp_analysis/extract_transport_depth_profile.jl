include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, DSP, PyCall
import PyPlot as plt 
include(srcdir("config_exp.jl"))

include(srcdir("MeshArraysPlots.jl"))

tecco = 1992+1/24:1/12:2018
nz = 50

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
λ_wrap = wrap_λ(λ)
ocean_mask = wet_pts(Γ)

region = "PAC"
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)

include(srcdir("plot_and_dir_config.jl"))

NW_PAC = region_mask(ocean_mask, λ_wrap, ϕ, (0, 90), (0, 190))
NE_PAC = region_mask(ocean_mask, λ_wrap, ϕ, (0, 90), (190, 360))

cs, sn = get_cs_and_sn(γ)

ϕ_avg = zonal_average(ϕ, area .* PAC_msk)
ϕ_avg = ϕ_avg[isfinite.(ϕ_avg)]
NPAC_boundidx = Base.findmin(abs.(ϕ_avg .- 23))[2]

function obtain_W(Ψ, NPAC_boundidx)
    Ws = zeros(50, 312)
    for k = 1:50
        Ws[k, :] = Ψ[k, NPAC_boundidx, :]
    end
    return -Ws
end

open_streamfunction(expname, type) = jldopen(datadir("Ψ_" *type * "_timeseries_PAC_" * expname *".jld2"))["Ψ_exp_timeseries"]
obtain_W(expname, type, NPAC_boundidx) = obtain_W(open_streamfunction(expname, type), NPAC_boundidx)


lin_exps = ["iter0_bulkformula", "only_init", "only_kappa", "only_sfc", "iter129_bulkformula", "SUM"]; nexps = length(lin_exps)
plot_labels_effect = ["CTRL (Iteration 0)", "INIT Effect", "MIXING Effect", "FORCING Effect", "FULL (Iteration 129)", "SUM"]
plot_labels_effect = Dict(lin_exps[i] => plot_labels_effect[i] for i in 1:nexps)
lin_exps = ["only_init", "iter0_bulkformula", "only_kappa", "iter129_bulkformula", "only_sfc", "SUM"]; nexps = length(lin_exps)

type = "EulBol"
Ws = Dict()
[Ws[ex] = 1e-6 .* obtain_W(ex, type, NPAC_boundidx) for ex in lin_exps[1:end-1]]
[Ws[ex] .= Ws[ex] .- Ws["iter0_bulkformula"] for ex in ["only_sfc", "only_init", "only_kappa"]]
Ws["SUM"] = Ws["iter0_bulkformula"] .+ Ws["only_sfc"] .+ Ws["only_init"] .+ Ws["only_kappa"]

fig, axs = plt.subplots(2, 3, figsize=(20,13), sharey = true)
CF = Any[]
for (i, exp) in enumerate(lin_exps)
    ax = axs[i]
    cb = ax.pcolormesh(tecco, z, Ws[exp], cmap = cmo.balance, vmin = -10, vmax = 10, shading = "gourad")
    ax.set_ylim(1000, 4000); ax.set_yticks([1000, 2000, 3000, 4000])
    ax.axhline(2000, c = "k", alpha = 0.35, linestyle = "--")
    ax.axhline(3000, c = "k", alpha = 0.35, linestyle = "--")
    ax.invert_yaxis()
    ax.set_title(plot_labels_effect[exp])
    ax.set_xlabel("time", fontweight = "bold") 
    push!(CF, cb)
end
[ax.set_ylabel("Depth [m]", fontweight = "bold") for ax in axs[1:2]]
fig.subplots_adjust(hspace=0.3)
fig.colorbar(CF[1], ax = axs[:], orientation = "horizontal", label = "[Sv]", 
fraction = 0.03, extend = "both")
fig.savefig(plotsdir("native/sensitivity_exps/EulBol_Vertical_Transport_FullDepth.png"), dpi= 400, bbox_inches = "tight")

type = "Eul"
Ws = Dict()
[Ws[ex] = 1e-6 .* obtain_W(ex, type, NPAC_boundidx) for ex in lin_exps[1:end-1]]
[Ws[ex] .= Ws[ex] .- Ws["iter0_bulkformula"] for ex in ["only_sfc", "only_init", "only_kappa"]]
Ws["SUM"] = Ws["iter0_bulkformula"] .+ Ws["only_sfc"] .+ Ws["only_init"] .+ Ws["only_kappa"]

fig, axs = plt.subplots(2, 3, figsize=(20,13), sharey = true)
CF = Any[]
for (i, exp) in enumerate(lin_exps)
    ax = axs[i]
    cb = ax.pcolormesh(tecco, z, Ws[exp], cmap = cmo.balance, vmin = -10, vmax = 10, shading = "gourad")
    ax.set_ylim(1000, 4000); ax.set_yticks([1000, 2000, 3000, 4000])
    ax.axhline(2000, c = "k", alpha = 0.35, linestyle = "--")
    ax.axhline(3000, c = "k", alpha = 0.35, linestyle = "--")
    ax.invert_yaxis()
    ax.set_title(plot_labels_effect[exp])
    ax.set_xlabel("time", fontweight = "bold") 
    push!(CF, cb)
end
[ax.set_ylabel("Depth [m]", fontweight = "bold") for ax in axs[1:2]]
fig.subplots_adjust(hspace=0.3)
fig.colorbar(CF[1], ax = axs[:], orientation = "horizontal", label = "[Sv]", 
fraction = 0.03, extend = "both")
fig.savefig(plotsdir("native/sensitivity_exps/Eul_Vertical_Transport_FullDepth.png"), dpi= 400, bbox_inches = "tight")
fig