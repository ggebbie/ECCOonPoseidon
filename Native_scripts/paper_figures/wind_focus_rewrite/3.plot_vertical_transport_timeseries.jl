using Pkg
Pkg.activate(".")

include("../../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson,
    BenchmarkTools, LaTeXStrings, PyCall, DataFrames
import NaNMath as nm
import PyPlot as plt
import NumericalIntegration

include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");

#blue, red, #green, orangeDataFrames

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)

region = "PAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "false", include_bering = true)

ϕ_avg = zonal_average(ϕ, area .* PAC_msk)
ϕ_avg = ϕ_avg[isfinite.(ϕ_avg)]

NPAC_boundidx = Base.findmin(abs.(ϕ_avg .- 23))[2]
lvls = findall( -3000 .<= -z[:].<= -2000)

obtain_Vin(Ψ, NPAC_boundidx) = Ψ[lvls[end] + 1, NPAC_boundidx, :] .- Ψ[lvls[1], NPAC_boundidx, :]
function obtain_Wtop(Ψ, NPAC_boundidx)
    Win = Ψ[lvls[1], NPAC_boundidx, :]
    return -Win
end

function obtain_Wbot(Ψ, NPAC_boundidx)
    Wout = Ψ[lvls[end] + 1, NPAC_boundidx, :]
    return -Wout
end

open_streamfunction(expname, type) = jldopen(datadir("Ψ_" *type * "_timeseries_PAC_" * expname *".jld2"))["Ψ_exp_timeseries"]

get_face_transports(Ψ, NPAC_boundidx) = (obtain_Vin(Ψ, NPAC_boundidx), obtain_Wtop(Ψ, NPAC_boundidx), obtain_Wbot(Ψ, NPAC_boundidx))


region = "NPAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)

cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)

lvls = findall( -3000 .<= -z[:].<= -2000)
V = sum(cell_volumes[:, lvls])
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
# lvls = findall( botlvl .<= z[:].<= uplvl)
tecco = 1992+1/24:1/12:2018

include(srcdir("plot_and_dir_config.jl"))

exps =  [ "mean_tau_noadjust_redo_bf", "mean_tau_yesadjust_redo_bf"]
Vin_d  = Dict()
Win_d  = Dict()
Wout_d  = Dict()

for (i, expname) in enumerate(exps)
    Ψ = open_streamfunction(expname, "EulBol")
    Vin, Wout, Win = get_face_transports(Ψ, NPAC_boundidx) 
    Vin_d[expname] = 1e-6 .* obtain_Vin(Ψ, NPAC_boundidx) 
    Wout_d[expname] = 1e-6 .* obtain_Wtop(Ψ, NPAC_boundidx)
    Win_d[expname] = 1e-6 .* obtain_Wbot(Ψ, NPAC_boundidx)
end

sns.set_style("darkgrid", Dict("axes.facecolor" => ".95", "grid.color" => "0.5"))

fig, axs = plt.subplots(1, 2, figsize = (12, 6), sharey = true)

axs[1].set_title("Top Boundary Residual Transport\n" * L"W^{res}_{\mathcal{T}}")
axs[2].set_title("Bottom Boundary Residual Transport\n" * L"W^{res}_{\mathcal{B}}")

titles = [L"\Delta^{{\tau}} \overline{W^{res}}}", 
L"\Delta^{\tau} \overline{W}}", L"\Delta^{\tau} \overline{W^{*}}}"]


axs[1].tick_params(which="both", bottom=true, left = true)

# axs[1].set_ylabel("[cK]", fontweight = "bold")
[ax.set_xlabel("time", fontweight = "bold") for ax in axs]
fig.tight_layout()

lw = 1.5; α = 0.95

exps =  [ "mean_tau_noadjust_redo_bf", "mean_tau_yesadjust_redo_bf", "Difference"]

# exps =  ["iter0_bulkformula", "iter129_bulkformula"]

plot_labels_list = [
"Seasonal First-Guess Wind", 
"Seasonal Adjusted Wind", "Response to Seasonal\nWind Stress Adjustment"]
alphas = [0.3, 0.3, 0.9]

exp_colors["mean_tau_noadjust_redo_bf"] = exp_colors["iter0_bulkformula"] 
exp_colors["mean_tau_yesadjust_redo_bf"] = exp_colors["only_wind"] 
exp_colors["Difference"] = "k"

Wout_d["Difference"] = Wout_d["mean_tau_yesadjust_redo_bf"] .- Wout_d["mean_tau_noadjust_redo_bf"]
Win_d["Difference"] = Win_d["mean_tau_yesadjust_redo_bf"] .- Win_d["mean_tau_noadjust_redo_bf"]


E, F = trend_matrices(tecco)
trend(x) = (F * x[:])[2]
trends_dict = Dict()
for (i, expname) in enumerate(exps)
    trends_dict[expname] = zeros(3)
    println(expname)
    println(mean(Wout_d[expname][96:end]))
    println(mean(Win_d[expname][96:end]))

    axs[1].plot(tecco, Wout_d[expname], label = plot_labels_list[i], 
                linewidth = lw, alpha = alphas[i], c = exp_colors[expname])
    axs[2].plot(tecco, Win_d[expname], label = plot_labels_list[i], 
                linewidth = lw, alpha = alphas[i], c = exp_colors[expname])
end

fig_labs = uppercase.(["a", "b", "c", "d", "e"])
for (i, a) in enumerate(axs)
    a.annotate(fig_labs[i], (0.90, 0.05), fontsize = 25, 
    xycoords="axes fraction", fontweight = "bold")
end
axs[1].set_ylabel("Upward Transport [Sv]")
[a.grid(0.3) for a in axs]
fig.tight_layout()
fig.subplots_adjust(wspace = 0.07)
axs[1].legend(loc = "lower center", fontsize=17,
              frameon = false, bbox_to_anchor = (1, -0.6), 
              ncols = 2)
axs[1].grid(alpha = 0.3)
axs[2].grid(alpha = 0.3)

fig.savefig(plotsdir("wind_rewrite/6.vertical_transports_all.png"), 
            bbox_inches = "tight", dpi = 400)
fig
