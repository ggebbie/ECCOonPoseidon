include("../../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson,
    BenchmarkTools, LaTeXStrings, PyCall, DataFrames
import NaNMath as nm
import PyPlot as plt
import NumericalIntegration

cumul_integrate = NumericalIntegration.cumul_integrate
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
sns.set_theme(context = "paper", style = "ticks",
              palette = colors, rc = custom_params);
exps =  ["iter0_bulkformula", "iter129_bulkformula", "only_buoyancy", 
              "only_init", "only_kappa", "only_wind"]
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

fig, axs = plt.subplots(2, 1, figsize = (8, 10), sharex = false, sharey = false)

# axs[1].set_ylabel("[cK]", fontweight = "bold")
[ax.set_xlabel("time", fontweight = "bold") for ax in axs]
fig.tight_layout()

lw = 1.0; α = 0.8

exps =  ["iter0_bulkformula", "only_buoyancy", 
"only_init", "only_kappa", "only_wind"]
# exps =  ["iter0_bulkformula", "iter129_bulkformula"]

plot_labels_list = ["Iteration 0", "Buoyancy Forcing", "Initial Condition", "Mixing Parameters", "Wind Stress"]
# plot_labels_list = ["Iteration 0", "Iteration 129"]

E, F = trend_matrices(tecco)
trend(x) = (F * x[:])[2]
trends_dict = Dict()
for (i, expname) in enumerate(exps)
    trends_dict[expname] = zeros(3)
    println(expname)
    axs[1].plot(tecco, Wout_d[expname] .-  Wout_d["iter0_bulkformula"], label = plot_labels_list[i], c = exp_colors[expname], linewidth = lw)
    axs[2].plot(tecco, Win_d[expname] .-  Win_d["iter0_bulkformula"], label = plot_labels_list[i], c = exp_colors[expname], linewidth = lw)
end
fig
Vin_d["only_init"]
fig_labs = uppercase.(["a", "b", "c", "d", "e"])
for (i, a) in enumerate(axs)
    a.annotate(fig_labs[i], (0.05, 0.05), fontsize = 15, 
    xycoords="axes fraction", fontweight = "bold")
end
trend(temps["iter0_bulkformula"])
trends_df = DataFrame(trends_dict)
trends_df = 1 .* round.(trends_df, digits = 3)
# trends_df.Difference .= trends_df.iter129_bulkformula - trends_df.iter0_bulkformula 
trends_df[!,:variable] = ["Top Advection", "Bottom Advection", "Horizontal Advection"]
trends_df
[a.grid() for a in axs]
axs[2].legend(loc = "lower center", frameon = false, bbox_to_anchor = (0.5, -0.35), ncols = 3)
fig.tight_layout()
fig.subplots_adjust(wspace = 0.1)
fig
fig.savefig(plotsdir("native/paper_figures/2.HeatBudgetTimeSeries_external_decomp_all.png"), bbox_inches = "tight", dpi = 1000)


fig, axs = plt.subplots(1, 3, figsize = (8, 5), sharex = false, sharey = true)
axs[1].set_title("Total Heat Flux By\nExternal Advection\nProcess\n" * L"\mathbf{A}^{top}")
axs[2].set_title("Vertical Advection \nContribution\n" * L"\mathbf{A}^{bot}")
axs[3].set_title("Zonal Advection \nContribution\n" * L"\mathbf{A}^{south}")

axs[1].set_ylabel("[cK]", fontweight = "bold")
[ax.set_xlabel("time", fontweight = "bold") for ax in axs]
fig.tight_layout()

integrate(start, x) = cumsum([start, x...])[1:end-1] .* (100 * 2.628e+6)

lw = 2.5; α = 0.8

exps =  ["iter0_bulkformula", "iter129_bulkformula", "only_buoyancy", 
"only_init", "only_kappa", "only_wind"]
plot_labels["Difference"] = "Difference"
exp_colors["Difference"] = "k"
plot_labels_list = ["Iteration 0", "Iteration 129", "Buoyancy Forcing", "Initial Condition", "Mixing Parameters", "Wind Stress"]
using LinearAlgebra, Statistics
E, F = trend_matrices(tecco)
trend(x) = (F * x[:])[2]
trends_dict = Dict()
for (i, expname) in enumerate(exps)
    println(expname)
    y = zeros(312)
    A = zeros(312, 2)
    y .= integrate(0, Vadvection[expname])
    # A[:, 1] .= integrate(0, Vin_d[expname])
    A[:, 1] .= integrate(0, Wout_d[expname])
    A[:, 2] .= integrate(0, Win_d[expname])
    coefs = pinv(A) * y
    println(coefs)
    println(cov(y, ( A * coefs)) / (std(y) * std( ( A * coefs))))

    axs[1].plot(tecco, coefs[1] * A[:, 1], label = plot_labels_list[i], c = exp_colors[expname], linewidth = lw)
    # trends_dict[expname][1] =  trend(integrate(0, Topadvection[expname]))
    axs[2].plot(tecco, coefs[2] * A[:, 2], label = plot_labels_list[i], c = exp_colors[expname], linewidth = lw)
    # trends_dict[expname][2] =  trend(integrate(0, Botadvection[expname]))

    # axs[3].plot(tecco, pinv(A[:, 1]) *integrate(0, Hadvection[expname]) * A[:, 1], label = plot_labels_list[i], c = exp_colors[expname])
    # trends_dict[expname][3] =  trend(integrate(0, Hadvection[expname]))

end

fig_labs = uppercase.(["a", "b", "c", "d", "e"])
for (i, a) in enumerate(axs)
    a.annotate(fig_labs[i], (0.05, 0.05), fontsize = 15, 
    xycoords="axes fraction", fontweight = "bold")
end
trend(temps["iter0_bulkformula"])
trends_df = DataFrame(trends_dict)
trends_df = 1 .* round.(trends_df, digits = 3)
# trends_df.Difference .= trends_df.iter129_bulkformula - trends_df.iter0_bulkformula 
trends_df[!,:variable] = ["Top Advection", "Bottom Advection", "Horizontal Advection"]
trends_df
[a.grid() for a in axs]
axs[2].legend(loc = "lower center", frameon = false, bbox_to_anchor = (0.5, -0.35), ncols = 3)
fig.tight_layout()
fig.subplots_adjust(wspace = 0.1)
fig
fig.savefig(plotsdir("native/paper_figures/2.HeatBudgetTimeSeries_external_decomp_all_reconstruct.png"), bbox_inches = "tight", dpi = 1000)



››