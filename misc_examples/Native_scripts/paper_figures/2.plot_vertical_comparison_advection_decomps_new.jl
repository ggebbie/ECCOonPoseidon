include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))
include(srcdir("plot_and_dir_config.jl"))

# cmo = pyimport("cmocean.cm");
@pyimport seaborn as sns;

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)
region = "NPAC"; 

uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= -z[:].<= uplvl)
tecco = 1992+1/24:1/12:2018

vars = ["iter0_bulkformula", "only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]

wθ_resid = Dict(); wθ_eul = Dict(); wθ_bol = Dict()

Δwθ = Dict()
wΔθ = Dict()
ΔwΔθ = Dict()

integrate(start, x) = cumsum([start, x...])[1:end-1]
include(srcdir("plot_and_dir_config.jl"))
sns.set_theme(context = "paper", style = "ticks",
              palette = colors, rc = custom_params);

wpθ = Dict()
wpΔθ = Dict()

function integration(t, x)
    int_x = 3.154e+7* cumul_integrate(tecco, x)
    int_x .-= int_x[1]
    return 100 * int_x
end

for (i, expname) in enumerate(vars)
    fname = datadir("native/" * expname * region * "_THETA_budget_Wθ_eul_bol_ΔDecomp" * suffix * ".jld2")
    print(expname)
    vars = jldopen(fname)["dθ"]

    wθfull = vars["wθBot"] .- vars["wθTop"]
    wθ_resid[expname] = 100 .* integrate(0, wθfull).* 2.628f+6
    Δwθ[expname] = 100 .* integrate(0, vars["ΔW_θBot"] .- vars["ΔW_θTop"]).* 2.628f+6
    wΔθ[expname] = 100 .* integrate(0, vars["W_ΔθBot"] .- vars["W_ΔθTop"]).* 2.628f+6
    ΔwΔθ[expname] = 100 .* integrate(0, vars["ΔW_ΔθBot"] .- vars["ΔW_ΔθTop"]).* 2.628f+6
end

vars =   [ "only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]

fig, axs = plt.subplots(1, 3, figsize = ( 8, 5), sharex = false, sharey = true)
# axs.set_xlabel("time", fontweight = "bold")
axs[1].set_title("Vertical Heat Flux By\nExternal Advection\nProcess\n" * L"\mathbf{A}_Z (w_0^{res} + \Delta w^{res}, \theta_0 + \Delta  \theta)")
axs[2].set_title("First-Guess Temperature \n Contribution \n" * L"\mathbf{A}_Z (w_0^{res} + \Delta w^{res}, \theta_0)")
axs[3].set_title("Temperature Adjustment \n Contribution \n" * L"\mathbf{A}_Z (w_0^{res} + \Delta w^{res}, \Delta  \theta)")

plot_labels = ["Initial Condition\nAdjustments", "Mixing\nAdjustment", "Wind Forcing\nAdjustment", 
"Buoyancy Forcing\nAdjustment", "Iteration 129\n(All Adjustments)"]

E, F = trend_matrices(tecco)
trend(x) = (F * x[:])[2]
trends_dict = Dict()

for (i, expname) in enumerate(vars)
    trends_dict[expname] = zeros(3)

    lw = 2.0

    axs[1].plot(tecco, wθ_resid[expname], label = L"\mathbf{A}_Z (w^{res}, \theta^\#)" * "(Iteration 129)", 
    linewidth = lw, color = exp_colors[expname])
    trends_dict[expname][1] =  trend(wθ_resid[expname])

    # axs.plot(tecco, wθ_resid["iter0_bulkformula"], label = L"\mathbf{A}_Z (w^{res}, \theta^\#)" * "(Iteration 0)", 
    # linewidth = lw, color = exp_colors["iter0_bulkformula"])

    axs[2].plot(tecco, Δwθ[expname] .+  wθ_resid["iter0_bulkformula"] .-  ΔwΔθ[expname], label = plot_labels[i], 
    linewidth = lw, color = exp_colors[expname], alpha = 1)
    trends_dict[expname][2] =  trend(Δwθ[expname] .+  wθ_resid["iter0_bulkformula"] .-  ΔwΔθ[expname])

    tmp =   wθ_resid[expname] .- (Δwθ[expname] .+  wθ_resid["iter0_bulkformula"] .-  ΔwΔθ[expname])

    axs[3].plot(tecco, tmp , label = " Heat flux due to  " *  L"\theta", 
    linewidth = lw, color = exp_colors[expname], alpha =1)
    trends_dict[expname][3] =  trend(tmp)


end
trends_df = DataFrame(trends_dict)
trends_df = 1 .* round.(trends_df, digits = 3)
# trends_df.Difference .= trends_df.iter129_bulkformula - trends_df.iter0_bulkformula 
trends_df[end, :] .= values(trends_df[1, :]) .- values(trends_df[2, :]) 
trends_df
[a.grid() for a in axs]
fig
fig_labs = uppercase.(["a", "b", "c", "d", "e"])
for (i, a) in enumerate(axs)
    a.annotate(fig_labs[i], (0.05, 0.05), fontsize = 15, 
    xycoords="axes fraction", fontweight = "bold")
end
axs[1].set_ylabel("[cK]", fontweight = "bold")
fig.subplots_adjust(hspace = 0.5)
fig
axs[2].legend(loc = "lower center", frameon = false, bbox_to_anchor = (0.5, -0.4), ncols = 3, alignment = "center")
fig.tight_layout()
fig.subplots_adjust(wspace = 0.1)

fig
fig.savefig(plotsdir("native/paper_figures/2.HeatBudgetVertΔWDecomp.png"), bbox_inches = "tight", dpi = 400)

fig


