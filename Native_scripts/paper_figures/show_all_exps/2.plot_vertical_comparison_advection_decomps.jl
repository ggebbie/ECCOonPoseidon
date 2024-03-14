include("../../../src/intro.jl")

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

vars =   ["iter0_bulkformula", "only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]

wθ_resid = Dict(); wθ_eul = Dict(); wθ_bol = Dict()

Δwθ = Dict()
wΔθ = Dict()
ΔwΔθ = Dict()

integrate(start, x) = cumsum([start, x...])[1:end-1]
include(srcdir("plot_and_dir_config.jl"))

sns.set_theme(context = "paper", style = "ticks",
              palette = colors, rc = custom_params);

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
E,F = trend_matrices(tecco)
get_trend(x) = (F * x)[2]
fig, axs = plt.subplots(figsize = ( 6, 4))
axs.set_xlabel("time", fontweight = "bold")
axs.set_title("Decomposition of " * L"\mathbf{A}_Z (w_{129}^{res}, \theta_{129})")
expname =   "iter129_bulkformula"
lw = 2.0

axs.plot(tecco, wθ_resid["iter129_bulkformula"], label = L"\mathbf{A}_Z (w_{129}^{res}, \theta_{129})", 
linewidth = lw, color = exp_colors[expname])
axs.plot(tecco, wθ_resid["iter0_bulkformula"], label = L"\mathbf{A}_Z (w_0^{res}, \theta_0)", 
linewidth = lw, color = exp_colors["iter0_bulkformula"])

tm1 = Δwθ[expname] .-ΔwΔθ[expname]
println((F * tm1)[2])
axs.plot(tecco, Δwθ[expname] .-ΔwΔθ[expname] , label = L"\mathbf{A}_Z (\Delta w^{res}, \theta_0)", 
linewidth = 3, color = 0.5 .* exp_colors[expname], alpha = 0.8)
tm2 = wΔθ[expname] .-ΔwΔθ[expname]
println((F * tm2)[2])

axs.plot(tecco, wΔθ[expname] .-ΔwΔθ[expname] , label = L"\mathbf{A}_Z (w_0^{res}, \Delta \theta)",
linewidth = lw, color = exp_colors[expname], alpha = 0.6)
tmp = (Δwθ[expname] .-ΔwΔθ[expname]) 
tmp = (wθ_resid["iter129_bulkformula"] .- wθ_resid["iter0_bulkformula"]) .- tmp
println((F * tmp)[2])

axs.plot(tecco, tmp , label = L"\mathbf{A}_Z (  \Delta w^{res}, \Delta \theta)", 
linewidth = lw, color = exp_colors[expname], alpha = 0.3)
axs.legend(frameon = false, ncols = 1, loc = "center right",  bbox_to_anchor=(1.35, 0.5))
axs.axhline(0, c = "k", linestyle = "--", alpha = 0.2)

fig
axs.set_ylabel("[cK]", fontweight = "bold")
# fig.subplots_adjust(hspace = 0.5)
axs.grid()
fig
fig.savefig(plotsdir("native/paper_figures/2.HeatBudgetVertΔDecomp_all.png"), bbox_inches = "tight", dpi = 400)

fig


vars =   [ "only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]

Δwθ0 = Dict()
[Δwθ0[nexpt] = Δwθ[nexpt] .- ΔwΔθ[nexpt] for nexpt in vars]

Δwθ0["SUM"] = sum([Δwθ0[nexpt] for nexpt in vars[1:end-1]])
Δwθ0["eps"] = Δwθ0["iter129_bulkformula"] .- Δwθ0["SUM"]
# Δwθ0["eps"]
vars = vcat(vars,"SUM")
plot_labels["eps"] = "ϵ"
exp_colors["eps"] = "grey"

fig, axs = plt.subplots(figsize = ( 7, 4))
axs.set_title("Effect of Control Adjustments on \n " * L" \mathbf{A_Z (\Delta w^{res}, \theta_0)}")
axs.set_xlabel("time", fontweight = "bold")
vars =   [ "only_init", "only_kappa", "only_wind", "only_buoyancy", "eps"]

axs.set_ylabel(" [cK]", fontweight = "bold") 

plot_labels_list = ["Initial Condition\nAdjustment Effect", "Mixing Parameter\nAdjustment Effect", 
"Wind Stress\nAdjustment Effect ", "Buoyancy Forcing\nAdjustment Effect ", "Non-Linear\nAdjustment Effects"]
for (i, expname) in enumerate(vars)
    println(expname)
    lw = 2.5
    tmp = Δwθ0[expname]
    println(round(get_trend(tmp), digits = 3))
    axs.plot(tecco, tmp, label = plot_labels_list[i], 
    linewidth = lw, color = exp_colors[expname])
    axs.axhline(0, c = "k", linestyle = "--")
end
axs.grid()
fig
fig.subplots_adjust(hspace = 0.5)
axs.legend(frameon = false, ncols = 3, loc = "lower center",  bbox_to_anchor=(0.5, -0.42))
fig
fig.savefig(plotsdir("native/paper_figures/2.HeatBudgetVertΔWDecomp_Adjustments_all.png"), bbox_inches = "tight", dpi = 400)

fig


