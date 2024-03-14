include("../../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall, Statistics
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

vars_names =   ["iter0_bulkformula", "only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]

wθ_resid_top = Dict(); wθ_eul = Dict(); wθ_bol = Dict()
wθ_resid_bot = Dict();
Δwθ_top = Dict()
wΔθ_top = Dict()
ΔwΔθ_top = Dict()

Δwθ_bot = Dict()
wΔθ_bot = Dict()
ΔwΔθ_bot = Dict()

integrate(start, x) = cumsum([start, x...])[1:end-1]
include(srcdir("plot_and_dir_config.jl"))

sns.set_theme(context = "paper", style = "ticks",
              palette = colors, rc = custom_params);

for (i, expname) in enumerate(vars_names)
    fname = datadir("native/" * expname * region * "_THETA_budget_Wθ_eul_bol_ΔDecomp" * suffix * ".jld2")
    print(expname)
    vars = jldopen(fname)["dθ"]

    wθ_resid_top[expname] = 100 .* integrate(0, - vars["wθTop"]).* 2.628f+6
    Δwθ_top[expname] = 100 .* integrate(0, - vars["ΔW_θTop"]).* 2.628f+6
    wΔθ_top[expname] = 100 .* integrate(0, - vars["W_ΔθTop"]).* 2.628f+6
    ΔwΔθ_top[expname] = 100 .* integrate(0, - vars["ΔW_ΔθTop"]).* 2.628f+6

    wθ_resid_bot[expname] = 100 .* integrate(0, vars["wθBot"]).* 2.628f+6
    Δwθ_bot[expname] = 100 .* integrate(0, vars["ΔW_θBot"] ).* 2.628f+6
    wΔθ_bot[expname] = 100 .* integrate(0, vars["W_ΔθBot"] ).* 2.628f+6
    ΔwΔθ_bot[expname] = 100 .* integrate(0, vars["ΔW_ΔθBot"] ).* 2.628f+6
end

vars =   [ "only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]
E,F = trend_matrices(tecco)
get_trend(x) = (F * x)[2]
var_explained(y, x) = 100 * ( 1 - (var(y - x) / var(y))) 
fig, ax = plt.subplots(1, 2, figsize = ( 8, 4), sharey = true)
expname =   "iter129_bulkformula"
lw = 1.5

wθ_resid_list = [wθ_resid_bot, wθ_resid_top]
Δwθ_list = [Δwθ_bot, Δwθ_top]
wΔθ_list = [wΔθ_bot, wΔθ_top]
ΔwΔθ_list = [ΔwΔθ_bot, ΔwΔθ_top]
ax[1].set_title("Decomposition of Bottom Boundary\n Advective Heat Flux Response to Optimization\n" * L"\Delta\mathbf{A}_{\mathcal{B}}")
ax[2].set_title("Decomposition of Upper Boundary\nAdvective Heat Flux Response to Optimization\n" * L"\Delta\mathbf{A}_{\mathcal{U}}")
fig
var_name = [L"\mathbf{A}^\mathcal{B}", L"\mathbf{A}^\mathcal{U}"]
for (i, axs) in enumerate(ax)
    wθ_resid = wθ_resid_list[i]
    Δwθ = Δwθ_list[i]
    wΔθ = wΔθ_list[i]
    ΔwΔθ = ΔwΔθ_list[i]
    axs.set_xlabel("time", fontweight = "bold")
    ΔA = wθ_resid["iter129_bulkformula"] .- wθ_resid["iter0_bulkformula"]
    axs.plot(tecco, ΔA, label = L"\Delta" * var_name[i], 
    linewidth = 1.5, color = exp_colors["iter129_bulkformula"])

    ΔA_Δw = Δwθ[expname] .-ΔwΔθ[expname]
    # println((F * tmp1)[2])
    axs.plot(tecco, ΔA_Δw, label =  var_name[i] * L"(\Delta w^{res}, \theta_0)", 
    linewidth = 1.5, color = "k", alpha = 0.8)
    ΔA_Δθ = wΔθ[expname] .-ΔwΔθ[expname]

    axs.plot(tecco, ΔA_Δθ , label = var_name[i] * L"(w_0^{res}, \Delta \theta)", 
    linewidth = lw, color = "blueviolet", alpha = 0.9)

    tmp = (ΔA_Δw) .+ (ΔA_Δθ)
    ΔA_Δw_Δθ = ΔA .- tmp

    axs.plot(tecco, ΔA_Δw_Δθ , label = var_name[i] * L"(  \Delta w^{res}, \Delta \theta)", 
    linewidth = lw, color = "hotpink", alpha = 0.9)
    axs.axhline(0, c = "k", linestyle = "--", alpha = 0.2)
    axs.grid()
    leg = axs.legend(frameon = false, ncols = 1, loc = "lower left", fontsize=12)

    println(var_explained(ΔA, ΔA_Δw))
    println(var_explained(ΔA, ΔA_Δθ))
    println(var_explained(ΔA, ΔA_Δw_Δθ))
end
ax[1].set_ylabel("[cK]", fontweight = "bold")
fig_labs = uppercase.(["a", "b", "c", "d", "e", "f"])
for (i, a) in enumerate(ax)
    a.annotate(fig_labs[i], (0.02, 0.92), fontsize = 20, color = "black", 
    xycoords="axes fraction", fontweight = "bold")
end
fig
# ax[1].legend(frameon = false, ncols = 2, loc = "center right",  bbox_to_anchor=(1.5, -0.2))
fig.subplots_adjust(wspace = 0.1)
fig
fig.savefig(plotsdir("native/paper_figures/2.HeatBudgetVertΔDecomp_faces_all.png"), bbox_inches = "tight", dpi = 400)

fig


vars =   [ "only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]

Δwθ0_Top = Dict()
Δwθ0_Bot = Dict()

[Δwθ0_Bot[nexpt] = Δwθ_bot[nexpt] .- ΔwΔθ_bot[nexpt] for nexpt in vars]
[Δwθ0_Top[nexpt] = Δwθ_top[nexpt] .- ΔwΔθ_top[nexpt] for nexpt in vars]

Δwθ0_Bot["SUM"] = sum([Δwθ0_Bot[nexpt] for nexpt in vars[1:end-1]])
Δwθ0_Bot["eps"] = Δwθ0_Bot["iter129_bulkformula"] .- Δwθ0_Bot["SUM"]

Δwθ0_Top["SUM"] = sum([Δwθ0_Top[nexpt] for nexpt in vars[1:end-1]])
Δwθ0_Top["eps"] = Δwθ0_Top["iter129_bulkformula"] .- Δwθ0_Top["SUM"]

# Δwθ0["eps"]
vars = vcat(vars,"SUM")
plot_labels["eps"] = "ϵ"
exp_colors["eps"] = "grey"

fig, axs = plt.subplots(1, 2, figsize = ( 10, 5), sharey = true)
axs[1].set_title("Effect of Control Adjustments on \n " * L"\mathbf{A}^\mathcal{B}(\Delta w^{res}, \theta_0)}")
axs[2].set_title("Effect of Control Adjustments on \n " * L"\mathbf{A}^\mathcal{U}(\Delta w^{res}, \theta_0)}")

# axs.set_xlabel("time", fontweight = "bold")
vars =   [ "only_init", "only_kappa", "only_wind", "only_buoyancy", "eps", "iter129_bulkformula"]

# axs.set_ylabel(" [cK]", fontweight = "bold") 

plot_labels_list = ["Initial Condition\nAdjustment Effect", "Mixing Parameter\nAdjustment Effect", 
"Wind Stress\nAdjustment Effect ", "Buoyancy Forcing\nAdjustment Effect ", "Non-Linear\nAdjustment Effects", "All Control\nAdjustments"]
for (i, expname) in enumerate(vars)
    println(expname)
    lw = 2.5
    println(round(get_trend(Δwθ0_Bot[expname]), digits = 3))
    axs[1].plot(tecco, Δwθ0_Bot[expname], label = plot_labels_list[i], 
    linewidth = lw, color = exp_colors[expname])
    axs[1].axhline(0, c = "k", linestyle = "--")

    axs[2].plot(tecco, Δwθ0_Top[expname], label = plot_labels_list[i], 
    linewidth = lw, color = exp_colors[expname])
    axs[2].axhline(0, c = "k", linestyle = "--")

end

var_explained(Δwθ0_Bot["iter129_bulkformula"], Δwθ0_Bot["only_wind"])
var_explained(Δwθ0_Bot["iter129_bulkformula"], Δwθ0_Bot["only_kappa"])
var_explained(Δwθ0_Bot["iter129_bulkformula"], Δwθ0_Bot["only_init"])
var_explained(Δwθ0_Bot["iter129_bulkformula"], Δwθ0_Bot["eps"])

var_explained(Δwθ0_Top["iter129_bulkformula"], Δwθ0_Top["only_wind"])
var_explained(Δwθ0_Top["iter129_bulkformula"], Δwθ0_Top["only_kappa"])
var_explained(Δwθ0_Top["iter129_bulkformula"], Δwθ0_Top["only_init"])
var_explained(Δwθ0_Top["iter129_bulkformula"], Δwθ0_Top["eps"])



# var_explained(Δwθ0_Top["iter129_bulkformula"], Δwθ0_Top["only_buoyancy"])

[a.grid for a in axs]
fig
fig.subplots_adjust(hspace = 0.5)
axs[1].legend(frameon = false, ncols = 3, loc = "lower center",  bbox_to_anchor=(1.15, -0.35), fontsize = 11)
fig
# fig.savefig(plotsdir("native/paper_figures/2.HeatBudgetVertΔWDecomp_Adjustments_all.png"), bbox_inches = "tight", dpi = 400)

fig


