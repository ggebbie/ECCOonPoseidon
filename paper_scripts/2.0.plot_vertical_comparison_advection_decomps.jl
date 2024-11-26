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



for (i, expname) in enumerate(["mean_tau_yesadjust_redo_bf", "mean_tau_noadjust_redo_bf"])
    fname = datadir("native/" * expname * region * "_THETA_budget_Wθ_eul_bol_ΔDecomp_climatological" * suffix * ".jld2")
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

exp_colors["mean_tau_noadjust_redo_bf"] = exp_colors["iter0_bulkformula"] 
exp_colors["mean_tau_yesadjust_redo_bf"] = exp_colors["only_wind"] 

sns.set_style("darkgrid", Dict("axes.facecolor" => ".95", "grid.color" => "0.5"))

vars =   [ "only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]
E,F = trend_matrices(tecco)
get_trend(x) = (F * x)[2]
var_explained(y, x) = 100 * ( 1 - (var(y - x) / var(y))) 
fig, ax = plt.subplots(1, 2, figsize = ( 12, 6), sharey = true)
expname =   "mean_tau_yesadjust_redo_bf"
lw = 1.5

wθ_resid_list = [wθ_resid_bot, wθ_resid_top]
Δwθ_list = [Δwθ_bot, Δwθ_top]
wΔθ_list = [wΔθ_bot, wΔθ_top]
ΔwΔθ_list = [ΔwΔθ_bot, ΔwΔθ_top]
ax[1].set_title("Decomposition of Bottom Boundary\n Advective Heat Flux Response to Optimization\n" * L"\Delta\mathbf{A}_{\mathcal{B}}")
ax[2].set_title("Decomposition of Upper Boundary\nAdvective Heat Flux Response to Optimization\n" * L"\Delta\mathbf{A}_{\mathcal{U}}")
fig
var_name = [L"\mathbf{A}^\mathcal{B}", L"\mathbf{A}^\mathcal{U}"]
MAE(x, y) = mean(abs.(y .- x))
for (i, axs) in enumerate(ax)
    wθ_resid = wθ_resid_list[i]
    Δwθ = Δwθ_list[i]
    wΔθ = wΔθ_list[i]
    wΔθ = ΔwΔθ_list[i]
    axs.set_xlabel("time", fontweight = "bold")
    ΔA = wθ_resid["mean_tau_yesadjust_redo_bf"] .- wθ_resid["mean_tau_noadjust_redo_bf"]
    axs.plot(tecco, ΔA, label = L"\Delta" * var_name[i], 
    linewidth = 2.5, color = exp_colors["mean_tau_yesadjust_redo_bf"])

    # ΔA_Δw = Δwθ[expname] .-ΔwΔθ[expname]
    # println((F * tmp1)[2])
    axs.plot(tecco, Δwθ["mean_tau_yesadjust_redo_bf"] , label =  var_name[i] * L"(\Delta w^{res}, \theta_0)", 
    linewidth = 1.5, color = "k", alpha = 0.8)
    # ΔA_Δθ = wΔθ[expname] .-ΔwΔθ[expname]
# 
    axs.plot(tecco, wΔθ["mean_tau_yesadjust_redo_bf"]  , label = var_name[i] * L"(w_0^{res}, \Delta \theta)", 
    linewidth = lw, color = "blueviolet", alpha = 0.9)

    # tmp = (ΔA_Δw) .+ (ΔA_Δθ)
    # ΔA_Δw_Δθ = ΔA .- tmp

    axs.plot(tecco, wΔθ["mean_tau_yesadjust_redo_bf"]  , label = var_name[i] * L"(  \Delta w^{res}, \Delta \theta)", 
    linewidth = lw, color = "hotpink", alpha = 0.9, linestyle = "--")
    axs.axhline(0, c = "k", linestyle = "--", alpha = 0.2)
    axs.grid()
    leg = axs.legend(frameon = false, ncols = 1, loc = "lower left", fontsize=12)
    # println("Var. Explained")
    # println(var_explained(ΔA, ΔA_Δw))
    # println(var_explained(ΔA, ΔA_Δθ))
    # println(var_explained(ΔA, ΔA_Δw_Δθ))

    # println("MAE")
    # println(MAE(ΔA, ΔA_Δw))
    # println(MAE(ΔA, ΔA_Δθ))
    # println(MAE(ΔA, ΔA_Δw_Δθ))
end
ax[1].set_ylabel("[cK]", fontweight = "bold")
fig_labs = uppercase.(["a", "b", "c", "d", "e", "f"])
for (i, a) in enumerate(ax)
    a.annotate(fig_labs[i], (0.02, 0.92), fontsize = 20, color = "black", 
    xycoords="axes fraction", fontweight = "bold")
end
fig