include("../../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))

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

for (i, expname) in enumerate(vars)
    fname = datadir("native/" * expname * region * "_THETA_budget_Wθ_eul_bol_ΔDecomp" * suffix * ".jld2")
    print(expname)
    vars = jldopen(fname)["dθ"]

    wθfull = vars["wθBot"] .- vars["wθTop"]
    wθ_resid[expname] = wθfull
    Δwθ[expname] = vars["ΔW_θBot"] .- vars["ΔW_θTop"]
    wΔθ[expname] = vars["W_ΔθBot"] .- vars["W_ΔθTop"]
    ΔwΔθ[expname] = vars["ΔW_ΔθBot"] .- vars["ΔW_ΔθTop"]

end

vars =   [ "only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]

[wθ_resid[nexpt] .-= wθ_resid["iter0_bulkformula"] for nexpt in vars]


fig, axs = plt.subplots(3, 1, figsize = (7.5, 21))
[ax.set_xlabel("time", fontweight = "bold") for ax in axs]

vars =   ["only_init", "only_kappa", "only_wind"]

for (i, expname) in enumerate(vars)
    # println(expname)
    lw = 3

    axs[i].plot(tecco, wθ_resid[expname], label = L"(w \frac{\partial}{\partial z} \theta)'", 
    linewidth = lw, color = exp_colors[expname])
    axs[i].plot(tecco, Δwθ[expname], label = L"w' \frac{\partial}{\partial z} \theta", 
    linewidth = lw, color = 0.7 .* exp_colors[expname], alpha = 0.9)
    axs[i].plot(tecco, wΔθ[expname], label = L"w \frac{\partial}{\partial z} \theta'", 
    linewidth = lw, color = exp_colors[expname], alpha = 0.5)
    axs[i].plot(tecco, -ΔwΔθ[expname], label = L"w' \frac{\partial}{\partial z} \theta'", 
    linewidth = lw, color = exp_colors[expname], alpha = 0.2)

    axs[i].axhline(0, c = "k", linestyle = "--")
    axs[i].legend(frameon = false, ncols = 1, loc = "upper left",  bbox_to_anchor=(1, 1))
    # axes[3].legend(markerscale = 0.1, ncols = 1, columnspacing = 0.5, frameon = false)

end
fig
