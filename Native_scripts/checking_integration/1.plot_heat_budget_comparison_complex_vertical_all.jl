include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson,
    BenchmarkTools, LaTeXStrings, PyCall, DataFrames
import NaNMath as nm
import PyPlot as plt
import NumericalIntegration

include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");

integrate(start, x) = cumsum([start, x...])[1:end-1]

function adams_bashforth_2nd_order(f, y0, t0, n_steps)
    h = 2.628e+6
    # Initialize arrays to store results
    t_values = zeros(n_steps + 1)
    y_values = zeros(n_steps + 1)
    
    # Set initial conditions
    t_values[1] = t0
    y_values[1] = y0
    
    # Use Euler's method to compute the first step
    y_values[2] = y0 + h * f[1]
    t_values[2] = t0 + h
    
    # Iterate using the Adams-Bashforth formula
    for n in 2:n_steps
        y_n = y_values[n]
        
        y_n_plus_1 = y_n + (h / 2) * (3 * f[n] - f[n-1])
        
        y_values[n + 1] = y_n_plus_1
    end
    
    return y_values
end

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

expname =   "iter129_bulkformula"


fig, axs = plt.subplots(figsize = ( 6, 4))
axs.set_xlabel("time", fontweight = "bold")
axs.set_title("Decomposition of " * L"\mathbf{A}_Z (w_{129}^{res}, \theta_{129})")
lw = 2.0

axs.plot(tecco, wθ_resid["iter129_bulkformula"], label = L"\mathbf{A}_Z (w_{129}^{res}, \theta_{129})", 
linewidth = lw, color = exp_colors[expname])
axs.plot(tecco, wθ_resid["iter0_bulkformula"], label = L"\mathbf{A}_Z (w_0^{res}, \theta_0)", 
linewidth = lw, color = exp_colors["iter0_bulkformula"])

tm1 = Δwθ[expname] .-ΔwΔθ[expname]
axs.plot(tecco, Δwθ[expname] .-ΔwΔθ[expname] , label = L"\mathbf{A}_Z (\Delta w^{res}, \theta_0)", 
linewidth = 3, color = 0.5 .* exp_colors[expname], alpha = 0.8)
tm2 = wΔθ[expname] .-ΔwΔθ[expname]

axs.plot(tecco, wΔθ[expname] .-ΔwΔθ[expname] , label = L"\mathbf{A}_Z (w_0^{res}, \Delta \theta)",
linewidth = lw, color = exp_colors[expname], alpha = 0.6)
tmp = (Δwθ[expname] .-ΔwΔθ[expname]) 
tmp = (wθ_resid["iter129_bulkformula"] .- wθ_resid["iter0_bulkformula"]) .- tmp
tmp = ΔwΔθ[expname]

axs.plot(tecco, tmp , label = L"\mathbf{A}_Z (  \Delta w^{res}, \Delta \theta)", 
linewidth = lw, color = exp_colors[expname], alpha = 0.3)
axs.legend(frameon = false, ncols = 1, loc = "center right",  bbox_to_anchor=(1.35, 0.5))
axs.axhline(0, c = "k", linestyle = "--", alpha = 0.2)

fig
axs.set_ylabel("[cK]", fontweight = "bold")
# fig.subplots_adjust(hspace = 0.5)
axs.grid()
fig
fig

fname = datadir("native/" * "iter129_bulkformula" * region * "_THETA_budget_Wθ_eul_bol_ΔDecomp" * suffix * ".jld2")
vars = jldopen(fname)["dθ"]

vars =   ["iter0_bulkformula", "only_init", "only_kappa", "only_wind", "only_buoyancy", "iter129_bulkformula"]

wθ_resid = Dict(); wθ_eul = Dict(); wθ_bol = Dict()

Δwθ = Dict()
wΔθ = Dict()
ΔwΔθ = Dict()

for (i, expname) in enumerate(vars)
    fname = datadir("native/" * expname * region * "_THETA_budget_Wθ_eul_bol_ΔDecomp" * suffix * ".jld2")
    print(expname)
    vars = jldopen(fname)["dθ"]

    wθfull = vars["wθBot"] .- vars["wθTop"]
    wθ_resid[expname] = 100 .* adams_bashforth_2nd_order(wθfull, 0, 0, 311)
    Δwθ[expname] = 100 .* adams_bashforth_2nd_order(vars["ΔW_θBot"] .- vars["ΔW_θTop"], 0, 0, 311)
    wΔθ[expname] = 100 .* adams_bashforth_2nd_order(vars["W_ΔθBot"] .- vars["W_ΔθTop"], 0, 0, 311)
    ΔwΔθ[expname] = 100 .* adams_bashforth_2nd_order(vars["ΔW_ΔθBot"] .- vars["ΔW_ΔθTop"], 0, 0, 311)
end

expname =   "iter129_bulkformula"
fig, axs = plt.subplots(figsize = ( 6, 4))
axs.set_xlabel("time", fontweight = "bold")
axs.set_title("Decomposition of " * L"\mathbf{A}_Z (w_{129}^{res}, \theta_{129})")
lw = 2.0

axs.plot(tecco, wθ_resid["iter129_bulkformula"], label = L"\mathbf{A}_Z (w_{129}^{res}, \theta_{129})", 
linewidth = lw, color = exp_colors[expname])
axs.plot(tecco, wθ_resid["iter0_bulkformula"], label = L"\mathbf{A}_Z (w_0^{res}, \theta_0)", 
linewidth = lw, color = exp_colors["iter0_bulkformula"])

tm1 = Δwθ[expname] .-ΔwΔθ[expname]
axs.plot(tecco, Δwθ[expname] .-ΔwΔθ[expname] , label = L"\mathbf{A}_Z (\Delta w^{res}, \theta_0)", 
linewidth = 3, color = 0.5 .* exp_colors[expname], alpha = 0.8)
tm2 = wΔθ[expname] .-ΔwΔθ[expname]

axs.plot(tecco, wΔθ[expname] .-ΔwΔθ[expname] , label = L"\mathbf{A}_Z (w_0^{res}, \Delta \theta)",
linewidth = lw, color = exp_colors[expname], alpha = 0.6)
tmp = (Δwθ[expname] .-ΔwΔθ[expname]) 

# tmp = (wθ_resid["iter129_bulkformula"] .- wθ_resid["iter0_bulkformula"]) .- tmp
tmp = ΔwΔθ[expname]
axs.plot(tecco, tmp , label = L"\mathbf{A}_Z (  \Delta w^{res}, \Delta \theta)", 
linewidth = lw, color = exp_colors[expname], alpha = 0.3)
axs.legend(frameon = false, ncols = 1, loc = "center right",  bbox_to_anchor=(1.35, 0.5))
axs.axhline(0, c = "k", linestyle = "--", alpha = 0.2)

fig
axs.set_ylabel("[cK]", fontweight = "bold")
# fig.subplots_adjust(hspace = 0.5)
axs.grid()
fig
fig