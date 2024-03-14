include("../../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools
import PyPlot as plt
using DSP
using LinearAlgebra, PyPlot, Statistics

using SavitzkyGolay

include(srcdir("plot_and_dir_config.jl"))
sns.set_theme(context = "paper", style = "ticks",
            palette = colors, rc = custom_params);
@pyimport matplotlib.patches as patches

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)
region = "NPAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)

cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)
ΔV = lateral_sum(cell_volumes)

lvls = findall( -3000 .<= -z[:].<= -2000)
mid_depths(x) = vec(sum(Float32, x[lvls, :] .* ΔV[lvls], dims = 1) / sum(Float32, ΔV[lvls]))

tecco = 1992+1/24:1/12:2018; nz = 50
get_datafiles(expname, key) = filter(x -> occursin("data",x),searchdir(diagpath[expname],key) )

adjust_exps =  jldopen(datadir(region * "_temperature_sens_exps.jld2"))["adjust_exps"]

adjust_exps["Difference"]= mid_depths(adjust_exps["iter129_bulkformula"] .- adjust_exps["iter0_bulkformula"])
adjust_exps["only_init"] = mid_depths(adjust_exps["only_init"] .- adjust_exps["iter0_bulkformula"])
adjust_exps["only_wind"] = mid_depths(adjust_exps["only_wind"] .- adjust_exps["iter0_bulkformula"])
adjust_exps["only_kappa"] = mid_depths(adjust_exps["only_kappa"] .- adjust_exps["iter0_bulkformula"])
adjust_exps["Difference"] = adjust_exps["only_init"] .+ adjust_exps["only_wind"] .+ adjust_exps["only_kappa"]
lw = 2

pad_arr(x) = vcat(reverse(x), x, reverse(x))
unpad_arr(x, arr_len) = x[arr_len+1:2*arr_len]
function low_pass(signal)
    nt = length(signal)
    signal_mean = mean(signal)
    ff = digitalfilter(Lowpass(1/(1*12), fs = 1),Butterworth(4))
    filtered_signal = filtfilt(ff, pad_arr(signal .- signal_mean))
    return unpad_arr(filtered_signal, nt) .+ signal_mean
end

using NLsolve, LsqFit

expnams = ["only_init", "only_wind", "only_kappa"]
coefs = zeros(3, 5)
fits = []
for (i, expt) in enumerate(expnams)
    data = low_pass(adjust_exps[expt])
    std_var = std(1  .* data)
    data = data ./ std_var
    t = collect(0:(length(data) - 1)) ./ 12
    lb = [-Inf, -Inf, -Inf, -Inf, -Inf]
    ub = [0, Inf, Inf, Inf, Inf]

    model(x, p) = p[4] .* exp.(p[1] .* x) .* cos.((x .* p[2]) .+ p[3]) .+ p[5]
    
    fit = curve_fit(model, t, data, [-0.01, 0.1, 0.1, 0.0, 0.001], lower=lb, upper=ub)
    coefs[i, :] = fit.param
    coefs[i, 4] = coefs[i, 4] .* std_var
    coefs[i, end] = coefs[i, 5] .* std_var
    push!(fits, fit)
end
coefs[:, 3] .= sign.(coefs[:, 2]) .* coefs[:, 3]
coefs[:, 2] .= sign.(coefs[:, 2]) .* coefs[:, 2]
t = collect(0:311) ./ 12
t_extrap = collect(0:1500) ./ 12

model(x, p) = p[4] .* exp.(p[1] .* x) .* cos.((x .* p[2]) .+ p[3]) .+ p[5]
using StatsBase

expnams = ["only_init", "only_wind", "only_kappa",]
fig, axs = plt.subplots(4, 1, sharey = true, sharex = false, figsize = (7, 7))
adjust_exps["iter129_bulkformula"] = adjust_exps["Difference"]
for (i, expt) in enumerate(expnams)
    ax = axs[i]
    ax.plot(t_extrap, 100 .* model(t_extrap, coefs[i, :]), 
            linewidth = lw, color = exp_colors[expt], alpha = 0.4)
    ax.plot(t, 100 .* adjust_exps[expt], linewidth = lw, color = exp_colors[expt], alpha = 1)
end

ax = axs[end]; expt = "iter129_bulkformula"
ax.plot(t_extrap, 100 .* (model(t_extrap, coefs[1, :]) .+ 
model(t_extrap, coefs[2, :]) .+ model(t_extrap, coefs[3, :])), 
linewidth = lw, color = exp_colors[expt], alpha = 0.4)
ax.plot(t, 100 .* adjust_exps[expt], linewidth = lw, color = exp_colors[expt], alpha = 1)
fig

fsize = 10
[a.set_ylabel(L"\theta'" * " [cK]") for a in axs[:]]
[a.spines["top"].set_visible(false) for a in axs[:]]
[a.spines["bottom"].set_visible(false) for a in axs[:]]
[a.spines["right"].set_visible(false) for a in axs[:]]
[a.set_ylim(-1.5, 1.5) for a in axs[:]]
[a.set_yticks(-1:0.5:1) for a in axs[:]]
[a.set_xticks([]) for a in axs[1:3]]

axs[1].annotate("Initial Condition\nAdjustment Response", (0.8, 0.65), fontsize = fsize, 
xycoords="axes fraction", ha="center", fontweight = "normal", color = exp_colors["only_init"])

axs[2].annotate("Mixing Parameter\nAdjustment Response", (0.8, 0.65), fontsize = fsize, 
xycoords="axes fraction", ha="center", fontweight = "normal", color = exp_colors["only_kappa"])

axs[3].annotate("Wind Stress\nAdjustment Response", (0.8, 0.65), fontsize = fsize, 
xycoords="axes fraction", ha="center", fontweight = "normal", color = exp_colors["only_wind"])

axs[4].annotate("All\nAdjustment Response", (0.8, 0.65), fontsize = fsize, 
xycoords="axes fraction", ha="center", fontweight = "normal", color = exp_colors["iter129_bulkformula"])

axs[end].spines["bottom"].set_visible(true)
fig.subplots_adjust(hspace = 0.0, wspace = 0.27)
axs[1].set_title("Projected Response of Mid-Depth North Pacific \n Temperature to ECCO Control Adjustments")
axs[end].set_xlabel("Years since 1992")
fig

P = inv.(coefs[:, 2]) .* 2π #period of oscillation

coefs