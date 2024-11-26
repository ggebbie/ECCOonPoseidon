include("../../src/intro.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
JLD2, Printf, LaTeXStrings, PyPlot, PyCall

cmo = pyimport("cmocean.cm");
@pyimport matplotlib.patches as patches

include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

region = "NPAC"
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region)

area_mask = area .* PAC_msk
runpath,diagpath = listexperiments(exprootdir());

tecco= 1992+1/24:1/12:2018 # ecco years
E,F = trend_matrices(tecco)
sig2grid = sigma2grid("NPAC")
nσ = length(sig1grid)
TSroot = "p_on_sigma2" 

P_dict = Dict();
vars = ["only_init", "only_kappa", "only_sfc",  "iter0_bulkformula", "iter129_bulkformula"]
nvars = length(vars)
for expname in vars
    Pσ = jldopen(datadir(expname * region * "_AVG_P_sigma2.jld2"))["P"]
    P_dict[expname] = Pσ .* 1

    #interpolate to z-axis for clarity
end

fig, axs = plt.subplots(1, 5, figsize = (20, 10))
vmax = 20
clevels = -vmax:2:vmax
baseline = P_dict["iter0_bulkformula"]
P_dict["SUM"] = P_dict["only_init"] .+ P_dict["only_kappa"] .+ P_dict["only_sfc"] .- (2 .* P_dict["iter0_bulkformula"])
for (i, expt) in enumerate(["only_init", "only_kappa", "only_sfc", "iter129_bulkformula", "SUM"])
    ax = axs[:][i]
    P = P_dict[expt] .- baseline; P_mean = mean(P, dims = 2)
    P_anom = P .- P_mean
    cm = ax.contourf(tecco, sig2grid, P .- P_mean, vmin = -vmax, vmax = vmax, cmap = cmo.balance, levels = clevels, extend = "both")
    fig.colorbar(cm, orientation = "horizontal", label = "Pa Anomaly")
    ax.axhline(36.7, color = "black", zorder = 10, alpha = 1, linewidth = 1)
    ax.axhline(37, color = "black", zorder = 10, alpha = 1, linewidth = 1)
    ax.set_title(expt)
    println(expt)
    println( P_anom[end-5, end] - P_anom[end-5, 1])
    ax.set_ylim(35.0, 37.1); ax.invert_yaxis()

end
fig
# cm = ax[1, 2].contourf(tecco, sig2grid, P .- P_mean, vmin = -vmax, vmax = vmax, cmap = cmo.balance, levels = clevels, extend = "both")
# fig.colorbar(cm, orientation = "horizontal", label = "Pa Anomaly")
fig

fig, ax = plt.subplots(2, 2, figsize = (10, 10), sharey = true)

fname = datadir("native/" * region * "_isopycnal_timeseries.jld2")
σz_mean = jldopen(fname)["σz_mean"]
vmin, vmax = extrema(sig2grid[34:end])

ax[2, 1].contourf(tecco, z, σz_mean["iter0_bulkformula"]', cmap = cmo.dense, vmin =  vmin, vmax = vmax + 0.02, 
levels = sig2grid[1:3:end], extend = "both")
cl = ax[2, 1].contour(tecco, z, σz_mean["iter0_bulkformula"]', colors = "black", linewidths = 4, 
levels = sig2grid[1:3:end])
ax[2, 1].clabel(cl)

ax[2, 2].contourf(tecco, z, σz_mean["iter129_bulkformula"]', cmap = cmo.dense, vmin =  vmin, vmax = vmax + 0.02, 
levels = sig2grid[1:3:end], extend = "both")
cl = ax[2, 2].contour(tecco, z, σz_mean["iter129_bulkformula"]', colors = "black", linewidths = 4, 
levels = sig2grid[1:3:end])
ax[2, 2].clabel(cl)
ax[2, 2].set_ylim(100, 3200)
ax[2, 2].invert_yaxis()

fig
