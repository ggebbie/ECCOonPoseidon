include("../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools
import PyPlot as plt
@pyimport cmocean.cm as cmos

include(srcdir("plot_and_dir_config.jl"))
@pyimport matplotlib.patches as patches

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)

cell_depths = get_cell_thickness(ocean_mask, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)

lvls = findall( -3000 .<= -z[:].<= -2000); suffix = "2to3"
tecco = 1992+1/24:1/12:2018; nz = 50

#load trends
fname = datadir("native/_THETA_spatial_trend_" * suffix * ".jld2")
theta_trends = jldopen(fname, "r")["theta_trends"]

bounds = 0.2 * 10

# theta_trends["SUM"] = theta_trends["Kappa"] .+ theta_trends["Initial"] .+ theta_trends["Forcing"] .+ (-2 .* theta_trends["CTRL"])
# plot_labels["SUM"] = "SUM"

lin_exps = ["CTRL", "Initial", "Kappa", "Forcing", "FULL", "SUM"]; nexps = length(lin_exps)
plot_labels = ["CTRL (Iteration 0)", "INIT", "MIXING", "FORCING", "FULL (Iteration 129)", "SUM"]
plot_labels = Dict(lin_exps[i] => plot_labels[i] for i in 1:nexps)


fig, axs = plt.subplots(2, 3, figsize=(17,12), subplot_kw=Dict("projection"=> proj0))
CF = Any[]

lin_exps = ["Initial", "CTRL", "Kappa", "FULL", "Forcing", "SUM"]; nexps = length(lin_exps)

for (i, exp) in enumerate(lin_exps)
    ax = axs[i]
    data = theta_trends[exp].* 100 * 10
    for ff = 1:5
        cf = ax.pcolormesh(λ[ff], ϕ[ff],  data[ff], transform=projPC, 
    cmap = cmos.balance, vmin = -bounds, vmax = bounds)   
    push!(CF, cf)

    end
    ax.coastlines(resolution="110m")
    ax.set_extent((120, 295, -70, 56),crs=projPC)
    ax.set_title(plot_labels[exp])
    gl = ax.gridlines(crs=projPC, draw_labels=true,
                        linewidth=2, color="gray", alpha=0, linestyle="--")
    gl.top_labels = false
    gl.right_labels = false

end
fig.colorbar(CF[1], ax = axs[:], orientation = "horizontal", extend = "both", 
fraction = 0.03, label = "cK per decade")
fig.savefig(plotsdir("native/sensitivity_exps/theta_trends.png"))


theta_trends_effect = jldopen(fname, "r")["theta_trends"]
theta_trends_effect["Initial"] = theta_trends["Initial"] .- theta_trends["CTRL"]
theta_trends_effect["Kappa"] = theta_trends["Kappa"] .- theta_trends["CTRL"]
theta_trends_effect["Forcing"] = theta_trends["Forcing"] .- theta_trends["CTRL"]

lin_exps = ["CTRL", "Initial", "Kappa", "Forcing", "FULL", "SUM"]; nexps = length(lin_exps)
plot_labels_effect = ["CTRL (Iteration 0)", "INIT Effect", "MIXING Effect", "FORCING Effect", "FULL (Iteration 129)", "SUM"]
plot_labels_effect = Dict(lin_exps[i] => plot_labels_effect[i] for i in 1:nexps)
lin_exps = ["Initial", "CTRL", "Kappa", "FULL", "Forcing", "SUM"]; nexps = length(lin_exps)


fig, axs = plt.subplots(2, 3, figsize=(17,12), subplot_kw=Dict("projection"=> proj0))
CF = Any[]

for (i, exp) in enumerate(lin_exps)
    ax = axs[i]
    data = (theta_trends_effect[exp]).* 100 * 10 
    for ff = 1:5
        cf = ax.pcolormesh(λ[ff], ϕ[ff],  data[ff], transform=projPC, 
    cmap = cmos.balance, vmin = -bounds, vmax = bounds)   
    push!(CF, cf)
    end
    ax.coastlines(resolution="110m")
    ax.set_extent((120, 295, -70, 56),crs=projPC)
    ax.set_title(plot_labels_effect[exp])
    gl = ax.gridlines(crs=projPC, draw_labels=true,
                        linewidth=2, color="gray", alpha=0, linestyle="--")
    gl.top_labels = false
    gl.right_labels = false

end
fig
fig.colorbar(CF[1], ax = axs[:], orientation = "horizontal", extend = "both", 
fraction = 0.03, label = "cK per decade")
fig.savefig(plotsdir("native/sensitivity_exps/theta_trends_effect.png"))
fig


