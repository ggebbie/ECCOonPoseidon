include("../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools
import PyPlot as plt
@pyimport cmocean.cm as cmos
@pyimport matplotlib.ticker as ticker
@pyimport matplotlib.colors as colors_plt
@pyimport matplotlib.tri as tri

include(srcdir("plot_and_dir_config.jl"))
@pyimport matplotlib.patches as patches

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)

cell_depths = get_cell_thickness(ocean_mask, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)

sns.set_theme(context = "notebook", style = "ticks",
              palette = colors, rc = custom_params);
lvls = findall( -3000 .<= -z[:].<= -2000); suffix = "2to3"
tecco = 1992+1/24:1/12:2018; nz = 50

#load trends
fname = datadir("native/_THETA_spatial_trend_" * suffix * ".jld2")
theta_trends = jldopen(fname, "r")["theta_trends"]

bounds = 0.15 * 100

lin_exps = ["FULL", "Initial", "Diff"]; nexps = length(lin_exps)
theta_trends["Diff"] = theta_trends["FULL"] .- theta_trends["Initial"]
plot_labels["Diff"] = "Difference"

expts = ["iter129_bulkformula", "iter0_bulkformula", "Diff"]

levels = ticker.MaxNLocator(nbins=11).tick_values(-bounds, bounds)
norm = colors_plt.BoundaryNorm(levels, ncolors=cmos.balance.N, clip=true)

# pick the desired colormap, sensible levels, and define a normalization
# instance which takes data values and translates those into levels.

fig, axs = plt.subplots(1, 3, figsize=(14,10), subplot_kw=Dict("projection"=> proj0))
CF = Any[]
for (i, exp) in enumerate(lin_exps)
    ax = axs[i]
    data = theta_trends[exp].* 100 * 100
    for ff = 1:5
        cf = ax.pcolormesh(λ[ff], ϕ[ff],  data[ff], transform=projPC, 
        cmap = cmos.balance, norm=norm, rasterized = true, shading = "nearest") 
        push!(CF, cf)
    end
    ax.coastlines(resolution="110m")
    ax.set_extent((120, 295, -45, 65),crs=projPC)
    ax.set_title(plot_labels[expts[i]])
    gl = ax.gridlines(crs=projPC, draw_labels=true,
                        linewidth=2, color="gray", alpha=0, linestyle="--")
    gl.top_labels = false
    gl.right_labels = false

end
fig.colorbar(CF[1], ax = axs[:], orientation = "horizontal", extend = "both", 
fraction = 0.027, label = "cK per century", pad = 0.06)
fig_labs = uppercase.(["a", "b", "c", "d", "e"])
for (i, a) in enumerate(axs)
    a.annotate(fig_labs[i], (0.03, 0.87), fontsize = 23, 
    xycoords="axes fraction", fontweight = "bold")
end
fig.suptitle("Mid-Depth Temperature Trends between 1992 and 2017", y = 0.47)

fig
fig.savefig(plotsdir("native/paper_figures/1.theta_trends.png"), bbox_inches = "tight", dpi = 800)