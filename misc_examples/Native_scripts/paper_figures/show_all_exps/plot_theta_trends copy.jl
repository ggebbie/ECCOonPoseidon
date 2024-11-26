include("../../../src/intro.jl")

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

bounds = 0.14 * 100
levels = collect(-bounds:2:bounds)
lin_exps = ["FULL", "Initial", "Diff"]; nexps = length(lin_exps)
theta_trends["Diff"] = theta_trends["FULL"] .- theta_trends["Initial"]
plot_labels["Diff"] = "Difference"

expts = ["iter129_bulkformula", "iter0_bulkformula", "Diff"]

fig, axs = plt.subplots(1, 3, figsize=(14,10), subplot_kw=Dict("projection"=> proj0))
CF = Any[]
for (i, exp) in enumerate(lin_exps)
    ax = axs[i]

    ax.coastlines(resolution="110m")
    ax.set_extent((110, 295, -65, 65),crs=projPC)
    ax.set_title(plot_labels[expts[i]])
    gl = ax.gridlines(crs=projPC, draw_labels=true,
                        linewidth=2, color="gray", alpha=0, linestyle="--")
    gl.top_labels = false
    gl.right_labels = false 

    data = theta_trends[exp].* 100 * 100

    x = vcat([λ.f[ff][:] for ff = 1:5]...);  x[x .< 0] .+= 360
    y = vcat([ϕ.f[ff][:] for ff = 1:5]...)
    z = vcat([data.f[ff][:] for ff = 1:5]...)
    triangle = tri.Triangulation(x, y)
    mask = isnan.(z[triangle.triangles[:, 2] .+ 1] .+ z[triangle.triangles[:, 1] .+ 1] .+ z[triangle.triangles[:, 3] .+ 1])
    triangle.set_mask(mask)
    z[isnan.(z)] .= 0.0
    tcf = ax.tricontourf(triangle, z, transform=projPC,     
    cmap = cmos.balance, levels = levels, 
    vmin = -bounds, vmax = bounds, extend = "both")
    CS = ax.tricontour(triangle, z, transform=projPC,     
    color = "k", levels = [0])
    ax.clabel(CS, CS.levels, fontsize = 12, inline = true, inline_spacing = 15)
    push!(CF, tcf)
   
end
cbar = fig.colorbar(CF[1], ax = axs[:], orientation = "horizontal", extend = "both", 
fraction = 0.027, label = "cK per century", pad = 0.06)
# cbar.set_ticks(levels)
fig_labs = uppercase.(["a", "b", "c", "d", "e"])
for (i, a) in enumerate(axs)
    a.annotate(fig_labs[i], (0.03, 0.87), fontsize = 20, 
    xycoords="axes fraction", fontweight = "bold")
end
fig.suptitle("Mid-Depth Temperature Trends between 1992 and 2017", y = 0.47)

fig
fig.savefig(plotsdir("native/paper_figures/1.theta_trends.png"), bbox_inches = "tight", dpi = 800)