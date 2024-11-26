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
@pyimport matplotlib.gridspec as gridspec

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
plot_labels["Diff"] = "Iteration 129 minus Iteration 0"

expts = ["iter129_bulkformula", "iter0_bulkformula", "Diff"]

fig = plt.figure(figsize=(10, 8))

# Define the grid spec layout with two rows and three columns
gs = gridspec.GridSpec(2, 3, width_ratios=[1, 1, 0.05])

# Add the first subplot (1st row, 2nd column)
ax1 = fig.add_subplot(py"$gs[0, 0]", projection= proj0)

# Add the second subplot (1st row, 3rd column)
ax2 = fig.add_subplot(py"$gs[0, 1]", projection= proj0)

# Add the third subplot (2nd row, 2nd and 3rd columns)
ax3 = fig.add_subplot(py"$gs[1, :2]", projection= proj0)
axs = [ax1, ax2, ax3]
CF = Any[]
for (i, exp) in enumerate(lin_exps)
    ax = axs[i]

    ax.coastlines(resolution="110m", color = "#949494")
    ax.set_extent((110, 295, -65, 65),crs=projPC)
    ax.set_title(plot_labels[expts[i]])
    gl = ax.gridlines(crs=projPC, draw_labels=true,
                        linewidth=2, color="gray", alpha=0, linestyle="--")
    gl.top_labels = false
    gl.right_labels = false 
    if i == 2
        gl.left_labels = false 
    end
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
    ax.add_feature( ECCOonPoseidon.cartopy.feature.LAND, facecolor="#949494")
    ax.clabel(CS, CS.levels, fontsize = 12, inline = true, inline_spacing = 15)
    push!(CF, tcf)
   
end

cbar_ax = fig.add_subplot(py"$gs[1, 2]")
cbar = fig.colorbar(CF[1],  cax=cbar_ax, orientation = "vertical", extend = "both", 
fraction = 0.027, label = "cK per century", pad = 0.02)
cbar_ax.set_position([ax3.get_position().x1 + 0.07, ax3.get_position().y0, 
                      0.02, ax3.get_position().height])
# cbar.set_ticks(levels)
fig_labs = uppercase.(["a", "b", "c", "d", "e"])
for (i, a) in enumerate(axs)
    a.annotate(fig_labs[i], (0.03, 0.87), fontsize = 20, 
    xycoords="axes fraction", fontweight = "bold")
end
fig.suptitle("Mid-Depth Temperature Trends (1992 and 2017)", y = 0.94)
fig.savefig(plotsdir("wind_rewrite/0.1.theta_trends.png"), bbox_inches = "tight", dpi = 400)

fig