using Pkg
Pkg.activate(".")
include("../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings
import NaNMath as nm
using PyCall
import PyPlot as plt
@pyimport cmocean.cm as cmos
@pyimport matplotlib.ticker as ticker
@pyimport matplotlib.colors as colors_plt
@pyimport matplotlib.tri as tri
@pyimport matplotlib.gridspec as gridspec

include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

include(srcdir("plot_and_dir_config.jl"))
 
ocean_mask = wet_pts(Γ)
region = "PAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "false", include_bering = true)

tecco = 1992+1/24:1/12:2018

cell_depths = get_cell_thickness(ocean_mask, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)


vars =  ["only_init", "only_kappa", "only_sfc", "iter129_bulkformula",  "iter0_bulkformula"]
vars =  ["only_buoyancy", "only_wind"]

uplvl = -2.0e3; botlvl = -3.0e3; suffix = "2to3"
lvls = findall( botlvl .<= -z[:].<= uplvl)

diagpath["mean_tau_noadjust_redo_bf"] = vastdiagdir("seasonalclimatology", "run_only_clim_iter0_tau_bf")
diagpath["mean_tau_yesadjust_redo_bf"] = vastdiagdir("seasonalclimatology", "run_only_clim_iter129_tau_bf2")

p₀ = 2000; nz = 50
P = MeshArray(γ,Float32, nz)
for ijk in eachindex(P)
    P[ijk] .= pstdz[ijk[2]]
end

θ_final_dict = Dict()
@time for expname in ["iter0_bulkformula", "only_wind", "mean_tau_noadjust_redo_bf", "mean_tau_yesadjust_redo_bf"]
    println(expname)

    ΔV = vertical_sum(cell_volumes[:, lvls]); 
    ΔV[findall(ΔV .== 0.0)] .= NaN

    filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
    datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    σ = MeshArray(γ,Float32,nz)

    fnameθ = datafilelist_θ[312]
    θS = γ.read(diagpath[expname]*fnameθ,MeshArray(γ,Float32,100))
    θz = θS[:, 1:50]
    Sz = θS[:, 51:end]

    for ijk in eachindex(P)
        σtemp = densityJMD95.(θz.f[ijk],Sz.f[ijk], P[ijk], p₀) #EOS from MITGCM 
        σ.f[ijk] .= σtemp .- 1000
    end 

    θ_final_dict[expname] = vertical_sum(θ[:, lvls] .* cell_volumes[:, lvls]); 
    θ_final_dict[expname] = θ_final_dict[expname] ./ ΔV
end

θ_final_dict["Difference"] = θ_final_dict["mean_tau_yesadjust_redo_bf"] .- θ_final_dict["mean_tau_noadjust_redo_bf"]

fig, axs = plt.subplots(1, 3, figsize=(15, 8), subplot_kw = Dict("projection"=>proj0))
CF = Any[]
boundss = [2, 2, 2.5]
mult = [1, 1, 100]
spacing = [0.5, 0.2, 0.5]

for (i, exp) in enumerate(["mean_tau_yesadjust_redo_bf", "mean_tau_noadjust_redo_bf", "Difference"])
    bounds = boundss[i]
    levels = collect(-bounds:spacing[i]:bounds)
    ax = axs[i]

    ax.coastlines(resolution="110m", color = "#949494")
    ax.set_extent((110, 295, -65, 65),crs=projPC)
    # ax.set_title(plot_labels[expts[i]])
    gl = ax.gridlines(crs=projPC, draw_labels=true,
                        linewidth=2, color="gray", alpha=0, linestyle="--")
    gl.top_labels = false
    gl.right_labels = false 
    if (i == 2) || (i == 3)
        gl.left_labels = false 
    end
    data = mult[i] .* θ_final_dict[exp]

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
cbar = fig.colorbar(CF[1],  ax=axs[1:2], orientation = "horizontal", extend = "both", 
fraction = 0.027, label = "deg C", pad = 0.1)

cbar = fig.colorbar(CF[3],  ax=axs[3], orientation = "horizontal", extend = "both", 
fraction = 0.027, label = "centikelvin", pad = 0.1)

fig
# cbar_ax.set_position([ax3.get_position().x1 + 0.07, ax3.get_position().y0, 
#                       0.02, ax3.get_position().height])
# cbar.set_ticks(levels)
fig_labs = uppercase.(["a", "b", "c", "d", "e"])
for (i, a) in enumerate(axs)
    a.annotate(fig_labs[i], (0.03, 0.87), fontsize = 20, 
    xycoords="axes fraction", fontweight = "bold")
end
fig.suptitle("Mid-Depth Temperature Trends between 1992 and 2017", y = 0.94)
# fig.savefig(plotsdir("wind_rewrite/0.1.theta_trends.png"), bbox_inches = "tight", dpi = 400)

fig