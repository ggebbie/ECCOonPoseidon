#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings,
Plots
gr()
import NaNMath as nm
using .OHC_helper
import PyPlot as plt
using PyCall

@pyimport seaborn as sns;
@pyimport pandas as pd;
colors =  sns.color_palette("colorblind")[1:4]
labels = ["ECCO", "ECCO Forcing", "ECCO Init. and κ", "CTRL"]
sns.set_theme(context = "poster", style = "ticks", font_scale = 1.2,
              palette = sns.color_palette("colorblind"));
pygui(false)
cm = pyimport("cmocean.cm");colorway = cm.balance;
include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
tecco = 1992+1/24:1/12:2018
runpath,diagpath = listexperiments(exprootdir());

#reshape \lambda for plotting 
for ff ∈ [3,4] 
    λ[ff][λ[ff] .<= 0 ] = λ[ff][λ[ff] .<= 0 ] .+ 360
end
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

#load in latitude mask 
ocean_mask = wet_pts(Γ)
region = "PAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
suffix = "2to3"
uplvl = -2e3; botlvl = -3e3
lvls = findall( botlvl .<= z[:].<= uplvl)

#create volume mask
area = readarea(γ)
cell_depths = get_cell_depths(ocean_mask, ΔzF, Γ.hFacC); 
H = smush(cell_depths[:, lvls]); 
H[findall(H .==0 )] = Inf
inv_depths = 1 ./ H
function depth_average(ds::MeshArray, Δh::MeshArray, H::MeshArray, γ)
    depth_avg = MeshArray(γ,Float32)
    fill!(depth_avg, 0.0)
    nz = size(ds, 2)
    for ff=1:5, k=1:nz
        depth_avg[ff] .+= (ds[ff, k] .* Δh[ff, k]) ./ H[ff]
    end
    return depth_avg
end
#load trend trend_matrices, F is LS estimator
E,F = trend_matrices(tecco)
proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
cf = Vector{Any}(undef ,1)

#load in temperatures
fig, axs = plt.subplots(2, 2, figsize=(15,10), subplot_kw=Dict("projection"=> proj0))
@time for (i, expname) in enumerate(keys(shortnames))
    ax = axs[i]
    filelist = searchdir(diagpath[expname], "state_2d_set1") # first filter for state_3d_set1
    datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
    datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    β = MeshArray(γ,Float32); fill!(β, 0.0)
    nt = length(datafilelist_S)
    @time for tt = 1:nt
        fnameS = datafilelist_S[tt]
        fnameθ = datafilelist_θ[tt]
        sθ = extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_depths)
        sθH = depth_average(sθ[:, lvls], cell_depths[:, lvls], H, γ)
        for ff in 1:5
            β[ff] .+= F[2,tt] .* sθH[ff] 
        end
    end
    #plotting climatology :-)
    bounds = [minimum(β.*PAC_msk), maximum(β.*PAC_msk)]
    b1, b2 = (-maximum(abs.(bounds)), maximum(abs.(bounds)))
    β_mask = deepcopy(β); β_mask[findall(β_mask.==0)] = NaN
    # ax.set_title("Depths in " * region)
    ax.coastlines(resolution="110m")
    ax.set_extent((-250, -80, -40, 70),crs=projPC)

    gl = ax.gridlines(crs=projPC, draw_labels=true,
                        linewidth=2, color="gray", alpha=0, linestyle="--")
    gl.xlabels_top = false
    gl.ylabels_right = false
    for ff in 1:5
        cf[1] = ax.pcolormesh(λ[ff], ϕ[ff],  β_mask[ff],
        vmin = b1, vmax = b2, shading="nearest", transform=projPC, 
        rasterized = true, cmap = colorway)               
    end
    ax.set_title(labels[i])
    ax.tick_params(bottom=true, left=true)
end


# ax.add_feature(ECCOonPoseidon.cartopy.feature.LAND)
fig.tight_layout()
fig
# ax.set_title(labels[i])


fig.suptitle("Temperature Trends between 2-3 km depth)")
# fig.savefig(plotsdir() * "/AGU_Plots/TempTrendClimatologyNoCbar_" * region * ".pdf",bbox_inches = "tight")
fig.tight_layout()
# cbar = fig.colorbar(cf[1], orientation="horizontal", fraction=0.04, extend = "both", label = L"^\circ" *"C per century")
# fig.savefig(plotsdir() * "/AGU_Plots/TempTrendClimatology_" * region * ".pdf",bbox_inches = "tight")
fig