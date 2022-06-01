#  First steps: 1. go into Reemergence project directory. 2. go into julia REPL package mode with `]`. 3. `activate .` 4. Backspace to return to command mode in REPL.
# Modeled after "filter_interannual.jl"
# Solve for linear trends at all points.

include("../src/intro.jl")
includet("../src/OHC_helper.jl")

using Revise 
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools,
    PyPlot, PyCall, JLD2, DrWatson, Statistics, JLD2
using Main.OHC_helper
include(srcdir("config_exp.jl"))
ccrs = pyimport("cartopy.crs")
cm = pyimport("cmocean.cm")
cfeature = pyimport("cartopy.feature")

colorway = cm.balance

do_constant_density=true 
do_total=true 
output2file = true
ocean_mask = Γ.hFacC[:,1] 
ocean_mask[findall(ocean_mask.>0.0)].=1.0
ocean_mask[findall(ocean_mask.<=0.0)].=NaN

basin_name="Pacific"
basinID=findall(basin_list.==basin_name)[1]

for ff in 1:length(area)
    # above_SO = (ϕ[ff] .> -56.0) #removes southern ocean 
    ocean_mask[ff] .= ocean_mask[ff].*(basins[ff].==basinID)
end

cell_depths = get_cell_depths(ocean_mask, Δz, Γ.hFacC)

(ϕ,λ) = latlonC(γ) #each tile has dims=(λ, ϕ)
area = readarea(γ)
for ff ∈ [3,4] 
    λ[ff][λ[ff] .<= 0 ] = λ[ff][λ[ff] .<= 0 ] .+ 360
end
runpath,diagpath = listexperiments(exprootdir())
shortnames = expnames()
marks = expsymbols()
nexp = length(shortnames) # number of experiments
fileroot = "state_3d_set1"
filedir = "ECCO_vars/"
filename = "THETAs"

# assume monthly averages, ECCOv4r4
tstart = 1992 + 1/24
tend = 2018
tecco = range(tstart,step=1/12,stop=2018)
nt = length(tecco)

# get weight matrix for finding trends
E,F = trend_matrices(tecco)
for (keys,values) in shortnames
    # pre-allocate β, linear trends
    nz = length(z)

    expname = keys
    @load datadir(filedir*filename*"_"*expname*".jld2") var_exp
    β = compute_β(var_exp, F, γ); var_exp = nothing 
    # lvl1 = -2300; lvl2 = -3300
    # lvl1 = -4000; lvl2 = -5000
    # lvl_idx = findall( lvl2 .<= z[:].<= lvl1)
    for lvl_idx in 1:nz
        β_sel = β[:, lvl_idx]
        var_sel = level_weighted_mean(β_sel, γ, cell_depths) .* ocean_mask  #trend over 26 years
        var_sel = var_sel .* 26
        fig, ax = plt.subplots(1, 1, figsize=(5,5), 
                            subplot_kw=Dict("projection"=> ccrs.PlateCarree(central_longitude = -150)))
        ax.set_extent((110, -70, -55, 60))
        # ax.add_feature(cfeature.COASTLINE)

        ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=true,
                            linewidth=2, color="gray", alpha=0, linestyle="--")
        cf = Vector{Any}(undef ,1)

        min_val = min_mesharray(var_sel)
        # max_val = max_mesharray(var_sel)
        max_val = -min_val
        for ff in 1:length(var_sel)
                cf[1] = ax.pcolormesh(λ[ff], ϕ[ff], var_sel[ff],shading="nearest", 
                transform=ccrs.PlateCarree(), rasterized = true, vmin = min_val, vmax = max_val,
                cmap = colorway)
        end
        z_string = string(abs(round(z[lvl_idx], digits = 2)))
        ax.set_title("∂θ/ (1992-2018) at " * z_string * "m \n " * values)
        cbar = fig.colorbar(cf[1])
        fig.tight_layout()

        #Saving the Figure
        folder = "/home/ameza/ECCOonPoseidon/plots/OHC Climatologies/" * expname * "/"
        mkpath(folder[1:end-1])
        # lvls_str = string(abs(lvl1)) * "to" * string(abs(lvl2))
        lvls_str = "lvl" * string(lvl_idx)
        outputfile = plotsdir(folder * "ΔΘ" *"climatology_" * expname * lvls_str * ".pdf")

        println("saving plots at.. " * outputfile)
        fig.savefig(outputfile)
        close("all")
    end
end