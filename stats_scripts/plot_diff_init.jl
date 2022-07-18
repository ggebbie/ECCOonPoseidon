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
cm = pyimport("cmocean.cm")

colorway = cm.balance

do_constant_density=true 
do_total=true 
output2file = true
ocean_mask = Γ.hFacC[:,1] 
ocean_mask[findall(ocean_mask.>0.0)].=1.0
ocean_mask[findall(ocean_mask.<=0.0)].=NaN

basin_name="Pacific"
basinID=findall(basin_list.==basin_name)[1]


cell_depths = get_cell_depths(ocean_mask, Δz, Γ.hFacC)
for ff in 1:length(ocean_mask)
    # above_SO = (ϕ[ff] .> -56.0) #removes southern ocean 
    ocean_mask[ff] .= ocean_mask[ff].*(basins[ff].==basinID)
end

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
uplvl = -1700; botlvl = -2300
# lvl1 = -4000; lvl2 = -5000
lvl_idx = findall( botlvl .<= z[:].<= uplvl)

#load iter0 
expname = "iter0_bulkformula"
expshort = "0bf"
@load datadir(filedir*filename*"_"*expname*".jld2") var_exp
iter0_init = var_exp[1][:, lvl_idx]
var_exp = nothing
for (keys,values) in shortnames
    # pre-allocate β, linear trends
    expname = keys
    @load datadir(filedir*filename*"_"*expname*".jld2") var_exp
    exp_init = var_exp[1][:, lvl_idx]
    var_exp = nothing
    diff = similar(exp_init)
    for lvl in 1:length(lvl_idx)
        diff[:, lvl] = exp_init[:, lvl] .- iter0_init[:, lvl]
    end
    var_sel = level_weighted_mean(diff, γ, cell_depths)  

    proj = ECCOonPoseidon.cartopy.crs.Robinson()
    proj0 = ECCOonPoseidon.cartopy.crs.Robinson(central_longitude=-150)
    projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
    fig, ax = plt.subplots(1, 1, figsize=(6,5), 
                           subplot_kw=Dict("projection"=> proj0))
    ax.set_global()
    ax.coastlines()
    ax.set_extent((110, -70, -55, 60))
    lvls_str = string(abs(uplvl)) * "to" * string(abs(botlvl))
    ax.set_title("θ trend difference \n at " * lvls_str * " m \n " * values * "-"* "0bf")

    min_val = minimum(MeshArrays.mask(var_sel  .* ocean_mask ,Inf))
    max_val = -min_val
    tight_layout()
    println(values)
    println(min_val)
    cf = Vector{Any}(undef ,1)
    for ff in 1:length(var_sel)
            cf[1] = ax.pcolormesh(λ[ff], ϕ[ff], var_sel[ff],shading="nearest", 
            transform=projPC, rasterized = true, vmin = min_val, vmax = max_val,
            cmap = colorway)
    end
    # depthlbl = string(abs(round(z[lvl_idx], digits = 2)))
    cbar = fig.colorbar(cf[1], label = L"^\circ"*"C")

    #Saving the Figure
    folder = "/home/ameza/ECCOonPoseidon/plots/OHC Climatologies/" * expname * "/"
    mkpath(folder[1:end-1])
    # lvls_str = "lvl" * string(lvl_idx)
    outputfile = plotsdir(folder * "Θ" *"climatology_diffinit_" * expname * lvls_str * ".pdf")

    savefig(outputfile)
    println("saving plots at.. " * outputfile)

    close("all")
end
