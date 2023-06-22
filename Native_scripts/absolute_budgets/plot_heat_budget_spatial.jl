#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
using .OHC_helper
import NaNMath as nm
import PyPlot as plt
include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");

(ϕ,λ) = latlonC(γ)
tecco = 1992+1/24:1/12:2018
runpath,diagpath = listexperiments(exprootdir());

ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

#load in latitude mask 
ocean_mask = OHC_helper.wet_pts(Γ)
suffix = "2to3"
uplvl = -2e3; botlvl = -3e3
_, lvls = findmin(abs.(z .- (2500)))

#create volume mask
area = readarea(γ)
cell_depths = OHC_helper.get_cell_depths(ocean_mask, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));

tecco = 1992+1/24:1/12:2018
nz = 50; tstart = 12*3; tstop = 312; nt = tstop - tstart + 1;
E,F = trend_matrices(tecco[tstart:tstop])

sum_fluxes(d) = sum(d[varn] for varn in ["κxyθ", "uvθ", "wθ", "κzθ"])
 
expname = "iter0_bulkformula"
# fname = datadir("THETA_budget_spatial" * "_" * expname * "_" * suffix *".jld2")
fname = datadir("native/THETA_budget_spatial" * "_" * suffix *"_1995_2017.jld2")
dθ = load(fname)["dθs"][expname]

dθ_plot = Dict()
dθ_plot["θ Trend"] = dθ["dθ"] 
dθ_plot["∇(⟨uθ⟩)"] = dθ["uvθ"] .+ dθ["wθ"]
dθ_plot["∇²(κ⟨θ⟩)"] = dθ["κxyθ"] .+ dθ["κzθ"]
dθ_plot["d⟨θ⟩/dt"] = sum_fluxes(dθ)
dθ_plot["Residual"] = dθ_plot["θ Trend"] .- dθ_plot["d⟨θ⟩/dt"]

ocn_reg = LLCcropC(ocean_mask,γ) 
reg_λ = LLCcropC(λ,γ); reg_ϕ = LLCcropC(ϕ,γ)
bounds = 0.3
levels = -bounds:0.05:bounds

proj0 = ECCOonPoseidon.cartopy.crs.Robinson(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

fig, axes = plt.subplots(1, 3, figsize=(30,10), subplot_kw=Dict("projection"=> proj0))
CF = Any[]
for (i, key) in enumerate(["d⟨θ⟩/dt", "∇(⟨uθ⟩)", "∇²(κ⟨θ⟩)"])
    ax = axes[i]
    data = LLCcropC(dθ_plot[key], γ)
    data[ocn_reg .== 0] .= NaN
    data .= data .* (100 * 3.154e+7)
    println(nm.maximum(data))
    cf = ax.pcolormesh(reg_λ, reg_ϕ,  data, transform=projPC, 
    cmap = cmo.balance, vmin = -bounds, vmax = bounds)   
    push!(CF, cf)
    ax.coastlines(resolution="110m")
    ax.set_extent((-180, 180, -70, 56),crs=projPC)
    gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
    gl.top_labels = false
    gl.bottom_labels = false

    gl.right_labels = false
    if i != 1
        gl.left_labels = false
    end
    ax.set_title(key)
end
fig.subplots_adjust(wspace=0.05)
fig.colorbar(CF[1], ax = axes, orientation = "horizontal", fraction = 0.04, 
label = L"[°C / century ]", pad=0.1, extend = "both")
fig.savefig(plotsdir("native/HeatBudgetSpatial_" * suffix * "_" *expname * ".png"), bbox_inches = "tight")
fig

fig, axes = plt.subplots(1, 3, figsize=(30,10), subplot_kw=Dict("projection"=> proj0))
CF = Any[]
for (i, key) in enumerate(["θ Trend", "d⟨θ⟩/dt", "Residual"])
    ax = axes[i]
    data = LLCcropC(dθ_plot[key], γ)
    data[ocn_reg .== 0] .= NaN
    data .= data .* (100 * 3.154e+7)
    println(nm.maximum(data))
    cf = ax.pcolormesh(reg_λ, reg_ϕ,  data, transform=projPC, 
    cmap = cmo.balance, vmin = -bounds, vmax = bounds)   
    push!(CF, cf)
    ax.coastlines(resolution="110m")
    ax.set_extent((-180, 180, -70, 56),crs=projPC)
    gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
    gl.top_labels = false
    gl.bottom_labels = false

    gl.right_labels = false
    if i != 1
        gl.left_labels = false
    end
    ax.set_title(key)
end
fig.subplots_adjust(wspace=0.05)
fig.colorbar(CF[1], ax = axes, orientation = "horizontal", fraction = 0.04, 
label = L"[°C / century ]", pad=0.1, extend = "both")
fig.savefig(plotsdir("native/HeatBudgetSpatialDifference_" * suffix * "_" *expname * ".png"), bbox_inches = "tight")
fig