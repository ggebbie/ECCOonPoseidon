#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
using .OHC_helper
import NaNMath as nm
import PyPlot as plt
include(srcdir("config_exp.jl"))
@pyimport seaborn as sns;
cmo = pyimport("cmocean.cm");

sns.set_theme(context = "notebook", style = "ticks",
              palette = sns.color_palette("colorblind"));
(ϕ,λ) = latlonC(γ)
tecco = 1992+1/24:1/12:2018
runpath,diagpath = listexperiments(exprootdir());

ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

#load in latitude mask 
ocean_mask = OHC_helper.wet_pts(Γ)
suffix = "2to3"
uplvl = -2e3; botlvl = -3e3
lvls = findall( botlvl .<= z[:].<= uplvl)

#create volume mask
area = readarea(γ)
cell_depths = OHC_helper.get_cell_depths(ocean_mask, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));

tecco = 1992+1/24:1/12:2018
nz = 50; tstart = 12*3; tstop = 312; nt = tstop - tstart + 1;
E,F = trend_matrices(tecco[tstart:tstop])
proj0 = ECCOonPoseidon.cartopy.crs.Robinson(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

sum_fluxes(d) = sum(d[varn] for varn in ["κxyθ", "uvθ", "wθ", "κzθ"])
fname = datadir("native/THETA_budget_spatial" * "_" * suffix * "_1995_2017.jld2")
dθ_iter129 = load(fname)["dθs"]["iter129_bulkformula"]
dθ_iter0 = load(fname)["dθs"]["iter0_bulkformula"]

dθ_plot = Dict()
dθ_plot["Iteration 129"] = dθ_iter129["dθ"] 
dθ_plot["Iteration 0"] = dθ_iter0["dθ"] 
dθ_plot["Iteration 129 minus Iteration 0"] = dθ_plot["Iteration 129"] .- dθ_plot["Iteration 0"]

ocn_reg = LLCcropC(ocean_mask,γ) 
reg_λ = LLCcropC(λ,γ); reg_ϕ = LLCcropC(ϕ,γ)
bounds = 0.3
levels = -bounds:0.05:bounds

fig, axes = plt.subplots(1, 3, figsize=(15,20), subplot_kw=Dict("projection"=> proj0))
CF = Any[]
for (i, key) in enumerate(["Iteration 129", "Iteration 0", "Iteration 129 minus Iteration 0"])
    ax = axes[i]
    data = LLCcropC(dθ_plot[key], γ)
    data[ocn_reg .== 0] .= NaN
    data .= data .* (100 * 3.154e+7)
    println(nm.maximum(data))
    cf = ax.pcolormesh(reg_λ, reg_ϕ,  data, transform=projPC, 
    cmap = cmo.balance, vmin = -bounds, vmax = bounds)   
    push!(CF, cf)
    gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
    gl.top_labels = false
    gl.bottom_labels = true

    gl.right_labels = false
    if i != 1
        gl.left_labels = false
    end
    ax.set_title(key)
end
[a.coastlines() for a in axes]
[a.set_extent((120, 285, -70, 56),crs=projPC) for a in axes]

fig.subplots_adjust(wspace=0.05)
fig.colorbar(CF[1], ax = axes, orientation = "horizontal", fraction = 0.012, 
label = L"[°C / century ]", pad=0.03, extend = "both")
fig.suptitle("Mid-Depth Temperature Trends in ECCO [z=2-3km, t = 1992-2017]", y=0.35)
fig
fig.savefig(plotsdir("native/HeatBudgetSpatialDifference_" * suffix * ".png"), bbox_inches = "tight", dpi = 500)
