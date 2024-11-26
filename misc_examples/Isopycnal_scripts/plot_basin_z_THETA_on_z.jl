#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings,
Plots, GibbsSeaWater

using .OHC_helper
using Interpolations
import NaNMath as nm
using Plots.PlotMeasures

include(srcdir("config_exp.jl"))
include("./IsopycnalHelpers.jl")

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

region = "PAC56"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = OHC_helper.get_cell_volumes(area, cell_depths);

where_dry = findall(x -> isapprox(x, 0), cell_volumes)

nz = length(z)
#computed on standard pressures, 
#may have a different value if using pressures diagnosed from ECCO
#https://ecco-v4-python-tutorial.readthedocs.io/VectorCalculus_ECCO_barotropicVorticity.html

tecco= 1992+1/24:1/12:2018 # ecco years
E,F = trend_matrices(tecco)
experiments = keys(shortnames)
plts = []
labels = ["Forcing + Init. Cond. ", "Forcing", "Init. Cond.", "Neither"]
σ1datadir(x="") = OHC_helper.σ1datadir(x)

!ispath(σ1datadir()) && mkpath(σ1datadir())

expname = "iter129_bulkformula"
print(expname)
#load in trends on sigma1 grid 
fname = expname * region * "_AVG_THETA_z.jld2"

!ispath(σ1datadir(fname)) && (include("get_THETA_on_z.jl"))

θz = load(σ1datadir(fname), "θ")
θzanom = θz .- mean(θz, dims = 2)

#plotting
p = plot([NaN], [NaN], label = nothing, ylabel = "Depth [m]", 
xlabel = "Year", xlims = extrema(tecco), title = labels[i])

contourf!(p, tecco, abs.(z), 100 .* θzanom, 
color = :balance, levels = -2.5:0.1:2.5,   
yflip = true, ylims = (500, 5000), lw=0.5, clabels=true, cbar=false)
push!(plts, p)


p = plot(title = region * " Temperature Anomaly [cK] in ECCO",
framestyle=nothing,showaxis=false,xticks=false,yticks=false,margin=0mm)

l = @layout [a{0.01h}; grid(1,1)]
p2 = plot(vcat(p, plts)..., layout=l,size = (600, 500))
savefig(p2, plotsdir(region * "_" * expname * "_TempAnom_Z.png"))
p2

