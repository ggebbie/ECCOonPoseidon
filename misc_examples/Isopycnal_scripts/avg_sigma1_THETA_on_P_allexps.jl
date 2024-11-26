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

region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = OHC_helper.get_cell_volumes(area, cell_depths);

where_dry = findall(x -> isapprox(x, 0), cell_volumes)

sig1grid = sigma1grid()
nσ = length(sig1grid); nz = length(z)

tecco= 1992+1/24:1/12:2018 # ecco years
E,F = trend_matrices(tecco)
experiments = keys(shortnames)
plts = []
labels = ["Forcing + Init. Cond. ", "Forcing", "Init. Cond.", "Neither"]

σ1datadir(x="") = OHC_helper.σ1datadir(x)
!ispath(σ1datadir()) && mkpath(σ1datadir())

global expname

for (i, expr) in enumerate(experiments)
    expname = expr
    println(expname)
    #load in trends on sigma1 grid 
    fname = expname * region * "_AVG_P_sigma1.jld2" 
    !ispath(σ1datadir(fname)) && (include("get_avgP_on_sigma1.jl"))
    Pσ = load(σ1datadir(fname), "P")

    #load in average sigma1 vals
    fname = expname * region * "_AVG_THETA_sigma1.jld2"
    !ispath(σ1datadir(fname)) && (include("get_THETA_on_sigma1.jl"))
    θσ = load(σ1datadir(fname), "θ")

    sigma2p = mean(abs.(Pσ), dims = 2) #quick fix for P not always being negative
    p2z = gsw_z_from_p.(sigma2p, 30., 0., 0.)
    #get time anomalies 
    θσanom = θσ .- mean(θσ, dims = 2)
    
    #plotting
    wherenan = (!isnan).(vec(p2z))

    println("percentage of valid pts: ", 100 * sum(wherenan) / length(p2z))
    p = plot([NaN], [NaN], label = nothing, ylabel = "Depth [m]", 
    xlabel = "Year", xlims = extrema(tecco), title = labels[i])

    contourf!(p, tecco, abs.(p2z[wherenan]), 100 .* θσanom[wherenan, :], 
    color = :balance, levels = -2.5:0.1:2.5,   
    yflip = true, ylims = (500, 3500), lw=0.5, clabels=true, cbar=false)
    push!(plts, p)

end

p = plot(title = region * " Temperature Anomaly [cK] in ECCO",framestyle=nothing,showaxis=false,xticks=false,yticks=false,margin=0mm)
plts2 = cat(p, plts, dims = 1)

l = @layout [a{0.01h}; grid(2,2)]
p2 = plot(plts2..., layout=l, size = (900, 650))
savefig(p2, plotsdir(region * "_AllExps_TempAnom_P2Z.png"))
p2