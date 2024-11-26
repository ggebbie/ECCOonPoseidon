#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2, 
Printf, DataFrames, LaTeXStrings,
Plots, GibbsSeaWater, PyCall,
Plots.PlotMeasures, Interpolations

import NaNMath as nm

@pyimport seaborn as sns;
@pyimport pandas as pd;

colors =  sns.color_palette("colorblind")[1:4]

include(srcdir("config_exp.jl"))
include("./IsopycnalHelpers.jl")

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

region = "NPAC30"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = OHC_helper.get_cell_volumes(area, cell_depths);

where_dry = findall(x -> isapprox(x, 0), cell_volumes)

sig1grid = sigma1grid()
nσ = length(sig1grid); nz = length(z)
#computed on standard pressures, 
#may have a different value if using pressures diagnosed from ECCO
#https://ecco-v4-python-tutorial.readthedocs.io/VectorCalculus_ECCO_barotropicVorticity.html

tecco= 1992+1/24:1/12:2018 # ecco years
E,F = trend_matrices(tecco)
experiments = keys(shortnames)
plts = []
labels = ["Forcing + Init. Cond. & κ", "Forcing", "Init. Cond. & κ", "None"]
diff(x) = x[2:end] - x[1:end-1]
derivative(y, x) = (diff(y) ./ diff(x), x[1:end-1] .+ (diff(x) ./ 2) )

σ1datadir(x="") = OHC_helper.σ1datadir(x)
!ispath(σ1datadir()) && mkpath(σ1datadir())
p1 = plot([NaN], [NaN], label = nothing)
p2 = plot([NaN], [NaN], label = nothing)

for (i, expr) in enumerate(experiments)

    expname =expr

    #load in Pressure  on sigma1 ßtimeseriesß
    fname = expname * region * "_AVG_P_sigma1.jld2" 
    !ispath(σ1datadir(fname)) && (include("get_avgP_on_sigma1.jl"))
    Pσ = load(σ1datadir(fname), "P")

    #load in temp. on sigma1 timeseries
    fname = expname * region * "_AVG_THETA_sigma1.jld2"
    !ispath(σ1datadir(fname)) && (include("get_THETA_on_sigma1.jl"))
    θσ = load(σ1datadir(fname), "θ")
    sigma2p = mean(Pσ, dims = 2)
    p2z = gsw_z_from_p.(abs.(sigma2p), 30., 0., 0.)

    # find maximum valid depth 
    _, idx = findmax(!isnan, p2z[:])
    p2z = p2z[idx:end]
    θσ = θσ[idx:end, :]
    θ̄= mean(θσ, dims = 2)[:]

    #take 1st derivative 
    dzθ, z_deriv = derivative(θ̄, p2z)
    # dz = p2z[2:end] .- p2z[1:end-1]
    # dzθ = (θ̄[2:end] .- θ̄[1:end-1]) ./ dz[:]
    # z_deriv = p2z[1:end-1] .+ (dz[1:end] ./2)

    #take 2nd derivative 
    dzθ, z_deriv = derivative(dzθ, z_deriv)

    # dz = z_deriv[2:end] .- z_deriv[1:end-1]
    # dzθ = (dzθ[2:end] .- dzθ[1:end-1]) ./ dz
    # z_deriv = z_deriv[1:end-1] .+ (dz[1:end] ./2)


    plot!(p1, dzθ[20:end], z_deriv[20:end], 
    title = L"\frac{d^2\bar{\theta}}{dz^2}", c = RGB(colors[i]...),
    label = nothing,
    xrotation = 20, ylim = (-3500, -500), xlim = (0, 1.2e-5))

    dtθ̄= (θσ[:, 2:end] .- θσ[:, 1:end-1]) ./ 2.628e+6
    dtθ̄= mean(dtθ̄, dims = 2)
    plot!(p2, dtθ̄[20:end], p2z[20:end], c = RGB(colors[i]...),
    title = L"\frac{d\bar{\theta}}{dt}", label = labels[i],
    xrotation = 20, ylim = (-3500, -500), xlim = (0, 1.2e-10),
    legend=:bottomright)
end

pT = plot(title = region * " Temperature Trend Analysis ",
framestyle=nothing,showaxis=false,xticks=false,yticks=false,margin=0mm)
plts2 = [pT, p1, p2]
l = @layout [a{0.01h}; grid(1,2)]
p2 = plot(plts2..., layout=l, size = (700, 500))
# savefig(p2, plotsdir(region * "_TempTrendAnalysis_P2Z.png"))
p2