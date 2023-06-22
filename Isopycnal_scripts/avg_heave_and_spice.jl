#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings, GibbsSeaWater, PyCall
using Interpolations
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))
include("./IsopycnalHelpers.jl")

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

@pyimport cmocean.cm as cmo
@pyimport seaborn as sns;

sns.set_theme(context = "poster", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));

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
#computed on standard pressures, 
#may have a different value if using pressures diagnosed from ECCO
#https://ecco-v4-python-tutorial.readthedocs.io/VectorCalculus_ECCO_barotropicVorticity.html

tecco= 1992+1/24:1/12:2018 # ecco years
E,F = trend_matrices(tecco)
experiments = keys(shortnames)
labels = ["Forcing + Init. Cond. ", "Forcing", "Init. Cond.", "Neither"]
titles = ["Net", "Heave Signal", "Water Mass Signal"]

!ispath(σ1datadir()) && mkpath(σ1datadir())
for (i, expt) in enumerate(keys(shortnames))
    expname = expt
    #load in trends on sigma1 grid 
    fname = expname * region * "_AVG_P_sigma1.jld2" 
    Pσ = load(σ1datadir(fname), "P")

    #load in average sigma1 vals
    fname = expname * region * "_AVG_THETA_sigma1.jld2"
    θσ = load(σ1datadir(fname), "θ")

    #load in average sigma1 vals
    fname = expname * region * "_THETA_budget_z.jld2"
    θz = load(datadir(fname), "dθ")["θ"]


    sigma2p = mean(Pσ, dims = 2)
    p2z = gsw_z_from_p.(abs.(sigma2p), 30., 0., 0.)

    θσ = θσ[(!isnan).(p2z[:, 1]), :]
    p2z = p2z[(!isnan).(p2z[:, 1])]


    iz = findall(-3500 .< z .< -400); nt = 312

    θσ_itp = zeros(length(iz), nt)
    #take βσ from sigma to depth 
    for tt = 1:nt
        itp = interpolate((-p2z,), θσ[:, tt], Gridded(Linear()))
        θσ_itp[:, tt] = itp(-z[iz][:])
    end

    #get time anomalies 
    θzanom = θz[iz, :] .- mean( θz[iz, :], dims = 2)
    θσanom = θσ_itp .- mean(θσ_itp, dims = 2)
    θhvanom = θzanom .- θσanom

    vars =  [θzanom, θhvanom, θσanom]
    bounds = 2.0
    levels = -bounds:0.5:bounds

    fig,axs=plt.subplots(1,3, sharey = true, figsize = (30, 7.5))
    cfs = []
    for (i, var) in enumerate(vars)
        cf = axs[i].contourf(collect(tecco), -z[iz][:], 100 .* var, 
        cmap = cmo.balance, vmin = -bounds, vmax = bounds, levels = levels, extend = "both")
        axs[i].set_title(titles[i])
        axs[i].set_xlabel("Year")
        push!(cfs, cf)
        cs = axs[i].contour(collect(tecco), -z[iz][:], 100 .* var, 
        colors="k", vmin = -bounds, vmax = bounds, levels = levels, linewidths = 0.75)
        axs[i].clabel(cs, fontsize=15, inline=true, fmt = "%.1f")
    end
    axs[1].set_ylabel("Depth [m]")
    axs[1].set_ylim(500, 3200)
    axs[1].invert_yaxis()
    fig.suptitle(region * " Temperature Anomaly [cK] in ECCO")
    fig.savefig(plotsdir(region * "_" * expname * "_HeaveandSpice.png"))
end 