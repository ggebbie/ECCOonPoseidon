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
using Interpolations
import NaNMath as nm
using Plots.PlotMeasures

include(srcdir("config_exp.jl"))
include("./IsopycnalHelpers.jl")

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = OHC_helper.get_cell_volumes(area, cell_depths);

where_dry = findall(x -> isapprox(x, 0), cell_volumes)

sig1grid = sigma1grid()
nσ = length(sig1grid)
nz = size(cell_volumes, 2)
#computed on standard pressures, 
#may have a different value if using pressures diagnosed from ECCO
#https://ecco-v4-python-tutorial.readthedocs.io/VectorCalculus_ECCO_barotropicVorticity.html
runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

tecco= 1992+1/24:1/12:2018 # ecco years
E,F = trend_matrices(tecco)

#read in the first time step of S and θ
expname = "nosfcadjust"

#load in theta on sigma1 timeseries 
θσ = load(OHC_helper.σ1datadir(expname * region * "_AVG_THETA_sigma1.jld2"), "θ")
#load in average sigma1 vals
σ_avg = γ.read(OHC_helper.σ1datadir(expname * "_AVG_THETA_sigma1.data"), MeshArray(γ,Float32,nz))
σ_avg[where_dry] = 0.0

#get temperature trend by depth 
σz_profile = volume_average_by_depth(σ_avg, cell_volumes, γ)
get_min_ind(x) = findmin(x)[2]
σstozs = [get_min_ind(abs.(σ .- sig1grid)) for σ in σz_profile]
(min_id, max_id) = extrema(σstozs) .+ (1, -1)

#take βσ from sigma to depth 
itp = interpolate((σz_profile[1:end-2],), z[1:end-2], Gridded(Linear()))
mapsigma2z = itp(sig1grid[min_id:max_id])

#get time anomalies 
θσanom = θσ .- mean(θσ, dims = 2)
#plotting
p2 = contourf(  tecco, abs.(mapsigma2z), 100 .* θσanom[min_id:max_id, :], 
color = :balance, levels = -2.:0.1:2.,   
yflip = true, ylims = (500, 3500), lw=0.5, clabels=true, cbar=false,
)

plot!(p2, [NaN], [NaN], label = nothing, ylabel = "Depth [m]", 
xlabel = "Year", xlims = extrema(tecco), title = region * " Temperature Anomaly [cK] in ECCO")
p2