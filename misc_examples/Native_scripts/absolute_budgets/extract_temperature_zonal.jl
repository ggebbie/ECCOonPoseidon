#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall, .OHC_helper
import NaNMath as nm
import PyPlot as plt
include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");

#using . is faster 
(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());

# abbreviations for each experiment for labels, etc.
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
marks = expsymbols()
nexp = length(shortnames) # number of experiments
 
ocean_mask = OHC_helper.wet_pts(Γ)
region = "PAC56"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));

tecco = 1992+1/24:1/12:2018
theta_adjust = OHC_helper.get_theta_init(γ)

X = OHC_helper.ma_zonal_sum(ϕ .* area) ./ OHC_helper.ma_zonal_sum(area)
Y = -z[:]
data = OHC_helper.ma_zonal_avg(theta_adjust, cell_volumes)
vmax = 7.5; levels = collect(-vmax:2.5:vmax)

fig, ax = plt.subplots(figsize = (10, 5))
cf = ax.contourf(X, Y,  100 .* data, cmap = cmo.balance, vmin = -vmax, vmax = vmax, levels = levels, extend = "both")
cs = ax.contour(X, Y,  100 .* data, colors = "black", levels = levels)
fig.colorbar(cf, label = "cK")
ax.set_title("ECCO Control Adjustment in Pacific Ocean")
ax.invert_yaxis(); 
ax.set_xlabel("Latitude")
ax.set_ylabel("Depth [m]")
ax.set_xticks(-40:10:60)
ax.set_xlim(-39, 59)
fig


uplvl = -2e3; botlvl = -3e3
lvls = findall( botlvl .<= z[:].<= uplvl)
cell_depths = OHC_helper.get_cell_depths(ocean_mask, ΔzF, Γ.hFacC); 
H = OHC_helper.sum_vertical(cell_depths[:, lvls], γ); [H.f[ijk][iszero.(H.f[ijk])] .= Inf for ijk in eachindex(H)]

data = OHC_helper.depth_average(theta_adjust, cell_depths, H, γ)
data = LLCcropC360(data, γ; modify_λ = false)
data[iszero.(data)] .= NaN
proj0 = ECCOonPoseidon.cartopy.crs.Robinson(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

fig, axes = plt.subplots(figsize=(10,10), subplot_kw=Dict("projection"=> proj0))
λ_crop = LLCcropC360(λ, γ; modify_λ = true)
ϕ_crop = LLCcropC360(ϕ, γ; modify_λ = false)
axes.contourf(λ_crop, ϕ_crop, 100 .* data, cmap = cmo.balance, vmin = -vmax, vmax = vmax, levels = levels, extend = "both", transform=projPC)
# cs = ax.contour(λ_crop, ϕ_crop, 100 .* data, colors = "black", levels = levels)
axes.coastlines("110m")
fig
axes.set_title("ECCO Control Adjustment in Pacific Ocean")
fig