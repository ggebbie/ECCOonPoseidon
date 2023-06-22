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

cmo = pyimport("cmocean.cm");

(ϕ,λ) = latlonC(γ)
tecco = 1992+1/24:1/12:2018
runpath,diagpath = listexperiments(exprootdir());

ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

#load in latitude mask 
region = "PAC"; 
ocean_mask = OHC_helper.wet_pts(Γ);
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not");
area = readarea(γ);
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));

tecco = 1992+1/24:1/12:2018;
nz = 50; tstart = 12*3; tstop = 312; nt = tstop - tstart + 1;

ocn_reg = LLCcropC(ocean_mask,γ);
reg_λ = LLCcropC(λ,γ); reg_ϕ = LLCcropC(ϕ,γ);

cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));
GTF_zonal = OHC_helper.ma_zonal_avg(GTF3d, cell_volumes)

X = OHC_helper.ma_zonal_sum(ϕ .* area) ./ OHC_helper.ma_zonal_sum(area)
Y = -z[:]

levels = 0:0.01:0.3

fig, ax = plt.subplots(figsize=(10,5), sharey = true)
data = (100 * 3.154e+7) .* deepcopy(GTF_zonal) #C/s -> C/ century 

cf = ax.contourf(X, Y,  data, cmap = cmo.thermal, extend = "both", 
vmin = 0, vmax = 0.3, levels = levels)   
ax.set_xticks(-40:10:60)
ax.set_xlim(-39, 59)
fig.colorbar(cf, ax = ax, orientation = "horizontal", fraction = 0.07, 
label = L"[°C / century ]",pad=0.2)
ax.set_ylim(1000, 5000); ax.invert_yaxis(); 
ax.set_ylabel("Depth [m]")
fig

sum(GTF3d .* cell_volumes) / sum(cell_volumes)


fig, ax = plt.subplots(figsize=(10,5), subplot_kw=Dict("projection"=> proj0));
H = OHC_helper.sum_vertical(cell_depths[:, 38:42], γ); 
da = (100 * 3.154e+7) .* OHC_helper.depth_average(GTF3d2[:, 38:42], cell_depths[:, 38:42], H, γ)
data = LLCcropC(da, γ);data[data .== 0.0] .= NaN;
bounds = nm.maximum(data);
ax.add_feature(ECCOonPoseidon.cartopy.feature.LAND)

cf = ax.pcolormesh(reg_λ, reg_ϕ,  data, transform=projPC, cmap = cmo.balance, 
vmin = -0.3, vmax = 0.3);
ax.coastlines(resolution="110m", facecolor = ECCOonPoseidon.cartopy.feature.COLORS["land"]);
ax.set_extent((-180, 180, -70, 56),crs=projPC);
ax.set_title("Geothermal Heating");
fig.colorbar(cf); 
fig



# fname = "/batou/eccodrive/files/Version4/Release4/input_init/geothermalFlux.bin";
# GTF = read_bin(fname,Float32,γ)

# fig, ax = plt.subplots(figsize=(10,5), subplot_kw=Dict("projection"=> proj0));
# data = LLCcropC(GTF .* PAC_msk, γ);data[data .== 0.0] .= NaN;
# bounds = nm.maximum(data);

# cf = ax.contourf(reg_λ, reg_ϕ,  data, transform=projPC, cmap = cmo.thermal);
# ax.coastlines(resolution="110m");
# ax.set_extent((-180, 180, -70, 56),crs=projPC);
# ax.set_title("Geothermal Heating");
# fig.colorbar(cf); 
# fig


# bathy_mask = deepcopy(Γ.hFacC)
# for a in eachindex(bathy_mask)
#     bathy_mask.f[a][bathy_mask.f[a] .> 0.0] .= 1.0
# end

# cell_depths = OHC_helper.get_cell_depths(ocean_mask, ΔzF, Γ.hFacC); 

# for ff in 1:5, k = 1:49
#     bathy_mask.f[ff, k] .= bathy_mask.f[ff, k] .- bathy_mask.f[ff, k+1]
# end

# rho = 1029; cp = 3994; 

# GTF3d = MeshArray(γ, Float32, 50)
# for a in eachindex(GTF3d)
#     inv_ρcpdz = cell_depths.f[a] .* (rho * cp); 
#     inv_ρcpdz[inv_ρcpdz .== 0] .= Inf
#     GTF3d.f[a] .= (GTF.f[a[1]] .* bathy_mask.f[a]) ./ (inv_ρcpdz)
# end

# GTF3d2 = OHC_helper.get_geothermalheating(Γ, γ)