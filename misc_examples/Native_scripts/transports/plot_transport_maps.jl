#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
using .OHC_helper
import NaNMath as nm
import PyPlot as plt
using ImageFiltering
include(srcdir("config_exp.jl"))
              
cm = pyimport("cmocean.cm");
sns = pyimport("seaborn");
sns.set_theme(context = "paper", style = "ticks",
              palette = sns.color_palette("colorblind"));

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

ocean_mask = OHC_helper.wet_pts(Γ)
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
cell_depths = OHC_helper.get_cell_depths(ocean_mask, ΔzF, Γ.hFacC); 
H = OHC_helper.smush(cell_depths, γ); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area .* PAC_msk, cell_depths));

runpath,diagpath = listexperiments(exprootdir());

#define control volume 
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= z[:].<= uplvl)
cell_volumes = cell_volumes .* PAC_msk
ctrl_vol = sum(cell_volumes[:, lvls])

meanvarname(x) = x * "_" * region * "_" * suffix *".data"
vars = ["iter129_bulkformula", "iter0_bulkformula"]
area32 = Float32.(area)
W_dict = Dict(); θ_dict = Dict(); Wθ_dict = Dict()
E_dict = Dict(); N_dict = Dict()
vol_avg(x) = sum(x .* cell_volumes[:, lvls[1]]) / sum(cell_volumes[:, lvls[1]])
inv_cell_volumes = deepcopy(cell_volumes); 
for a in eachindex(inv_cell_volumes)
    inv_cell_volumes.f[a][cell_volumes.f[a] .== 0.0] .= Inf; 
    inv_cell_volumes.f[a] = 1 ./ inv_cell_volumes.f[a]
end
cs, sn = OHC_helper.get_cs_and_sn(γ)

for expname in vars
    Ubar = γ.read(datadir(expname * "/" * meanvarname("U_EB_mean")),MeshArray(γ,Float32,50))
    Vbar = γ.read(datadir(expname * "/" * meanvarname("V_EB_mean")),MeshArray(γ,Float32,50))
    Wbar = γ.read(datadir(expname * "/" * meanvarname("W_EB_mean")),MeshArray(γ,Float32,50))
    sθbar = γ.read(datadir(expname * "/" * meanvarname("THETA_ref_mean")),MeshArray(γ,Float32,50))
    Ubar, Vbar = UVtoTrsp(Ubar, Vbar, Γ)
    Ebar, Nbar = rotate_UV_native(Ubar, Vbar, cs, sn)
    N_dict[expname] = OHC_helper.sum_vertical(Nbar[:, lvls], γ)
    E_dict[expname] = OHC_helper.sum_vertical(Ebar[:, lvls], γ)
    OHC_helper.interpolate_to_vertical_faces!(sθbar, sθbar)

    W_dict[expname] = Wbar[:, lvls[1]].* area32;
    θ_dict[expname] = sθbar[:, lvls[1]] .* PAC_msk; 
    Wθ_dict[expname] = W_dict[expname] .* θ_dict[expname]; 
    Wθ_dict[expname] = Wθ_dict[expname] .* inv_cell_volumes[:, lvls[1]]
end

diff_dict(d, x, y) = d[x] .- d[y]
diff_label = shortnames[vars[1]] * " minus " * shortnames[vars[2]]
E_dict[diff_label] = diff_dict(E_dict, vars[1], vars[2])
N_dict[diff_label] = diff_dict(N_dict, vars[1], vars[2])
Wθ_dict[diff_label] = diff_dict(Wθ_dict, vars[1], vars[2])
W_dict[diff_label] = diff_dict(W_dict, vars[1], vars[2])
θ_dict[diff_label] = diff_dict(θ_dict, vars[1], vars[2])
push!(vars, diff_label) 

proj0 = ECCOonPoseidon.cartopy.crs.Robinson(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
λ_wrap = OHC_helper.wrap_λ(λ) 
CF = Any[1, 2, 1]

bounds = [1, 1,1]
#plot τx
fig, ax = plt.subplots(1, 3, figsize=(15,20), 
                       subplot_kw=Dict("projection"=> proj0), sharey = true)
[a.coastlines() for a in ax]
[a.set_extent((120, 250, -56, 56),crs=projPC) for a in ax]
filterfunc(x) = imfilter(x, Kernel.gaussian(1))
for (i, var) in enumerate(vars)
    CF[i] = OHC_helper.pcolormesh_ma_filt(ax[i], λ_wrap, ϕ, -3.154e+9 .* Wθ_dict[var], 
                        bounds[i], cm.balance, projPC, true, filterfunc)
    ax[i].set_title(var)
end
# ax.set_extent((-180, 180, -70, 56),crs=projPC)
fig.colorbar(CF[1], ax = ax[1:2], fraction = 0.01, pad = 0.03, orientation = "horizontal", label = "[Sv]")
fig.colorbar(CF[3], ax = ax[3], fraction = 0.01,  pad = 0.03, orientation = "horizontal", label = "[Sv]")
fig

maximum(θ_dict["iter129_bulkformula"])

fig, ax = plt.subplots(1, 3, figsize=(15,20), 
                       subplot_kw=Dict("projection"=> proj0), sharey = true)
bounds = [1, 1,0.5] ./2

[a.coastlines() for a in ax]
[a.set_extent((120, 250, 0, 56),crs=projPC) for a in ax]
filterfunc(x) = imfilter(x, Kernel.gaussian(1))
for (i, var) in enumerate(vars)
    CF[i] = OHC_helper.pcolormesh_ma_filt(ax[i], λ_wrap, ϕ, -θ_dict[var], 
                        bounds[i], cm.balance, projPC, true, filterfunc)
    ax[i].set_title(var)
end
# ax.set_extent((-180, 180, -70, 56),crs=projPC)
fig.colorbar(CF[1], ax = ax[1:2], fraction = 0.01, pad = 0.03, orientation = "horizontal", label = "[Sv]")
fig.colorbar(CF[3], ax = ax[3], fraction = 0.01,  pad = 0.03, orientation = "horizontal", label = "[Sv]")
fig


fig, ax = plt.subplots(1, 3, figsize=(15,20), 
                       subplot_kw=Dict("projection"=> proj0), sharey = true)
bounds = [2, 2,0.5]

[a.coastlines() for a in ax]
[a.set_extent((120, 250, 0, 56),crs=projPC) for a in ax]
for (i, var) in enumerate(vars)
    E = LLCcropC360(E_dict[var], γ)[1:1:end, 1:1:end]; N = LLCcropC360(N_dict[var], γ)[1:1:end, 1:1:end];
    x = LLCcropC360(λ_wrap, γ)[1:1:end, 1:1:end]; y = LLCcropC360(ϕ, γ)[1:1:end, 1:1:end]
    E[E .== 0] .= NaN; N[N .== 0] .= NaN

    ax[i].streamplot(x=x, y=y, u=1e-6 .* E, v=1e-6 .* N, transform = projPC)
    ax[i].set_title(var)
end
fig
# ax.set_extent((-180, 180, -70, 56),crs=projPC)
fig.colorbar(CF[1], ax = ax[1:2], fraction = 0.01, pad = 0.03, orientation = "horizontal", label = "[Sv]")
fig.colorbar(CF[3], ax = ax[3], fraction = 0.01,  pad = 0.03, orientation = "horizontal", label = "[Sv]")
fig


fig, ax = plt.subplots(1, 3, figsize=(15,20), 
                       subplot_kw=Dict("projection"=> proj0), sharey = true)
bounds = [2, 2,0.5]

[a.coastlines() for a in ax]
[a.set_extent((120, 250, 0, 56),crs=projPC) for a in ax]
for (i, var) in enumerate(vars)
    CF[i] = OHC_helper.pcolormesh_ma(ax[i], λ_wrap, ϕ, 1e-6 .* E_dict[var], 
                        bounds[i], cm.balance, projPC, true)
    ax[i].set_title(var)
end
# ax.set_extent((-180, 180, -70, 56),crs=projPC)
fig.colorbar(CF[1], ax = ax[1:2], fraction = 0.01, pad = 0.03, orientation = "horizontal", label = "[Sv]")
fig.colorbar(CF[3], ax = ax[3], fraction = 0.01,  pad = 0.03, orientation = "horizontal", label = "[Sv]")
fig

maximum(N_dict["iter129_bulkformula"])