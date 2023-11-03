include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, DSP, PyCall
import PyPlot as plt 
include(srcdir("config_exp.jl"))

include(srcdir("MeshArraysPlots.jl"))

tecco = 1992+1/24:1/12:2018
nz = 50

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
λ_wrap = wrap_λ(λ)
ocean_mask = wet_pts(Γ)

region = "NPAC"
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)

include(srcdir("plot_and_dir_config.jl"))

NW_PAC = region_mask(ocean_mask, λ_wrap, ϕ, (0, 90), (0, 190))
NE_PAC = region_mask(ocean_mask, λ_wrap, ϕ, (0, 90), (190, 360))

cs, sn = get_cs_and_sn(γ)

reg_mask = LLCcropC(PAC_msk,γ)

function get_transports(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid)

    filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
    datafilelist_τ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    nt = length(datafilelist);

    ws = zeros(size(reg_mask)..., nt)
    curlτ = zeros(size(reg_mask)..., nt)

    @time for tt = 1:nt
        println(tt)
        fname = datafilelist[tt]
        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        _, _, w = extract_eulerian_velocities(diagpath, expname, fname, γ)

        ws[:, :, tt] .= LLCcropC(w[:, 43], γ)


        Tname = datafilelist_τ[tt]
        τx, τy = extract_ocnTAU(diagpath, expname , Tname, γ)
        τcurl = MeshArrays.curl(τx .* 10000, τy .* 10000, Γ)
        curlτ[:, :, tt].= LLCcropC(τcurl, γ) ./ 10000

    end

    return ws, curlτ
end
        
expname = "only_kappa"
W_mean_sfc, tau_mean_sfc= get_transports(diagpath, expname, γ); 

expname = "iter0_bulkformula"
W_mean_0, tau_mean_0 = get_transports(diagpath, expname, γ); 

force_eff_w = W_mean_sfc .- W_mean_0; force_eff_w = 1e-4 .* force_eff_w .* LLCcropC(area, γ)
force_eff_tau = tau_mean_sfc .- tau_mean_0;

reg_ϕ = LLCcropC(ϕ, γ); reg_λ = LLCcropC(λ, γ)

fig, axs = plt.subplots(2, 1, figsize=(20,10), subplot_kw=Dict("projection"=> proj0))
[ax.set_extent((120, 240, 20, 56),crs=projPC) for ax in axs]
[ax.coastlines()  for ax in axs]
[ax.gridlines(crs=projPC, draw_labels=true, linewidth=2, color="gray", alpha=0, linestyle="--") for ax in axs]

ax = axs[1]
data = mean(force_eff_w, dims = 3)[:, :, 1]
data[data .== 0.0] .= NaN
CM = ax.pcolormesh(reg_λ, reg_ϕ, data .* 1e-2, cmap = cmo.balance, transform = projPC, 
vmin = -2 * 1e-2, vmax = 2* 1e-2)
ax.scatter(reg_λ[330, end - 30], reg_ϕ[330, end - 30], transform = projPC, color = "red")
# ax.scatter(reg_λ[15, end - 18], reg_ϕ[15, end - 18], transform = projPC, color = "red")
ax.scatter(reg_λ[34, end - 13], reg_ϕ[34, end - 13], transform = projPC, color = "red")
fig.colorbar(CM, ax = ax, orientation = "vertical", label = "Vertical Transport Anomaly [Sv]")

ax = axs[2]
data = mean(force_eff_tau, dims = 3)[:, :, 1]
data = data .* reg_mask; data[data .== 0.0] .= NaN
CM = ax.pcolormesh(reg_λ, reg_ϕ, data, cmap = cmo.curl, transform = projPC, 
vmin = -1e-12, vmax = 1e-12)
ax.scatter(reg_λ[330, end - 30], reg_ϕ[330, end - 30], transform = projPC, color = "red")
# ax.scatter(reg_λ[15, end - 18], reg_ϕ[15, end - 18], transform = projPC, color = "red")
ax.scatter(reg_λ[34, end - 13], reg_ϕ[34, end - 13], transform = projPC, color = "red")
fig.colorbar(CM, ax = ax, orientation = "vertical", label = "Wind Stress Curl Anomaly[N/m]")

# correlation_map = zeros(360, 180)
# cov(force_eff_w, force_eff_tau, dims = 3)
# data1 = mean(force_eff_w, dims = 3)[:, :, 1]
# data2 = mean(force_eff_tau, dims = 3)[:, :, 1]

fig


fig, axs = plt.subplots(2, 1, figsize=(15,12), sharex = true)
coords = string(round.(reg_λ[330, end - 30], digits = 2)) *" E , " * 
string(round.(reg_ϕ[330, end - 30], digits = 2)) *" N"
data1 = vec(force_eff_w[330, end - 30, :])
data2 =  vec(force_eff_tau[330, end - 30, :])
axs[1].plot(tecco, 1e-2 .* data1); 
axs[1].set_title("Vertical Transport Anomaly (z = 3000 m) Due to Wind Adjustment \n " * coords, fontweight = "bold")
axs[1].set_ylabel("[Sv]", fontweight = "bold")
axs[2].plot(tecco, data2)
axs[2].set_title("Curl Anomaly Anomaly Due to Wind Adjustment \n " * coords, fontweight = "bold")
axs[2].set_ylabel("[N/m³]", fontweight = "bold")
axs[2].set_xlabel("time", fontweight = "bold")
fig.tight_layout()
fig

fig, axs = plt.subplots(2, 1, figsize=(15,12), sharex = true)
coords = string(round.(reg_λ[327, end - 28], digits = 2)) *" E , " * 
string(round.(reg_ϕ[327, end - 28], digits = 2)) *" N"
axs[1].plot(tecco, 1e-2 .* vec(force_eff_w[328, end - 28, :])); 
axs[1].set_title("Vertical Transport Anomaly (z = 3000 m) Due to Wind Adjustment \n " * coords, fontweight = "bold")
axs[1].set_ylabel("[Sv]", fontweight = "bold")
axs[2].plot(tecco, vec(force_eff_tau[328, end - 28, :]))
axs[2].set_title("Curl Anomaly Anomaly Due to Wind Adjustment \n " * coords, fontweight = "bold")
axs[2].set_ylabel("[N/m³]", fontweight = "bold")
axs[2].set_xlabel("time", fontweight = "bold")
fig.tight_layout()
fig


fig, axs = plt.subplots(2, 1, figsize=(15,12), sharex = true)
coords = string(round.(reg_λ[327, end - 28], digits = 2)) *" E , " * 
string(round.(reg_ϕ[327, end - 28], digits = 2)) *" N"
axs[1].plot(tecco, 1e-2 .* vec(force_eff_w[328, end - 28, :])); 
axs[1].set_title("Vertical Transport Anomaly Due to Wind Adjustment \n " * coords, fontweight = "bold")
axs[1].set_ylabel("[Sv]", fontweight = "bold")
data = vec(data[328, end - 28, :])
axs[2].plot(tecco, vec(data[328, end - 28, :]))
axs[2].set_title("Curl Anomaly Anomaly Due to Wind Adjustment \n " * coords, fontweight = "bold")
axs[2].set_ylabel("[kg/m³]", fontweight = "bold")
axs[2].set_xlabel("time", fontweight = "bold")
fig.tight_layout()
fig

# fig, axs = plt.subplots(2, 1, figsize=(15,12))
# coords = string(round.(reg_λ[15, end - 18], digits = 2)) *" E , " * 
# string(round.(reg_ϕ[15, end - 18], digits = 2)) *" N"
# axs[1].plot(tecco, 1e-2 .* vec(force_eff_w[15, end - 18, :])); 
# axs[1].set_title("Vertical Transport at 3000 meters Anomaly Due to Wind Adjustment \n " * coords)
# axs[1].set_ylabel("[Sv]")
# axs[2].plot(tecco, vec(force_eff_tau[15, end - 18, :]))
# axs[2].set_title("Curl Anomaly Anomaly Due to Wind Adjustment \n " * coords)
# axs[2].set_ylabel("[kg/m³]")
# axs[2].set_xlabel("time")
# fig


using Statistics
correlation_map = zeros(360, 180)

for ix in 1:360, iy = 1:180
    xx = force_eff_w[ix, iy, :][:]; yy = force_eff_tau[ix, iy, :][:]
    corrxy = cov(xx, yy) / (std(xx) * std(yy))
    correlation_map[ix, iy] = corrxy
end

fig, axs = plt.subplots(2, 1, figsize=(20,10), subplot_kw=Dict("projection"=> proj0))
[ax.set_extent((120, 240, 20, 56),crs=projPC) for ax in axs]
[ax.coastlines()  for ax in axs]
[ax.gridlines(crs=projPC, draw_labels=true, linewidth=2, color="gray", alpha=0, linestyle="--") for ax in axs]

ax = axs[1]
data = correlation_map .* 1
# data[abs.(data) .<= 0.2] .= NaN
CM = ax.pcolormesh(reg_λ, reg_ϕ, data, cmap = cmo.balance, transform = projPC, 
vmin = -1, vmax = 1)
fig.colorbar(CM, ax = axs[1], orientation = "vertical", label = "correlation")
ax.scatter(reg_λ[327, end - 28], reg_ϕ[327, end - 28], transform = projPC, color = "red")
# ax.scatter(reg_λ[15, end - 18], reg_ϕ[15, end - 18], transform = projPC, color = "red")
ax.scatter(reg_λ[34, end - 13], reg_ϕ[34, end - 13], transform = projPC, color = "red")
fig