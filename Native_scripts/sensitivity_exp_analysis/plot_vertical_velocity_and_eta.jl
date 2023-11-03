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
@pyimport cmocean.cm as cmo

NW_PAC = region_mask(ocean_mask, λ_wrap, ϕ, (0, 90), (0, 190))
NE_PAC = region_mask(ocean_mask, λ_wrap, ϕ, (0, 90), (190, 360))

cs, sn = get_cs_and_sn(γ)

reg_mask = LLCcropC(PAC_msk,γ)

function get_transports(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid)

    filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    datafilelist_η = filter(x -> occursin("data",x),filelist) # second filter for "data"

    nt = length(datafilelist); 

    ηreg = zeros(size(reg_mask)..., nt)
    ws = zeros(size(reg_mask)..., nt)
    wtop = zeros(size(reg_mask)..., nt)
    ma_template = MeshArray(γ,Float32)

    @time for tt = 1:nt
        println(tt)
        fname = datafilelist[tt]
        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        _, _, w = extract_eulerian_velocities(diagpath, expname, fname, γ)

        wtop[:, :, tt] .= LLCcropC(w[:, 1], γ)

        ws[:, :, tt] .= LLCcropC(w[:, 43], γ)

        fname_η = datafilelist_η[tt]
        η = γ.read(diagpath[expname]*fname_η,ma_template)

        ηreg[:, :, tt] .= LLCcropC(η, γ)


    end
    var_dict = Dict()
    var_dict["Wbot"] = ws
    var_dict["Wsfc"] = wtop
    var_dict["η"] = ηreg

    return var_dict
end
expname = "only_kappa"
fname = datadir("native/w_oceτ_timeseries_regular_grid_" * expname * ".jld2")
vars = jldopen(fname)["vars"]

W_sfc, Wtop_sfc,η_sfc = (var_dict["Wbot"], var_dict["Wsfc"], var_dict["η"])
vars = get_transports(diagpath, expname, γ); 
jldsave(fname, vars = vars)


jldsave
expname = "iter0_bulkformula"
vars_iter0 = get_transports(diagpath, expname, γ); 
fname = datadir("native/w_oceτ_timeseries_regular_grid_" * expname * ".jld2")
jldsave(fname, vars = vars_iter0)

kappa_eff_w = W_sfc .- W_0; kappa_eff_w = 1e-4 .* kappa_eff_w .* LLCcropC(area, γ)
kappa_eff_wtop = Wtop_sfc .- Wtop_0; kappa_eff_wtop = 1e-4 .* kappa_eff_wtop .* LLCcropC(area, γ)

kappa_eff_η = η_sfc .- η_0;

reg_ϕ = LLCcropC(ϕ, γ); reg_λ = LLCcropC(λ, γ)

fig, axs = plt.subplots(3, 1, figsize=(20,10), subplot_kw=Dict("projection"=> proj0))
[ax.set_extent((120, 240, 20, 56),crs=projPC) for ax in axs]
[ax.coastlines()  for ax in axs]
[ax.gridlines(crs=projPC, draw_labels=true, linewidth=2, color="gray", alpha=0, linestyle="--") for ax in axs]

ax = axs[1]
data = mean(kappa_eff_w, dims = 3)[:, :, 1]; bounds = maximum(abs.(data)) / 4
data[data .== 0.0] .= NaN; 
ax.pcolormesh(reg_λ, reg_ϕ, data, cmap = cmo.balance, transform = projPC, 
vmin = -bounds, vmax = bounds)
ax.scatter(reg_λ[330, end - 30], reg_ϕ[327, end - 30], transform = projPC, color = "red")

ax = axs[2]
data = mean(kappa_eff_wtop, dims = 3)[:, :, 1]; bounds = maximum(abs.(data)) / 4
data[data .== 0.0] .= NaN; 
ax.pcolormesh(reg_λ, reg_ϕ, data, cmap = cmo.balance, transform = projPC, 
vmin = -bounds, vmax = bounds)

ax = axs[3]
data = mean(kappa_eff_η, dims = 3)[:, :, 1]; 
data[2:end, :] .= data[2:end, :] .- data[1:end-1, :]; data[1, :] .= 0.0
bounds = maximum(abs.(data)) / 10
data = data .* reg_mask; data[data .== 0.0] .= NaN
ax.pcolormesh(reg_λ, reg_ϕ, data, cmap = cmo.curl, transform = projPC, 
vmin = -bounds, vmax = bounds)
ax.scatter(reg_λ[330, end - 30], reg_ϕ[327, end - 30], transform = projPC, color = "red")
# ax.scatter(reg_λ[15, end - 18], reg_ϕ[15, end - 18], transform = projPC, color = "red")
# ax.scatter(reg_λ[34, end - 11], reg_ϕ[34, end - 11], transform = projPC, color = "red")

fig


fig, axs = plt.subplots(2, 1, figsize=(15,12))
coords = string(round.(reg_λ[330, end - 30], digits = 2)) *" E , " * 
string(round.(reg_ϕ[330, end - 30], digits = 2)) *" N"
data = 
axs[1].plot(tecco, 1e-2 .* vec(kappa_eff_w[330, end - 30, :])); 
axs[1].set_title("Vertical Transport at 3000 meters Anomaly Due to Wind Adjustment \n " * coords)
axs[1].set_ylabel("[Sv]")

datax = kappa_eff_η .* 1; 
datax[2:end, :, :] .= datax[2:end, :, :] .- datax[1:end-1, :, :]; datax[1, :, :] .= NaN
datax[2:end, :, :] .= datax[2:end, :, :] .- datax[1:end-1, :, :]; datax[1, :, :] .= NaN

datay = kappa_eff_η .* 1; diffy = reg_ϕ .* 1; diffy[:, 2:end, :] .= diffy[:, 2:end, :] .- diffy[:, 1:end-1, :]; datay[:, 1, :] .= NaN

datay[:, 2:end, :] .= (datay[:, 2:end, :] .- datay[:, 1:end-1, :] ); 
datay[:, 1, :] .= NaN;
datay[:, 2:end, :] .= (datay[:, 2:end, :] .- datay[:, 1:end-1, :] ); 
datay[:, 1, :] .= NaN


axs[2].plot(tecco,-vec(datax[330, end - 30, :] .+ datay[330, end - 30, :]))
axs[2].set_title("Curl Anomaly Anomaly Due to Wind Adjustment \n " * coords)
axs[2].set_ylabel("[kg/m³]")
axs[2].set_xlabel("time")
fig