include("../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools
import PyPlot as plt

include(srcdir("plot_and_dir_config.jl"))
@pyimport matplotlib.patches as patches

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)
region = "NPAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)

cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)
ΔV = lateral_sum(cell_volumes)

area_mask = sum(area .* PAC_msk)
sum(cell_volumes[:, lvls])
inv_D = 1 .* cell_depths; inv_D[findall(inv_D .== 0.0)] = Inf
inv_D = 1 ./ inv_D

lvls = findall( -3000 .<= -z[:].<= -2000)
mid_depths(x) = vec(sum(Float32, x[lvls, :] .* ΔV[lvls], dims = 1) / sum(Float32, ΔV[lvls]))

tecco = 1992+1/24:1/12:2018; nz = 50
get_datafiles(expname, key) = filter(x -> occursin("data",x),searchdir(diagpath[expname],key) )

function get_temperature(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid, cell_volumes)

    datafilelist_θ  = get_datafiles(expname, "state_3d_set1")
    nt = 312; nz = 50
    println(nt, " months available")
    θ_avg = zeros(nt); 
    ma_template = MeshArray(γ,Float32,50)
    θ_interp = MeshArray(γ,Float32,50); fill!(θ_interp, 0.0)
    θ_conv = MeshArray(γ,Float32,50); fill!(θ_conv, 0.0)

    @time for tt = 1:nt
        println(tt)
        fnameθ = datafilelist_θ[tt]
        @time θ = γ.read(diagpath[expname]*fnameθ,ma_template)
        interpolate_to_vertical_faces!(θ, θ_interp)
        calc_W_conv3D!(θ_interp, θ_conv)
        θ_conv = -θ_conv .* inv_D
        vol_avg_dθ = sum(θ_conv[:, lvls] .* cell_volumes[:, lvls]) / sum(cell_volumes[:, lvls])
        θ_avg[tt] = vol_avg_dθ 
    end

    return θ_avg
end

adjust_exps = Dict()
adjust_exps["iter129_bulkformula"] = get_temperature(diagpath, "iter129_bulkformula", γ, cell_volumes)
adjust_exps["only_kappa"] = get_temperature(diagpath, "only_kappa", γ, cell_volumes)
adjust_exps["only_init"] = get_temperature(diagpath, "only_init", γ, cell_volumes)
adjust_exps["only_sfc"] = get_temperature(diagpath, "only_sfc", γ, cell_volumes)
adjust_exps["only_wind"] = get_temperature(diagpath, "only_wind", γ, cell_volumes)
adjust_exps["only_buoyancy"] = get_temperature(diagpath, "only_buoyancy", γ, cell_volumes)
adjust_exps["iter0_bulkformula"] = get_temperature(diagpath, "iter0_bulkformula", γ, cell_volumes)

jldsave(datadir(region * "_temperature_gradient_sens_exps.jld2"), adjust_exps= adjust_exps)

fig, ax = plt.subplots(figsize = (10, 7.5))
ax.plot(tecco, adjust_exps["iter0_bulkformula"], label = "CTRL (Iteration 0)")
ax.plot(tecco, adjust_exps["only_init"], label = "INIT")
ax.plot(tecco, adjust_exps["only_kappa"], label = "MIXING")
ax.plot(tecco, adjust_exps["only_sfc"], label = "FORCING")
ax.plot(tecco, adjust_exps["iter129_bulkformula"], label = "FULL (Iteration 129)", color = "k")
ax.legend(frameon = false)
ax.set_title("Mid-Depth North Pacific Temperature (z = 2-3km)", fontweight = "bold")
ax.set_xlabel("time", fontweight = "bold"); ax.set_ylabel("deg C", fontweight = "bold")
fig
# fig.savefig(plotsdir("native/sensitivity_exps/North_Pacific_Sens_Exps.png"))
