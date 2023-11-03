include("../../src/intro.jl")

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

lvls = findall( -3000 .<= -z[:].<= -2000)
mid_depths(x) = vec(sum(Float32, x[lvls, :] .* ΔV[lvls], dims = 1) / sum(Float32, ΔV[lvls]))

tecco = 1992+1/24:1/12:2018;
get_datafiles(expname, key) = filter(x -> occursin("data",x),searchdir(diagpath[expname],key) )

function get_temperature(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid, cell_volumes)

    datafilelist_θ  = get_datafiles(expname, "state_3d_set1")
    ΔV = lateral_sum(cell_volumes)
    nt = length(datafilelist_θ); nz = 50
    println(nt, " months available")
    θ_avg = zeros(Float32, nz, nt)
    ma_template = MeshArray(γ,Float32,50)
    @time for tt = 1:nt
        println(tt)
        fnameθ = datafilelist_θ[tt]
        @time θ = γ.read(diagpath[expname]*fnameθ,ma_template)

        θ_avg[:, tt] .= lateral_sum(θ .* cell_volumes); 
        θ_avg[:, tt] .= θ_avg[:, tt] ./ ΔV
    end

    return θ_avg
end

adjust_exps = Dict()

adjust_exps["only_buoyancy"] = get_temperature(diagpath, "only_buoyancy", γ, cell_volumes)
adjust_exps["only_sfc"] = get_temperature(diagpath, "only_sfc", γ, cell_volumes)
adjust_exps["iter0"] = get_temperature(diagpath, "iter0_bulkformula", γ, cell_volumes)
adjust_exps["only_wind"] = get_temperature(diagpath, "only_wind", γ, cell_volumes)

fig, ax = plt.subplots(figsize = (12.5, 10))
ft = size(adjust_exps["only_wind"], 2)
ax.plot(tecco[1:ft], mid_depths(adjust_exps["only_wind"])[1:ft], label = "WIND Adjust", color = "darkblue", linestyle = "--", linewidth = 2.5)
ax.plot(tecco[1:ft], mid_depths(adjust_exps["only_buoyancy"])[1:ft], label = "BUOYANCY Adjust", color = "darkred", linestyle = "--", linewidth = 2.5)
ax.plot(tecco[1:ft], mid_depths(adjust_exps["only_sfc"])[1:ft], label = "FORCING Adjust", color = exp_colors["only_sfc"])
ax.plot(tecco[1:ft], mid_depths(adjust_exps["iter0"])[1:ft], label = "CTRL", color = exp_colors["iter0_bulkformula"])
ax.legend(frameon = false)
ax.set_title("Mid-Depth North Pacific Temperature (z = 2-3km)", fontweight = "bold")
ax.set_xlabel("time", fontweight = "bold"); ax.set_ylabel("deg C", fontweight = "bold")
fig
fig.savefig(plotsdir("native/sensitivity_exps/North_Pacific_Atm_Exps_Wind.png"))
