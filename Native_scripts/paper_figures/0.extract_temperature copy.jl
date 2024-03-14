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
region = "PAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)

cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)
ΔV = lateral_sum(cell_volumes)

lvls = findall( -3000 .<= -z[:].<= -2000)
mid_depths(x) = vec(sum(Float32, x[lvls, :] .* ΔV[lvls], dims = 1) / sum(Float32, ΔV[lvls]))

tecco = 1992+1/24:1/12:2018; nz = 50
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
adjust_exps["iter0_bulkformula"] = get_temperature(diagpath, "iter0_bulkformula", γ, cell_volumes)
adjust_exps["only_kappa"] = get_temperature(diagpath, "only_kappa", γ, cell_volumes)
adjust_exps["only_init"] = get_temperature(diagpath, "only_init", γ, cell_volumes)
adjust_exps["only_wind"] = get_temperature(diagpath, "only_wind", γ, cell_volumes)
adjust_exps["only_buoyancy"] = get_temperature(diagpath, "only_buoyancy", γ, cell_volumes)
adjust_exps["iter129_bulkformula"] = get_temperature(diagpath, "iter129_bulkformula", γ, cell_volumes)

jldsave(datadir(region * "_temperature_sens_exps.jld2"), adjust_exps= adjust_exps)

adjust_exps_PAC = jldopen(datadir(region * "_temperature_sens_exps.jld2"))["adjust_exps"]
adjust_exps_NPAC = jldopen(datadir("NPAC" * "_temperature_sens_exps.jld2"))["adjust_exps"]


adjust_exps_trends = Dict()
adjust_exps_trends = Dict(key => compute_depth_trends(adjust_exps[key]) for key in keys(adjust_exps))

fig, ax = plt.subplots(figsize = (10, 10))
ax.plot(adjust_exps_trends["iter0_bulkformula"], z, label = "CTRL")
ax.plot(adjust_exps_trends["iter129_bulkformula"], z, label = "FULL")

ax.plot(adjust_exps_trends["only_init"], z, label = "INIT")
ax.plot(adjust_exps_trends["only_kappa"], z, label = "MIX")
ax.plot(adjust_exps_trends["only_wind"], z, label = "WIND")
ax.plot(adjust_exps_trends["only_buoyancy"], z, label = "BUOYANCY")

rect = patches.Rectangle((-10*0.1, 2000), 10*0.2, 1000, linewidth=3, edgecolor="none",facecolor="black", alpha = 0.05)
ax.add_patch(rect)
ax.axvline(0, color = "black", zorder = 0, alpha = 0.55, linewidth = 1)
ax.legend(frameon = false)
ax.set_ylim(1500, 5000); ax.invert_yaxis()
ax.set_xlim(-1, 1)

fig

region = "NPAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)

cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)
NPAC_ΔV = lateral_sum(cell_volumes)

Δθ = adjust_exps_PAC["iter129_bulkformula"] .- adjust_exps_PAC["iter0_bulkformula"]
Δθ_NPAC = adjust_exps_NPAC["iter129_bulkformula"] .- adjust_exps_PAC["iter0_bulkformula"]
Δθ_NPAC = Δθ_NPAC .* NPAC_ΔV ./ ΔV
Δθ_ESPAC = Δθ .- Δθ_NPAC

E,F = trend_matrices(Float32.(tecco))
compute_depth_trends(x) = (100 * 10) .* (x * F[2, :])

adjust_exps_trends = Dict()
adjust_exps_trends = Dict(key => compute_depth_trends(adjust_exps[key]) for key in keys(adjust_exps))

fig, ax = plt.subplots(figsize = (10, 10))
ax.plot(compute_depth_trends(Δθ), z, label = "Δθ")
ax.plot(compute_depth_trends(Δθ_NPAC), z, label = "Δθ_NPAC")
ax.plot(compute_depth_trends(Δθ_ESPAC), z, label = "Δθ_ESPAC")
rect = patches.Rectangle((-10*0.1, 2000), 10*0.2, 1000, linewidth=3, edgecolor="none",facecolor="black", alpha = 0.05)
ax.add_patch(rect)
ax.axvline(0, color = "black", zorder = 0, alpha = 0.55, linewidth = 1)
fig
ax.legend(frameon = false)
ax.set_ylim(1500, 5000); ax.invert_yaxis()
ax.set_xlim(-1, 1)

fig
# fig.savefig(plotsdir("native/sensitivity_exps/North_Pacific_Temp_Trend.png"))
# fig.savefig(plotsdir("native/"), bbox_inches = "tight")



compute_depth_diff(x) = 100 * (x[:, end] .- x[:, 1])
adjust_exps_trends = Dict()
adjust_exps_trends = Dict(key => compute_depth_diff(adjust_exps[key]) for key in keys(adjust_exps))
exps = ["Initial", "Kappa", "Forcing"]
[adjust_exps_trends[key] = adjust_exps_trends[key] .- adjust_exps_trends["CTRL"] for key in exps] #remove the effect of the baseline
adjust_exps_trends["SUM"] = adjust_exps_trends["CTRL"] .+ adjust_exps_trends["Initial"] .+ 
                            adjust_exps_trends["Kappa"] .+ adjust_exps_trends["Forcing"]

fig, ax = plt.subplots(figsize = (10, 10))
ax.plot(adjust_exps_trends["CTRL"], z, label = "CTRL (Iteration 0)")
ax.plot(adjust_exps_trends["Initial"], z, label = "INIT Effect")
ax.plot(adjust_exps_trends["Kappa"], z, label = "MIXING Effect")
ax.plot(adjust_exps_trends["Forcing"], z, label = "FORCING Effect")
ax.plot(adjust_exps_trends["SUM"], z, label = "SUM", color = "k", linestyle = "--", alpha = 0.6)
ax.plot(adjust_exps_trends["FULL"], z, label = "FULL (Iteration 129)", color = "k", linewidth = 3)
rect = patches.Rectangle((-3, 2000), 6, 1000, linewidth=3, edgecolor="none",facecolor="black", alpha = 0.05)
ax.add_patch(rect)
ax.axvline(0, color = "black", zorder = 0, alpha = 0.55, linewidth = 1)
ax.legend(frameon = false)
ax.set_title("North Pacific Temperature Change (2017 minus 1992)", fontweight = "bold")
ax.set_xlabel("cK", fontweight = "bold"); ax.set_ylabel("Depth [m]", fontweight = "bold")
ax.set_xlim(-3, 3)
ax.set_ylim(1500, 5000); ax.invert_yaxis()
fig
fig.savefig(plotsdir("native/sensitivity_exps/North_Pacific_Temp_Diff.png"))




fig, ax = plt.subplots(figsize = (10, 10))
ax.plot((adjust_exps_trends["Initial"] .+ adjust_exps_trends["CTRL"])./10, z, label = "INIT")
ax.plot(adjust_exps_trends["FULL"] ./10, z, label = "INIT")

rect = patches.Rectangle((-0.1, 2000), 0.2, 1000, linewidth=3, edgecolor="none",facecolor="black", alpha = 0.05)
ax.add_patch(rect)
ax.axvline(0, color = "black", zorder = 0, alpha = 0.55, linewidth = 1)
ax.legend(frameon = false)
ax.set_title("North Pacific Temperature Trend Profile", fontweight = "bold")
ax.set_xlabel("cK per decade", fontweight = "bold"); ax.set_ylabel("Depth [m]", fontweight = "bold")
ax.set_xlim(-.1, .1)
ax.set_ylim(1500, 5000); ax.invert_yaxis()
fig