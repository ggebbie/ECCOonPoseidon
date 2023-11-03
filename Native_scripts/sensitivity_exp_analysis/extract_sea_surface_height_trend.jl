include("../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools
import PyPlot as plt

include(srcdir("plot_and_dir_config.jl"))
@pyimport matplotlib.patches as patches
@pyimport cmocean.cm as cmo

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)
region = "NPAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)

cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)
ΔV = lateral_sum(cell_volumes)

tecco = 1992+1/24:1/12:2018; nz = 50
E, F = trend_matrices(tecco)
get_datafiles(expname, key) = filter(x -> occursin("data",x),searchdir(diagpath[expname],key) )

function get_eta(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid)

    datafilelist_θ  = get_datafiles(expname, "state_2d_set1")
    nt = length(datafilelist_θ);
    println(nt, " months available")
    η_avg = MeshArray(γ,Float32); fill!(η_avg, 0.0)
    ma_template = MeshArray(γ,Float32)
    @time for tt = 1:nt
        println(tt)
        fnameθ = datafilelist_θ[tt]
        @time η = γ.read(diagpath[expname]*fnameθ,ma_template)
        for ff in 1:5
            η_avg.f[ff] .+= F[2, tt] .* η.f[ff]
        end
    end

    return η_avg
end

adjust_exps = Dict()
adjust_exps["CTRL"] = get_eta(diagpath, "iter0_bulkformula", γ)
adjust_exps["Kappa"] = get_eta(diagpath, "only_kappa", γ)
adjust_exps["Initial"] = get_eta(diagpath, "only_init", γ)
adjust_exps["Forcing"] = get_eta(diagpath, "only_sfc", γ)
adjust_exps["FULL"] = get_eta(diagpath, "iter129_bulkformula", γ)

eta_trends_effect = deepcopy(adjust_exps)
eta_trends_effect["Initial"] = adjust_exps["Initial"] .- adjust_exps["CTRL"]
eta_trends_effect["Kappa"] = adjust_exps["Kappa"] .- adjust_exps["CTRL"]
eta_trends_effect["Forcing"] = adjust_exps["Forcing"] .- adjust_exps["CTRL"]

eta_trends_effect["SUM"] = eta_trends_effect["Kappa"] .+ eta_trends_effect["Initial"] .+ eta_trends_effect["Forcing"] .+ eta_trends_effect["CTRL"]


lin_exps = ["CTRL", "Initial", "Kappa", "Forcing", "FULL", "SUM"]; nexps = length(lin_exps)
plot_labels_effect = ["CTRL (Iteration 0)", "INIT Effect", "MIXING Effect", "FORCING Effect", "FULL (Iteration 129)", "SUM"]
plot_labels_effect = Dict(lin_exps[i] => plot_labels_effect[i] for i in 1:nexps)
lin_exps = ["Initial", "CTRL", "Kappa", "FULL", "Forcing", "SUM"]; nexps = length(lin_exps)

fig, axs = plt.subplots(2, 3, figsize=(17,12), subplot_kw=Dict("projection"=> proj0))
CF = Any[]
tmp = adjust_exps["FULL"] 
bnds = maximum(abs.(tmp))

for (i, exp) in enumerate(lin_exps)
    ax = axs[i]
    CF = Any[]
    data = eta_trends_effect[exp];
    for ff = 1:5
        cf = ax.pcolormesh(λ[ff], ϕ[ff],  1000 .* data[ff], transform=projPC, 
        cmap = cmo.balance, vmin = -1, vmax = 8)   
        push!(CF, cf)
    end
    ax.coastlines(resolution="110m")
    ax.set_extent((120, 295, -70, 56),crs=projPC)
    ax.set_title(plot_labels_effect[exp])
    gl = ax.gridlines(crs=projPC, draw_labels=true,
                        linewidth=2, color="gray", alpha=0, linestyle="--")
    gl.top_labels = false
    gl.right_labels = false
end
fig

# close("all")