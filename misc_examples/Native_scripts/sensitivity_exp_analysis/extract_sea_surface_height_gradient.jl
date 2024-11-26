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
cell_depths = get_cell_thickness(ocean_mask, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)
tecco = 1992+1/24:1/12:2018; nz = 50
get_datafiles(expname, key) = filter(x -> occursin("data",x),searchdir(diagpath[expname],key) )

inv_depths = vertical_sum(cell_depths) .* 1; inv_depths = 1 ./ inv_depths
[inv_depths.f[ff][(!isfinite).(inv_depths.f[ff])] .= 0.0 for ff = 1:5]
bottom_mask = vertical_sum(cell_volumes[:, 38:42])
for ff = 1:5 
    bottom_mask.f[ff][bottom_mask.f[ff] .> 0] .= 1.0
end

height_3000 = vertical_sum(cell_volumes[:, 38:end]) ./ area
function get_eta(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid, cell_volumes)

    datafilelist_θ  = get_datafiles(expname, "state_2d_set1")
    nt = length(datafilelist_θ); nz = 50;
    println(nt, " months available")
    η_avg = MeshArray(γ,Float32); fill!(η_avg, 0.0)
    ma_template = MeshArray(γ,Float32)
    @time for tt = 1:nt
        println(tt)
        fnameθ = datafilelist_θ[tt]
        @time η = γ.read(diagpath[expname]*fnameθ,ma_template)
        η_avg .+= η ./nt
    end

    return η_avg
end

adjust_exps = Dict()
adjust_exps["CTRL"] = get_eta(diagpath, "iter0_bulkformula", γ, cell_volumes)
adjust_exps["Kappa"] = get_eta(diagpath, "only_kappa", γ, cell_volumes)
adjust_exps["Initial"] = get_eta(diagpath, "only_init", γ, cell_volumes)
adjust_exps["Forcing"] = get_eta(diagpath, "only_sfc", γ, cell_volumes)
adjust_exps["FULL"] = get_eta(diagpath, "iter129_bulkformula", γ, cell_volumes)

# eta_trends_effect["SUM"] = eta_trends_effect["Kappa"] .+ eta_trends_effect["Initial"] .+ eta_trends_effect["Forcing"] .+ eta_trends_effect["CTRL"]
conv_adjust_exps = Dict()
DXCsm=1.5*Γ.DXC; DYCsm=1.5*Γ.DYC;
for expt in keys(adjust_exps)
    η = adjust_exps[expt]

    s = η .* inv_depths
    s .+= 1
    dataU, dataV = gradient(s .* η, Γ, true)

    conv_adjust_exps[expt] = UVtoUEVN(dataU, dataV, Γ)[1]
    conv_adjust_exps[expt] = conv_adjust_exps[expt] .* bottom_mask
end

conv_adjust_exps["Initial"] = conv_adjust_exps["Initial"] .- conv_adjust_exps["CTRL"]
conv_adjust_exps["Kappa"] = conv_adjust_exps["Kappa"] .- conv_adjust_exps["CTRL"]
conv_adjust_exps["Forcing"] = conv_adjust_exps["Forcing"] .- conv_adjust_exps["CTRL"]
conv_adjust_exps["SUM"] =  conv_adjust_exps["Initial"] .+ conv_adjust_exps["Kappa"] .+ conv_adjust_exps["Forcing"] .+ conv_adjust_exps["CTRL"]
# conv_adjust_exps["Kappa"] = MeshArrays.smooth(conv_adjust_exps["Kappa"] ,DXCsm,DYCsm,Γ)
# conv_adjust_exps["Forcing"] = MeshArrays.smooth(conv_adjust_exps["Forcing"] ,DXCsm,DYCsm,Γ)
# conv_adjust_exps["Initial"] = MeshArrays.smooth(conv_adjust_exps["Initial"] ,DXCsm,DYCsm,Γ)

lin_exps = ["CTRL", "Initial", "Kappa", "Forcing", "FULL", "SUM"]; nexps = length(lin_exps)
plot_labels_effect = ["CTRL (Iteration 0)", "INIT Effect", "MIXING Effect", "FORCING Effect", "FULL (Iteration 129)", "SUM"]
plot_labels_effect = Dict(lin_exps[i] => plot_labels_effect[i] for i in 1:nexps)
lin_exps = ["Initial", "CTRL", "Kappa", "FULL", "Forcing", "SUM"]; nexps = length(lin_exps)

fig, axs = plt.subplots(2, 3, figsize=(17,12), subplot_kw=Dict("projection"=> proj0))
CF = Any[]
tmp = maximum(conv_adjust_exps["CTRL"][4])
bnds = 86400 * 365
for (i, exp) in enumerate(lin_exps)
    ax = axs[i]
    data = conv_adjust_exps[exp]
    # mult = 1000 * 86400 * 365
    for ff = 1:5
        cf = ax.pcolormesh(λ[ff], ϕ[ff],  -data[ff], transform=projPC, 
        cmap = cmo.balance, vmin = -tmp/30, vmax = tmp/30)   
        push!(CF, cf)
    end
    ax.coastlines(resolution="110m")
    ax.set_extent((120, 295, -40, 70),crs=projPC)
    ax.set_title(plot_labels_effect[exp])
    gl = ax.gridlines(crs=projPC, draw_labels=true,
                        linewidth=2, color="gray", alpha=0, linestyle="--")
    gl.top_labels = false
    gl.right_labels = false
end
fig.colorbar(CF[1], ax = axs[:], fraction = 0.04, 
label = "-∂η/∂λ [unitless]", orientation = "horizontal")
fig
fig.savefig(plotsdir("native/sensitivity_exps/minus∂η∂x.png"), bbox_inches = "tight")
