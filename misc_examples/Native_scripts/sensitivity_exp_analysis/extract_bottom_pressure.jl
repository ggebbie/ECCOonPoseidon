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

inv_depths = Γ.Depth .* 1; inv_depths = 1 ./ inv_depths
[inv_depths.f[ff][(!isfinite).(inv_depths.f[ff])] .= 0.0 for ff = 1:5]
tecco = 1992+1/24:1/12:2018; nz = 50
get_datafiles(expname, key) = filter(x -> occursin("data",x),searchdir(diagpath[expname],key) )

function get_P(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid, cell_volumes)

    datafilelist_θ  = get_datafiles(expname, "state_2d_set1")
    nt = length(datafilelist_θ); nz = 50;
    println(nt, " months available")
    η_avg = MeshArray(γ,Float32); fill!(η_avg, 0.0)
    ma_template = MeshArray(γ,Float32, 6)
    @time for tt = 1:nt
        println(tt)
        fnameθ = datafilelist_θ[tt]
        @time Pbot = γ.read(diagpath[expname]*fnameθ,ma_template)[:, 6]
        η_avg .+= Pbot ./nt
    end

    return η_avg
end

adjust_exps = Dict()
adjust_exps["CTRL"] = get_P(diagpath, "iter0_bulkformula", γ, cell_volumes)
adjust_exps["Kappa"] = get_P(diagpath, "only_kappa", γ, cell_volumes)
adjust_exps["Initial"] = get_P(diagpath, "only_init", γ, cell_volumes)
adjust_exps["Forcing"] = get_P(diagpath, "only_sfc", γ, cell_volumes)
adjust_exps["FULL"] = get_P(diagpath, "iter129_bulkformula", γ, cell_volumes)

# eta_trends_effect["SUM"] = eta_trends_effect["Kappa"] .+ eta_trends_effect["Initial"] .+ eta_trends_effect["Forcing"] .+ eta_trends_effect["CTRL"]
lin_exps = ["CTRL", "Initial", "Kappa", "Forcing", "FULL"]; nexps = length(lin_exps)

conv_adjust_exps = Dict()
for expt in lin_exps
    η = adjust_exps[expt]

    # s = η .* inv_depths
    # s .+= 1

    dataU, dataV = gradient(η, Γ, true)

    # dataUU, dataUV = gradient(dataU, Γ, true)
    # dataVU, dataVV = gradient(dataV, Γ, true)

    # conv_adjust_exps[expt] = dataUU .+ dataVV

    conv_adjust_exps[expt] = -MeshArrays.convergence(-dataV.* bottom_mask, dataU .* bottom_mask)


end

conv_adjust_exps["Initial"] = conv_adjust_exps["Initial"] .- conv_adjust_exps["CTRL"]
# conv_adjust_exps["Kappa"] = conv_adjust_exps["Kappa"] .- conv_adjust_exps["CTRL"]
# conv_adjust_exps["Forcing"] = conv_adjust_exps["Forcing"] .- conv_adjust_exps["CTRL"]

plot_labels_effect = ["CTRL (Iteration 0)", "INIT Effect", "MIXING Effect", "FORCING Effect", "FULL (Iteration 129)"]
plot_labels_effect = Dict(lin_exps[i] => plot_labels_effect[i] for i in 1:nexps)
lin_exps = ["Initial", "CTRL", "Kappa", "FULL", "Forcing"]; nexps = length(lin_exps)

fig, axs = plt.subplots(2, 3, figsize=(17,12), subplot_kw=Dict("projection"=> proj0))
CF = Any[]
tmp = conv_adjust_exps["Initial"] .* PAC_msk
bnds = maximum(abs.(tmp[5]))

for (i, exp) in enumerate(lin_exps)
    ax = axs[i]
    CF = Any[]
    data = conv_adjust_exps[exp];
    for ff = 1:5
        cf = ax.pcolormesh(λ[ff], ϕ[ff],  data[ff], transform=projPC, 
        cmap = cmo.balance, vmin = -bnds / 1e6, vmax = bnds / 1e6)   
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
