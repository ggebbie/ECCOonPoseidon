include("../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools
import PyPlot as plt
@pyimport cmocean.cm as cmos

include(srcdir("plot_and_dir_config.jl"))
@pyimport matplotlib.patches as patches

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)

cell_depths = get_cell_thickness(ocean_mask, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)

lvls = findall( -3000 .<= -z[:].<= -2000); suffix = "2to3"

tecco = 1992+1/24:1/12:2018; nz = 50
E,F = trend_matrices(tecco)

get_datafiles(expname, key) = filter(x -> occursin("data",x),searchdir(diagpath[expname],key) )

function get_temperature_trends(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid, cell_volumes)

    datafilelist_θ  = get_datafiles(expname, "state_3d_set1")
    nt = length(datafilelist_θ); nz = 50
    println(nt, " months available")

    β = MeshArray(γ,Float32); fill!(β, 0.0)
    flat_volumes = vertical_sum(cell_volumes[:, lvls])
    for ff=1:5
        flat_volumes.f[ff][iszero.(flat_volumes.f[ff])] .= NaN
    end
    ma_template = MeshArray(γ,Float32,50)
    @time for tt = 1:nt
        println(tt)
        fnameθ = datafilelist_θ[tt]
        @time θ = γ.read(diagpath[expname]*fnameθ,ma_template)
        sθH = vertical_sum(θ[:, lvls] .* cell_volumes[:, lvls]) ./ flat_volumes
        for ff in 1:5
            β[ff] .+= F[2,tt] .* sθH[ff] 
        end
    end
    return β
end


adjust_exps = Dict()
adjust_exps["CTRL"] = get_temperature_trends(diagpath, "iter0_bulkformula", γ, cell_volumes)
adjust_exps["Kappa"] = get_temperature_trends(diagpath, "only_kappa", γ, cell_volumes)
adjust_exps["Initial"] = get_temperature_trends(diagpath, "only_init", γ, cell_volumes)
adjust_exps["Forcing"] = get_temperature_trends(diagpath, "only_sfc", γ, cell_volumes)
adjust_exps["FULL"] = get_temperature_trends(diagpath, "iter129_bulkformula", γ, cell_volumes)

savename = datadir("native/_THETA_spatial_trend_" * suffix * ".jld2")
jldsave(savename, theta_trends = adjust_exps)
