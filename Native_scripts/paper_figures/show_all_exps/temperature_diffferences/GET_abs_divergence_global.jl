include("../../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools
import DataFrames as DF
import PyPlot as plt

include(srcdir("plot_and_dir_config.jl"))
@pyimport matplotlib.patches as patches

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)
PAC_msk = ocean_mask

cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)
ΔV = lateral_sum(cell_volumes)
lvls = findall( -3000 .<= -z[:].<= -2000)

tecco = 1992+1/24:1/12:2018; nz = 50
get_datafiles(expname, key) = filter(x -> occursin("data",x),searchdir(diagpath[expname],key) )

γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
Γ=GridLoad(γ; option="full")

lvls = findall( -3000 .<= -z[:].<= -2000)
sum_volumes = sum( cell_volumes[:, lvls])

function get_temperature(diagpath::Dict{String, String}, 
    expbase::String,  expcomp::String, γ::gcmgrid, cell_volumes)

    datafilelist_θ  = get_datafiles(expbase, "state_3d_set1")

    ΔV = lateral_sum(cell_volumes)
    nt = length(datafilelist_θ); nz = 50
    println(nt, " months available")
    θ_avg = zeros(Float32, nt)
    ma_template = MeshArray(γ,Float32,50)
    ts = Any[]
    @time for tt = 1:nt
        println(tt)
        fnameθ = datafilelist_θ[tt]
        @time θ_base = γ.read(diagpath[expbase]*fnameθ,ma_template)
        @time θ_comp = γ.read(diagpath[expcomp]*fnameθ,ma_template)
        Δθ = θ_base .- θ_comp
        Δθ = Δθ[:, lvls]
        Δθ = abs.(Δθ)
        tmp = sum(Δθ .* cell_volumes[:, lvls]) /  sum_volumes
        θ_avg[tt] = 1 * tmp
    end

    return θ_avg
end
[println(k) for k in keys(diagpath)]
adjust_exps = Dict()

adjust_exps["only_kappa"] = get_temperature(diagpath, "iter0_bulkformula", "only_kappa", γ, cell_volumes)

adjust_exps["only_init"] = get_temperature(diagpath, "iter0_bulkformula", "only_init", γ, cell_volumes)
adjust_exps["noinitadjust"] = get_temperature(diagpath, "iter129_bulkformula", "noinitadjust", γ, cell_volumes)

adjust_exps["only_sfc"] = get_temperature(diagpath, "iter0_bulkformula", "only_sfc", γ, cell_volumes)
adjust_exps["nosfcadjust"] = get_temperature(diagpath, "iter129_bulkformula", "nosfcadjust", γ, cell_volumes)


jldsave(datadir("θ_ABS_DIV_GLOBAL_all_times.jld2"), adjust_exps = adjust_exps)
