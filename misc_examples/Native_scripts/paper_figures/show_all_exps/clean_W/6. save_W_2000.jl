include("../../../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, DSP, PyCall, 
    RollingFunctions, LinearAlgebra, Statistics, FFTW, LaTeXStrings
import PyPlot as plt 
import NaNMath as nm

include(srcdir("config_exp.jl"))
@pyimport matplotlib.patches as patches

include(srcdir("MeshArraysPlots.jl"))

tecco = 1992+1/24:1/12:2018
nz = 50

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
λ_wrap = wrap_λ(λ)
ocean_mask = wet_pts(Γ)
cmo = pyimport("cmocean.cm")

region = "PAC"
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)

cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)

region = "NPAC"
NPAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)
mskC, mskW, mskS = get_msk(Γ)
include(srcdir("plot_and_dir_config.jl"))
function get_transports(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid)

    filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
    datafilelist_τ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    nt = length(datafilelist_τ);
    W_avg = MeshArray(γ, Float32, 312)
    fill!(W_avg, 0.0)
    @time for tt = 1:312
        println(tt)
        fnameuvw = datafilelist_τ[tt]

        u, v, w, Ub, Vb, Wb = extract_eulerian_and_bolus_velocities(diagpath, expname, 
        fnameuvw, γ, Γ, mskC, mskW, mskS)
        W = w .+ Wb
        W_avg.f[:, tt] .= 1 .* W.f[:, 38]
    end

    return W_avg
end

W_mean = Dict()
W_mean["only_wind"]   = get_transports(diagpath, "only_wind", γ)
W_mean["only_init"]   = get_transports(diagpath, "only_init", γ)
W_mean["only_kappa"]   = get_transports(diagpath, "only_kappa", γ)
W_mean["iter0_bulkformula"]   = get_transports(diagpath, "iter0_bulkformula", γ)

# W_mean["diff_wind"] = W_mean["only_wind"] .- W_mean["iter0_bulkformula"]
# W_mean["diff_init"] = W_mean["only_init"] .- W_mean["iter0_bulkformula"]
# W_mean["diff_kappa"] = W_mean["only_kappa"] .- W_mean["iter0_bulkformula"]

jldsave(datadir("W_2000_all_times.jld2"), W = W_mean)
# jldopen(datadir("W_2000_all_times.jld2"))["W"]