include("../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools
using NaNMath
import PyPlot as plt

include(srcdir("plot_and_dir_config.jl"))

interp_facs = jldopen(datadir("0.5deg_interpolation_factors.jld2"))

lon=[i for i=-179.:2:179., j=-89.:2:89.]
lat=[j for i=-179.:2:179., j=-89.:2:89.]
(f,i,j,w)=InterpolationFactors(Γ,vec(lon),vec(lat))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)
region = "NPAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)

mskC, mskW, mskS = get_msk(Γ)

adjust_exps = Dict()
interp_mask = Interpolate(PAC_msk,f,i,j,w) #interpolate using half-degree resolution
interp_mask = reshape(interp_mask, size(lon));
interp_mask[interp_mask .> 0] .= 1
interp_mask[interp_mask .== 0] .= NaN
adjust_exps["RegularMask"] = interp_mask
adjust_exps["lat"] = lat
adjust_exps["lon"] = lon

diagpath["mean_tau_noadjusts"] = vastdiagdir("seasonalclimatology_iter0", "run")
diagpath["mean_tau_adjusts"] = vastdiagdir("seasonalclimatology", "run_only_clim_tau_adjusts")

for expname in ["only_wind", "mean_tau_adjusts", "iter0_bulkformula", "mean_tau_noadjusts"]
    println(expname)
    filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
    datafilelist_τ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    nt = 312
    W_reg = zeros(size(lon)..., nt)
    for tt in 1:nt
        U, V, W, Ub, Vb, Wb = extract_eulerian_and_bolus_velocities(diagpath, expname, 
        datafilelist_τ[tt], γ, Γ, mskC, mskW, mskS)

        W_interp = Interpolate(W[:, 40],f,i,j,w) #interpolate using half-degree resolution
        W_interp = reshape(W_interp, size(lon));
        W_interp[W_interp .== 0.0] .= NaN
        W_reg[:, :, tt] .= W_interp
    end
    adjust_exps[expname] = W_reg
end

jldsave(datadir(region * "W_interpolated_lev40.jld2"), adjust_exps= adjust_exps)
