using Pkg
Pkg.activate(".")
include("../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings
import NaNMath as nm
using PyCall

include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

include(srcdir("plot_and_dir_config.jl"))
 
ocean_mask = wet_pts(Γ)
region = "PAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "false", include_bering = true)

tecco = 1992+1/24:1/12:2018

vars =  ["only_init", "only_kappa", "only_sfc", "iter129_bulkformula",  "iter0_bulkformula"]
vars =  ["only_buoyancy", "only_wind"]

uplvl = -2.4e3; botlvl = -2.6e3; suffix = "2to3"
lvls = findall( botlvl .<= -z[:].<= uplvl)

diagpath["mean_tau_noadjust_redo_bf"] = vastdiagdir("seasonalclimatology", "run_only_clim_iter0_tau_bf")
diagpath["mean_tau_yesadjust_redo_bf"] = vastdiagdir("seasonalclimatology", "run_only_clim_iter129_tau_bf2")

@time for expname in ["mean_tau_noadjust_redo_bf", "mean_tau_yesadjust_redo_bf"]
    println(expname)
    w_timeseries = MeshArray(γ,Float32,312)
    w_eul_timeseries = MeshArray(γ,Float32,312)
    w_bol_timeseries = MeshArray(γ,Float32,312)

    filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
    datafilelist_uvw  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    mskC, mskW, mskS = get_msk(Γ)
    for tt = 1:312
        println(tt)
        u, v, w, Ub, Vb, Wb = extract_eulerian_and_bolus_velocities(diagpath, expname, 
        datafilelist_uvw[tt], γ, Γ, mskC, mskW, mskS)
        w_eul_timeseries.f[:, tt] .= w.f[:, lvls[1]]
        w_bol_timeseries.f[:, tt] .= Wb.f[:, lvls[1]]

        W = w .+ Wb
        w_timeseries.f[:, tt] .= W.f[:, lvls[1]]
    end
    savename = datadir("native/" * expname * "_W_residual_2500m.jld2")
    jldsave(savename, Wres = w_timeseries)
    savename = datadir("native/" * expname * "_W_eul_2500m.jld2")
    jldsave(savename, Weul = w_eul_timeseries)
    savename = datadir("native/" * expname * "_W_bol_2500m.jld2")
    jldsave(savename, Wbol = w_bol_timeseries)

end