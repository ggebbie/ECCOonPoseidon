include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, DSP, PyCall
import PyPlot as plt 
include(srcdir("config_exp.jl"))

tecco = 1992+1/24:1/12:2018
nz = 50

@pyimport cmocean.cm as cmo


(ϕ,λ) = latlonC(γ)
area = readarea(γ)
λ_wrap = wrap_λ(λ)
ocean_mask = wet_pts(Γ)

include(srcdir("plot_and_dir_config.jl"))

function get_transports(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid, Γ)

    mskC, mskW, mskS = get_msk(Γ)

    fileroot = "trsp_3d_set1"
    filelist = searchdir(diagpath[expname],fileroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    nt = length(datafilelist);
    w_mean = MeshArray(γ, Float32, 50); fill!(w_mean, 0.0)
    v_mean = MeshArray(γ, Float32, 50); fill!(v_mean, 0.0)
    u_mean = MeshArray(γ, Float32, 50); fill!(u_mean, 0.0)

    wb_mean = MeshArray(γ, Float32, 50); fill!(wb_mean, 0.0)
    vb_mean = MeshArray(γ, Float32, 50); fill!(vb_mean, 0.0)
    ub_mean = MeshArray(γ, Float32, 50); fill!(ub_mean, 0.0)

    @time for tt = 1:nt
        println(tt)
        fname = datafilelist[tt]
        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        u, v, w, Ub, Vb, Wb = extract_eulerian_and_bolus_velocities(diagpath, expname, 
                                fname, γ, Γ, mskC, mskW, mskS)

        wb_mean .+= Wb ./nt
        vb_mean .+= Vb ./nt
        ub_mean .+= Ub ./nt

        w_mean .+= w ./nt
        v_mean .+= v ./nt
        u_mean .+= u ./nt

    end
    
    Eul_Vel = Dict(); Bol_Vel = Dict()
    Eul_Vel["U"] = u_mean
    Eul_Vel["V"] = v_mean
    Eul_Vel["W"] = w_mean

    Bol_Vel["U"] = ub_mean
    Bol_Vel["V"] = vb_mean
    Bol_Vel["W"] = wb_mean

    return Eul_Vel, Bol_Vel
end



vars =  ["only_init", "only_kappa", "only_sfc", "iter129_bulkformula",  "iter0_bulkformula"]
@time for expname in vars
    Eul_Vel, Bol_Vel = get_transports(diagpath, expname, γ, Γ); 
    svename = datadir("native/"* expname * "_mean_Eul_Bol_Trsp.jld2")
    jldsave(svename; Eul_Vel = Eul_Vel, Bol_Vel = Bol_Vel)
end

