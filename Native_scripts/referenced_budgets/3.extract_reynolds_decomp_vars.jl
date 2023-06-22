#this script references θ to some control volume and 
#also computes the mean fields of U, V, W, θ
#where U, V, W are the eulerian + bolus velocities

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
using .OHC_helper
import NaNMath as nm
include(srcdir("config_exp.jl"))


(ϕ,λ) = latlonC(γ)
area = readarea(γ)
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

ocean_mask = OHC_helper.wet_pts(Γ)
cell_depths = OHC_helper.get_cell_depths(ocean_mask, ΔzF, Γ.hFacC); 
H = OHC_helper.smush(cell_depths, γ); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H

region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= z[:].<= uplvl)

runpath,diagpath = listexperiments(exprootdir());

mskC, mskW, mskS = OHC_helper.get_msk(Γ)

meanvarname(x) = x * "_" * region * "_" * suffix *".data"
varname(x, timestamp) = x * "_" * region * "_" * suffix * "_" * timestamp *".data"

function filter_advection_budget(diagpath::Dict{String, String}, expname::String, γ::gcmgrid, 
    mskC::MeshArrays.gcmarray{T,2,Matrix{T}}, 
    mskW::MeshArrays.gcmarray{T,2,Matrix{T}}, 
    mskS::MeshArrays.gcmarray{T,2,Matrix{T}}) where T<:AbstractFloat

    mkpath(datadir(expname))
    filelist = searchdir(diagpath[expname], "state_2d_set1") # first filter for state_3d_set1
    datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
    datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
    datafilelist_uvw  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    sθbar = MeshArray(γ,Float32,50); fill!(sθbar, 0.0)
    Ubar = MeshArray(γ,Float32,50); fill!(Ubar, 0.0)
    Vbar = MeshArray(γ,Float32,50); fill!(Vbar, 0.0)
    Wbar = MeshArray(γ,Float32,50); fill!(Wbar, 0.0)
    nt = 312
    ctrl_vol = sum(cell_volumes[:, lvls])
    for tt = 1:nt
        println(tt)
        sθ = OHC_helper.extract_sθ(expname,diagpath, γ, datafilelist_S[tt], datafilelist_θ[tt], inv_H)
        sθ_ref = sum(sθ[:, lvls] .* cell_volumes[:, lvls]) / ctrl_vol
        sθ = sθ .- sθ_ref #reference the defined control volume
        u, v, w, Ub, Vb, Wb = OHC_helper.extract_velocities_and_bolus(diagpath, expname , 
        datafilelist_uvw[tt], γ, Γ, mskC, mskW, mskS)
        timestamp = datafilelist_θ[tt][15:end-5]
        write(datadir(expname * "/" * varname("THETA_ref", timestamp)),sθ)

        #take the time average of velocity and referenced faces
        for a in eachindex(sθbar)
            sθbar.f[a] .+= (sθ.f[a]) ./ nt
            Ubar.f[a] .+= (u.f[a] .+ Ub.f[a]) ./ nt
            Vbar.f[a] .+= (v.f[a] .+ Vb.f[a]) ./ nt
            Wbar.f[a] .+= (w.f[a] .+ Wb.f[a]) ./ nt
        end
    end

    write(datadir(expname * "/" * meanvarname("U_EB_mean")),Ubar)
    write(datadir(expname * "/" * meanvarname("V_EB_mean")),Vbar)
    write(datadir(expname * "/" * meanvarname("W_EB_mean")),Wbar)
    write(datadir(expname * "/" * meanvarname("THETA_ref_mean")),sθbar)

    println("saving files ")
end 

vars = ["seasonalclimatology", "seasonalclimatology_iter0"]

for expname in vars
    filter_advection_budget(diagpath, expname, γ, mskC, mskW, mskS)
end