#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
using .OHC_helper
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

ocean_mask = OHC_helper.wet_pts(Γ)
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
H = OHC_helper.smush(cell_depths, γ); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));

runpath,diagpath = listexperiments(exprootdir());

ϕ_min_mask, ϕ_max_mask = OHC_helper.get_ϕ_max_min_mask(region, Γ, λ, ϕ, basins, basin_list)
mskC, mskW, mskS = OHC_helper.get_msk(Γ)

#define control volume 
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= z[:].<= uplvl)
ctrl_vol = sum(cell_volumes[:, lvls])

meanvarname(x) = x * "_" * region * "_" * suffix *".data"
varname(x, timestamp) = x * "_" * region * "_" * suffix * "_" * timestamp *".data"
function reynolds_products(xbar::MeshArrays.gcmarray{T, 2, Matrix{T}}, xprime::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    ybar::MeshArrays.gcmarray{T, 2, Matrix{T}}, yprime::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:Real
    xbyb = T.(similar(xbar)); xpyp = T.(similar(xbar))
    xbyp = T.(similar(xbar)); xpyb = T.(similar(xbar))
    for a in eachindex(xbyb)
        xbyb.f[a] .= xbar.f[a] .* ybar.f[a]
        xpyp.f[a] .= xprime.f[a] .* yprime.f[a]
        xbyp.f[a] .= xbar.f[a] .* yprime.f[a]
        xpyb.f[a] .= xprime.f[a] .* ybar.f[a]
    end
    return xbyb, xbyp, xpyb, xpyp
end

function interpolate_to_vertical_faces!(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    ma_interp::MeshArrays.gcmarray{T, 2, Matrix{T}})  where T<:Real
    for k = 1:49 
        ma_interp[:, k+1] = (ma[:, k] .+ ma[:, k+1]) ./2 #linear interpolate to faces
    end
end

area32 = Float32.(area)

function filter_heat_budget(diagpath::Dict{String, String}, expname::String, γ::gcmgrid)

    filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
    datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
    datafilelist_uvw  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"trsp_3d_set2")
    datafilelist_H  = filter(x -> occursin("data",x),filelist) 
    filelist = searchdir(diagpath[expname],"trsp_3d_set3")
    datafilelist_R  = filter(x -> occursin("data",x),filelist)

    Ubar = γ.read(datadir(expname * "/" * meanvarname("U_EB_mean")),MeshArray(γ,Float32,50))
    Vbar = γ.read(datadir(expname * "/" * meanvarname("V_EB_mean")),MeshArray(γ,Float32,50))
    Wbar = γ.read(datadir(expname * "/" * meanvarname("W_EB_mean")),MeshArray(γ,Float32,50))

    Wbar =  Wbar .* area32

    sθbar = γ.read(datadir(expname * "/" * meanvarname("THETA_ref_mean")),MeshArray(γ,Float32,50))
    sθzprime = MeshArray(γ,Float32,50); 
    sθzbar = MeshArray(γ,Float32,50); 
    interpolate_to_vertical_faces!(sθbar, sθzbar)
    nt = 312;
    var_names = ["Wbθb_Bottom", "Wbθp_Bottom", "Wpθb_Bottom", "Wpθp_Bottom", 
                "Wbθb_Top", "Wbθp_Top", "Wpθb_Top", "Wpθp_Top"]
    vars = Dict(varname => zeros(Float32, nt) for varname in var_names)

    @time for tt = 1:nt
        println(tt)
        timestamp = datafilelist_θ[tt][15:end-5]

        sθ = γ.read(datadir(expname * "/" * varname("THETA_ref", timestamp)),MeshArray(γ,Float64,50))
        sθ = Float32.(sθ)
        sθprime = sθ .- sθbar
        fill!(sθzprime, 0.0)
        interpolate_to_vertical_faces!(sθprime, sθzprime)
        u, v, w, Ub, Vb, Wb = OHC_helper.extract_velocities_and_bolus(diagpath, expname, datafilelist_uvw[tt], 
                                                                      γ, Γ, mskC, mskW, mskS)
        for i in eachindex(u)
          u.f[i] .= u.f[i] .+ Ub.f[i]
          v.f[i] .= v.f[i] .+ Vb.f[i]
          w.f[i] .= w.f[i] .+ Wb.f[i]
        end  

        # Utrsp, Vtrsp = OHC_helper.UVtoTrsp(u, v, Γ) not needed for now
        Wprime = w .* area32; Wprime = Wprime .- Wbar; 
        Wbθb, Wbθp, Wpθb, Wpθp = reynolds_products(Wbar, Wprime, sθbar, sθzprime)

        vars["Wbθb_Bottom"][tt] = sum(Wbθb[:, lvls[end] + 1] .* PAC_msk )/ ctrl_vol
        vars["Wbθp_Bottom"][tt] = sum(Wbθp[:, lvls[end] + 1] .* PAC_msk )/ ctrl_vol
        vars["Wpθb_Bottom"][tt] = sum(Wpθb[:, lvls[end] + 1] .* PAC_msk )/ ctrl_vol
        vars["Wpθp_Bottom"][tt] = sum(Wpθp[:, lvls[end] + 1] .* PAC_msk )/ ctrl_vol

        vars["Wbθb_Top"][tt] = sum(Wbθb[:, lvls[1]] .* PAC_msk )/ ctrl_vol
        vars["Wbθp_Top"][tt] = sum(Wbθp[:, lvls[1]] .* PAC_msk )/ ctrl_vol
        vars["Wpθb_Top"][tt] = sum(Wpθb[:, lvls[1]] .* PAC_msk )/ ctrl_vol
        vars["Wpθp_Top"][tt] = sum(Wpθp[:, lvls[1]] .* PAC_msk )/ ctrl_vol

    end
    return vars
end

vars = ["iter129_bulkformula", "iter0_bulkformula"]
for expname in vars
    vars = filter_heat_budget(diagpath, expname, γ)
    savename = datadir("native/" * expname * region * "_THETA_budget_ref_Reynolds_" * suffix * ".jld2")
    jldsave(savename, dθ = vars)
end
