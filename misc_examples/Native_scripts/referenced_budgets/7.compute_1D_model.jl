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
cell_depths = OHC_helper.get_cell_depths(ocean_mask, ΔzF, Γ.hFacC); 
H = OHC_helper.smush(cell_depths, γ); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area .* PAC_msk, cell_depths));

runpath,diagpath = listexperiments(exprootdir());

ϕ_min_mask, ϕ_max_mask = OHC_helper.get_ϕ_max_min_mask(region, Γ, λ, ϕ, basins, basin_list)
mskC, mskW, mskS = OHC_helper.get_msk(Γ)

#define control volume 
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= z[:].<= uplvl)
cell_volumes = cell_volumes .* PAC_msk
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
ho = 1029; cp = 3994
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

    Wbar =  Wbar .* area32 #volume flux 

    sθbar = γ.read(datadir(expname * "/" * meanvarname("THETA_ref_mean")),MeshArray(γ,Float32,50))
    sθzprime = MeshArray(γ,Float32,50); 
    sθzbar = MeshArray(γ,Float32,50); 
    interpolate_to_vertical_faces!(sθbar, sθzbar)
    nt = 312;
    var_names = ["Wθ_Bottom", "W_Bottom", "Wθ_Top", "W_Top"]
    vars = Dict(varname => zeros(Float32, 1) for varname in var_names)
    Wbθb = Wbar .* sθzbar
    vars["Wθ_Bottom"] .= sum(Wbθb[:, lvls[end] + 1] .* PAC_msk )
    vars["W_Bottom"] .= sum(Wbar[:, lvls[end] + 1] .* PAC_msk ) 

    vars["Wθ_Top"] .= sum(Wbθb[:, lvls[1]] .* PAC_msk )
    vars["W_Top"] .= sum(Wbar[:, lvls[1]] .* PAC_msk )

    return vars
end

experiments = ["iter129_bulkformula", "iter0_bulkformula", "seasonalclimatology", "seasonalclimatology_iter0"]
for expname in experiments
    vars = filter_heat_budget(diagpath, expname, γ)
    savename = datadir("native/" * expname * region * "_THETA_vertical_fluxes_mean_" * suffix * ".jld2")
    jldsave(savename, dθ = vars)
end

using DataFrames

ne = length(experiments)
df = DataFrame(Wθ_bot=zeros(ne), W_bot=zeros(ne), θ_bot=zeros(ne), 
Wθ_top=zeros(ne), W_top=zeros(ne), θ_top=zeros(ne), names = experiments)
df
for (i, expname) in enumerate(experiments)
    fname = datadir("native/" * expname * region * "_THETA_vertical_fluxes_mean_" * suffix * ".jld2")
    v = load(fname)["dθ"]
    df.Wθ_bot[i] = v["Wθ_Bottom"][1] * 1e-6
    df.W_bot[i] = v["W_Bottom"][1] * 1e-6
    df.θ_bot[i] = v["Wθ_Bottom"][1] / (v["W_Bottom"][1])
    df.Wθ_top[i] = -v["Wθ_Top"][1] * 1e-6
    df.W_top[i] = -v["W_Top"][1] * 1e-6
    df.θ_top[i] = v["Wθ_Top"][1] /(v["W_Top"][1])
end

df