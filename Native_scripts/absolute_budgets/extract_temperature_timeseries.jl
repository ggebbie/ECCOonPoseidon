include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson
using .OHC_helper
import PyPlot as plt 

include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());
diagpath["clim_tau_iter0"] = "/vast/ECCOv4r4/exps/clim_tau_iter0/run/diags/"
 
ocean_mask = OHC_helper.wet_pts(Γ)
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = (OHC_helper.get_cell_volumes(area, cell_depths));

H = OHC_helper.smush(cell_depths, γ); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H
tecco = 1992+1/24:1/12:2018
nz = 50

ΔV = zeros(Float32, nz)
[ΔV[k] = Float32(sum(cell_volumes[:, k])) for k=1:nz]
H = OHC_helper.smush(cell_depths, γ)

function filter_heat_budget(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid)

    filelist = searchdir(diagpath[expname],"state_3d_set1") 
    datafilelist_θ  = filter(x -> occursin("data",x),filelist)
    filelist = searchdir(diagpath[expname], "state_2d_set1") # first filter for state_3d_set1
    datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
 
    nt = length(datafilelist_θ); nz = 50
    vars = zeros(Float32, nz, nt)

    @time for tt = 1:nt
        println(tt)
        fnameθ = datafilelist_θ[tt]
        fnameS = datafilelist_S[tt]

        sθ = OHC_helper.extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_H)
        vars[:, tt]  = OHC_helper.sum_heat_flux_profile(sθ .* cell_volumes, ΔV); 
    end

    return vars
end

vars =  ["climatological_tau", "nosfcadjust"]

for expname in vars
    θ = filter_heat_budget(diagpath, expname, γ)
    svename = datadir(region * "_" * expname * "_THETA_levels" * ".jld2")
    jldsave(svename, θ = θ)
end


labels = ["Full Fluxes", "Climatological Fluxes", "Climatological Wind Stress"]
sumV = sum(ΔV[lvls])
mean_state(x) = sum(x[lvls, :] .* ΔV[lvls], dims = 1) / sumV
fig, ax = plt.subplots()
for (i, expname) in enumerate(keys(θ))
    ax.plot(tecco, mean_state(θ[expname])[:], label = expname)
end
ax.legend()
fig