include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson
using .OHC_helper
import PyPlot as plt 

include(srcdir("config_exp.jl"))
include("sections_helper.jl")

tecco = 1992+1/24:1/12:2018
nz = 50

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());
diagpath["clim_tau_iter0"] = "/vast/ECCOv4r4/exps/clim_tau_iter0/run/diags/"
 
ocean_mask = OHC_helper.wet_pts(Γ)
region = "NPAC"; 
abs_dist(x, r) = abs(x) < r

P02_ϕ = 30; P01_ϕ = 45;
P16_λ = -150

region = "NPAC"
PAC_mask = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region, extent = "full")

region2 = "PAC56"
PAC56_mask = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region2, extent = "full")
sections = Dict()
sections["P01"] = latitude_section_mask(ϕ, P01_ϕ, PAC56_mask)
sections["P02"] = latitude_section_mask(ϕ, P02_ϕ, PAC56_mask)
sections["P16"] = longitude_section_mask(λ, P16_λ, PAC56_mask)
sections["P16"] = sections["P16"] .* PAC_mask

proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

fig, ax = plt.subplots( figsize=(15,10), subplot_kw=Dict("projection"=> proj0))
ax.coastlines(resolution="110m")
gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
gl.top_labels = false; gl.right_labels = false
plot_basin_mask!(λ, ϕ, sections["P02"], projPC)    
fig

function filter_heat_budget(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid, mask)
    nz = 50
    cell_depths = OHC_helper.get_cell_depths(mask, ΔzF, Γ.hFacC); 
    cell_volumes = (OHC_helper.get_cell_volumes(area, cell_depths));
    ΔV = zeros(Float32, nz)
    [ΔV[k] = Float32(sum(cell_volumes[:, k])) for k=1:nz]
    filelist = searchdir(diagpath[expname],"state_3d_set1") 
    datafilelist_θ  = filter(x -> occursin("data",x),filelist)

    nt = length(datafilelist_θ); 
    θ_vars = zeros(Float32, nz, nt)
    S_vars = zeros(Float32, nz, nt)

    @time for tt = 1:nt
        println(tt)
        θname = datafilelist_θ[tt]
        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        Sθ = γ.read(diagpath[expname]*θname,MeshArray(γ,Float32,2*nz))
        θ = Sθ[:, 1:50]; S = Sθ[:, 51:end]
        θ_vars[:, tt]  = OHC_helper.sum_heat_flux_profile(θ .* cell_volumes, ΔV); 
        S_vars[:, tt]  = OHC_helper.sum_heat_flux_profile(S .* cell_volumes, ΔV); 
    end

    return θ_vars, S_vars, ΔV
end

expname = "iter129_bulkformula"
for (section, section_mask) in sections
    svename = datadir(region * "_" * expname * "_" * section * "_THETA_levels" * ".jld2")
    θ, S, ΔV = filter_heat_budget(diagpath, expname, γ, section_mask)
    jldsave(svename, θ = θ, S = S, ΔV = ΔV)
    println(svename)
end

SA = 
using Pkg; Pkg.add