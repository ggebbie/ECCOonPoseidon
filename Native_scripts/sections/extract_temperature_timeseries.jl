include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson
import PyPlot as plt 

include(srcdir("config_exp.jl"))
include("sections_helper.jl")

tecco = 1992+1/24:1/12:2018
nz = 50

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());
 
ocean_mask = wet_pts(Γ)
abs_dist(x, r) = abs(x) < r

region2 = "PAC56"
PAC56_mask = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region2, extent = "full")
sections = Dict()

sections["P01"] = latitude_section_mask(ϕ, 45.6, PAC56_mask)
sections["P02"] = latitude_section_mask(ϕ, 30.3, PAC56_mask)
sections["P03"] = latitude_section_mask(ϕ, 25.3, PAC56_mask)
sections["P06"] = latitude_section_mask(ϕ, -31.9, PAC56_mask)
sections["P21"] = latitude_section_mask(ϕ, -17.8, PAC56_mask)

sections["P09"] = longitude_section_mask(λ, 138.2, PAC56_mask)
sections["P10"] = longitude_section_mask(λ, 147.4, PAC56_mask)
sections["P14"] = longitude_section_mask(λ, 178.3, PAC56_mask)
sections["P16"] = longitude_section_mask(λ, -151.1, PAC56_mask)
sections["P17"] = longitude_section_mask(λ, -139.9, PAC56_mask)
sections["P18"] = longitude_section_mask(λ, -106.2, PAC56_mask)

function filter_heat_budget(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid, mask)
    nz = 50
    cell_depths = get_cell_thickness(mask, ΔzF, Γ.hFacC); 
    cell_volumes = get_cell_volumes(area, cell_depths);
    ΔV = lateral_sum(cell_volumes)

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

        θ_vars[:, tt] .= lateral_sum(θ .* cell_volumes); 
        θ_vars[:, tt] .= θ_vars[:, tt] ./ ΔV

        S_vars[:, tt] .= lateral_sum(S .* cell_volumes); 
        S_vars[:, tt] .= S_vars[:, tt] ./ ΔV
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

proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

fig, ax = plt.subplots( figsize=(15,10), subplot_kw=Dict("projection"=> proj0))
ax.coastlines(resolution="110m")
gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
gl.top_labels = false; gl.right_labels = false
plot_basin_mask!(λ, ϕ, sections["P02"], projPC)    

fig