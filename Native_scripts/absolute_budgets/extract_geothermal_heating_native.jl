#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, Plots
using .OHC_helper

include(srcdir("config_exp.jl"))

#using . is faster 
(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());
 
ocean_mask = OHC_helper.wet_pts(Γ)
cell_depths = OHC_helper.get_cell_depths(ocean_mask, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));

inv_H = smush(cell_depths, γ); 
inv_cv = deepcopy(cell_volumes); 
for a in eachindex(H) 
    inv_H.f[a][inv_H.f[a].== 0.0] .= Inf
    inv_H.f[a] .= 1 ./ inv_H.f[a]
end
for a in eachindex(inv_cv) 
    inv_cv.f[a][inv_cv.f[a].== 0.0] .= Inf
    inv_cv.f[a] .= 1 ./ inv_cv.f[a]
end
tecco = 1992+1/24:1/12:2018

function filter_heat_budget(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid)

    filelist = searchdir(diagpath[expname],"trsp_3d_set2")
    datafilelist_H  = filter(x -> occursin("data",x),filelist) 
    filelist = searchdir(diagpath[expname],"trsp_3d_set3")
    datafilelist_R  = filter(x -> occursin("data",x),filelist)
    filelist = searchdir(diagpath[expname],"state_3d_set1") 
    datafilelist_θ  = filter(x -> occursin("data",x),filelist)
    filelist = searchdir(diagpath[expname], "state_2d_set1") # first filter for state_3d_set1
    datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
 
    κθ_conv3D = MeshArray(γ,Float32,50); uθ_conv3D = MeshArray(γ,Float32,50);
    κzθ_conv = MeshArray(γ,Float32,50); wθ_conv = MeshArray(γ,Float32,50);
    mean_fluxes = MeshArray(γ,Float32,50); fill!(mean_fluxes, 0.0)

    nt = 312;
    @time for tt = 1:nt
        println(tt)

        κθx, κθy, uθx, uθy = OHC_helper.extract_heatbudgetH(diagpath, expname, datafilelist_H[tt], γ)
        κzθ, wθ = OHC_helper.extract_heatbudgetR(diagpath, expname, datafilelist_R[tt], γ)

        OHC_helper.calc_UV_conv3D!(κθx, κθy, κθ_conv3D); 
        OHC_helper.calc_UV_conv3D!(uθx, uθy, uθ_conv3D);
        OHC_helper.calc_Wconv3D!(wθ, wθ_conv);
        OHC_helper.calc_Wconv3D!(κzθ, κzθ_conv);
        
        #take the average of the fluxes should be equal to (sθend - sθstart) / end
        for a in eachindex(mean_fluxes)
            tmp = κθ_conv3D.f[a] .+ uθ_conv3D.f[a] .+ wθ_conv.f[a] .+ κzθ_conv.f[a]
            tmp .= tmp .* inv_cv.f[a]
            mean_fluxes.f[a] .+= tmp ./ nt
        end
    end

    sθstart = OHC_helper.extract_sθ(expname,diagpath, γ, datafilelist_S[1], datafilelist_θ[1], inv_H)
    sθend = OHC_helper.extract_sθ(expname,diagpath, γ, datafilelist_S[nt], datafilelist_θ[nt], inv_H)
    dt = 3.154e+7 * (tecco[nt] -tecco[1])

    GTF = MeshArray(γ,Float32, 50)
    for a in eachindex(mean_fluxes)
        sθdiff = (sθend.f[a] .- sθstart.f[a])./ dt
        GTF.f[a]    .= sθdiff .- mean_fluxes.f[a]
    end
    
    return GTF
end

# for expname in keys(shortnames)
expname = "iter129_bulkformula"
GTF = filter_heat_budget(diagpath, expname, γ)

#plot GTF to check 
GTF2d = smush(GTF[:, 30:end], γ); 
import PyPlot as plt
proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
fig, ax = plt.subplots( figsize=(15,10), subplot_kw=Dict("projection"=> proj0))
ax.coastlines(resolution="110m")
ax.set_extent((-250, -65, -56, 60),crs=projPC)

for ff in 1:5
    ax.pcolormesh(λ[ff], ϕ[ff],  GTF2d[ff],
    vmin = 0, vmax = 1e-9, shading="nearest", transform=projPC, 
    rasterized = 10, cmap = "magma")               
end
fig

write(datadir("native/GeothermalHeating.data"),GTF)

# end
