#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../../src/intro.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
Printf, LaTeXStrings, PyCall, GibbsSeaWater, Dierckx, Interpolations

import PyPlot as plt, NaNMath as nm
@pyimport seaborn as sns
@pyimport pandas as pd
@pyimport cmocean.cm as cmo

include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());

include(srcdir("plot_and_dir_config.jl"))
ocean_mask = wet_pts(Γ)
cell_depths = get_cell_thickness(ocean_mask, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths);


tecco= 1992+1/24:1/12:2018 # ecco years
i100m = Base.findmin(abs.(z .- 100))[2]
i300m = Base.findmin(abs.(z .- 300))[2]

index_means = Dict() 

vertical_average(x) = vertical_sum(x) ./ 2
z_deriv(x, dz) = (x[:, 2] .- x[:, 1]) ./ dz
expname = "iter129_bulkformula"
for iidx in [i100m, i300m]

    dz = z[iidx + 1] .- z[iidx]
    vol_100m = cell_volumes[:, iidx:iidx + 1]
    ps = pstdz[iidx:iidx+1]
    filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
    datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    nt = 312; nz = 50
    index_mean = MeshArray(γ,Float32); fill!(index_mean, 0.0)

    for tt = 1:nt
        θname = datafilelist_θ[tt]
        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        @time θSz = γ.read(diagpath[expname]*θname,MeshArray(γ,Float32,2*nz))
        θtmp = θSz[:, iidx:iidx + 1]; Stmp = θSz[:, nz+iidx:nz+iidx+1]

        for a in eachindex(θtmp)
            θtmp.f[a][iszero.(vol_100m.f[a])] .= NaN
            Stmp.f[a][iszero.(vol_100m.f[a])] .= NaN
        end

        CTtmp = gsw_ct_from_pt.(Stmp,θtmp)
        SAtmp = MeshArray(γ,Float32, 2)
        for a in eachindex(SAtmp)
            SAtmp.f[a] .= gsw_sa_from_sp.(Stmp.f[a],ps[a[2]],
                                            λ.f[a[1]],ϕ.f[a[1]])
        end

        θavg = vertical_average(CTtmp); Savg = vertical_average(SAtmp)

        α = gsw_alpha.(Savg,θavg,mean(pstdz[iidx:iidx+1]))
        β = gsw_beta.(Savg,θavg,mean(pstdz[iidx:iidx+1]))

        # N² = -9.81 .* z_deriv(rhotmp) ./ 1043
        αN² = 9.81 .* α .* z_deriv(CTtmp, dz)
        βN² = -9.81 .* β .* z_deriv(SAtmp, dz)
        N² = αN² .+ βN²
        index = (αN² .- βN²) ./ N²
        index_mean .+= index ./ nt
    end
    index_means[round(abs.(z[iidx]))] = index_mean
end

#plotting 100 meters
proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

#plotting climatology :-)
fig, ax = plt.subplots( figsize=(15,10), subplot_kw=Dict("projection"=> proj0))
ax.coastlines(resolution="110m")
gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
gl.top_labels = false; gl.right_labels = false
cms = []
for ff in 1:5
    c = ax.pcolormesh(λ[ff], ϕ[ff],  index_means[95.0][ff],
    vmin = -1, vmax = 1, shading="nearest", transform=projPC, 
    rasterized = true, cmap = cmo.balance)          
    push!(cms, c)     
end
ax.set_title("z = 100 meters, t = 1992 - 2017")
fig.colorbar(cms[1], ax = ax, orientation = "horizontal", fraction = 0.045, label = L"\kappa^2 / N^2")
fig



#plotting 100 meters
proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

#plotting climatology :-)
fig, ax = plt.subplots( figsize=(15,10), subplot_kw=Dict("projection"=> proj0))
ax.coastlines(resolution="110m")
gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
gl.top_labels = false; gl.right_labels = false
cms = []
for ff in 1:5
    c = ax.pcolormesh(λ[ff], ϕ[ff],  index_means[300.0][ff],
    vmin = -1, vmax = 1, shading="nearest", transform=projPC, 
    rasterized = true, cmap = cmo.balance)          
    push!(cms, c)     
end
ax.set_title("z = 300 meters, t = 1992 - 2017")
fig.colorbar(cms[1], ax = ax, orientation = "horizontal", fraction = 0.045, label = L"\kappa^2 / N^2")
fig
