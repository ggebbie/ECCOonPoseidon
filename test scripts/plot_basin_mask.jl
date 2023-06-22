include("../src/intro.jl")
include("../src/OHC_helper.jl")
using Revise 
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools,
    PyPlot, JLD2, DrWatson, Statistics, JLD2
import CairoMakie as Mkie
using .OHC_helper

import GeoMakie
include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
area = readarea(γ)
# basin_name = "Atlantic"

plot_basin = true
output2file = true

# Γ.hFacC[:,1] can be used as an indicator for wet points 
# (there might be a better way to do this)
region = "NPAC"
wet = OHC_helper.wet_pts(Γ)

basin_mask = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; region, extent = "full")

abs_dist(x, r) = abs(x) < r
ϕ_mask_min = OHC_helper.get_min_lat(ϕ, PAC_msk); ϕ_mask = ϕ .> -Inf
[ϕ_mask.f[ff] .= abs_dist.(ϕ.f[ff] .- ϕ_mask_min, 0.1) .* PAC_msk.f[ff] for ff in 1:5]

basin_mask[findall(basin_mask.==0.0)].=NaN
fig = Mkie.Figure(resolution = (900,600), backgroundcolor = :grey95)
ax = Mkie.Axis(fig[1,1],xlabel="longitude",ylabel="latitude",title=region* " (shown in red)")

mx=maximum(basins)
function adjust_lons(x)
    temp = x
    temp[findall(x .< 0)] = x[findall(x .< 0)] .+ 360
    return temp
end
#used to plot basins
for ff in 1:length(Γ.RAC)
    col=basin_mask[ff][:]
    kk=findall(col.>0.0)
    !isempty(kk) ? Mkie.scatter!(ax,adjust_lons(Γ.XC[ff][kk]),Γ.YC[ff][kk],color=:red,markersize=2.0) : nothing
    kk=findall(map(isnan, basin_mask[ff][:]) .* Bool.(wet[ff][:]))
    !isempty(kk) && Mkie.scatter!(ax,adjust_lons(Γ.XC[ff][kk]),Γ.YC[ff][kk],color=:grey,markersize=2.0 )
end

#used to plot basins
for ff in 1:length(Γ.RAC)
    col=ϕ_mask[ff][:]
    kk=findall(col.>0.0)
    !isempty(kk) ? Mkie.scatter!(ax,adjust_lons(Γ.XC[ff][kk]),Γ.YC[ff][kk],color=:green,markersize=2.0) : nothing
    kk=findall(map(isnan, basin_mask[ff][:]) .* Bool.(wet[ff][:]))
    !isempty(kk) && Mkie.scatter!(ax,adjust_lons(Γ.XC[ff][kk]),Γ.YC[ff][kk],color=:grey,markersize=2.0 )
end

# Mkie.Colorbar(fig[1,2], height = Mkie.Relative(0.65))
Mkie.save(plotsdir("/GRID_LLC90_plotted_"*region*".png"), fig)

