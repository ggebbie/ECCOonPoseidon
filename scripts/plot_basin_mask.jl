
using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools,
    PyPlot, JLD2, DrWatson, Statistics, JLD2
import CairoMakie as Mkie
import GeoMakie
using MAT

pth = MeshArrays.GRID_LLC90
γ = GridSpec("LatLonCap",pth)
Γ = GridLoad(γ;option="full")
basins=read(joinpath(pth,"v4_basin.bin"),MeshArray(γ,Float32))

basin_list=["Pacific","Atlantic","indian","Arctic","Bering Sea",
            "South China Sea","Gulf of Mexico","Okhotsk Sea",
            "Hudson Bay","Mediterranean Sea","Java Sea","North Sea",
            "Japan Sea", "Timor Sea","East China Sea","Red Sea",
            "Gulf","Baffin Bay","GIN Seas","Barents Sea"];

basin_name="Pacfic"
# basin_name = "Atlantic"

plot_basin = true
output2file = true

# Γ.hFacC[:,1] can be used as an indicator for wet points
# (there might be a better way to do this)
ocean_mask = Γ.hFacC[:,1]
ocean_mask[findall(ocean_mask.>0.0)].=1.0

basinID=findall(basin_list.==basin_name)[1]
basin_mask=similar(basins)
for ff in 1:5
    basin_mask[ff] .= ocean_mask[ff].*(basins[ff].==basinID) #put this into mat file
end

basin_name_nospace = replace(basin_name, " " => "")
file = matopen("GRID_LLC90_"*basin_name_nospace*".mat","w")
for ff in 1:5
    write(file,basin_name_nospace*"_mask"*string(ff),basin_mask[ff])
end
close(file)

if plot_basin == true
    ocean_mask[findall(ocean_mask.==0.0)].=NaN #replace()
    fig = Mkie.Figure(resolution = (900,600), backgroundcolor = :grey95)
    ax = Mkie.Axis(fig[1,1],xlabel="longitude",ylabel="latitude",title=basin_name* " (shown in red)")

    mx=maximum(basins)

    #used to plot basins
    for ff in 1:length(Γ.RAC)
        col=ocean_mask[ff][:].*(basins[ff][:].==basinID)
        kk=findall(col.>0.0)
        !isempty(kk) ? Mkie.scatter!(ax,Γ.XC[ff][kk],Γ.YC[ff][kk],color=:red,markersize=2.0) : nothing
        kk=findall((col.==0.0).*(!isnan).(ocean_mask[ff][:]))
        !isempty(kk) ? Mkie.scatter!(ax,Γ.XC[ff][kk],Γ.YC[ff][kk],color=basins[ff][kk],
            colorrange=(0.0,mx),markersize=2.0,colormap=:lisbon) : nothing
#    end
#    Mkie.Colorbar(fig[1,2], colormap=:lisbon, colorrange=(0.0, mx), height = Mkie.Relative(0.65))
#    Mkie.save("GRID_LLC90_plotted_"*basin_name*".png", fig)
#    ocean_mask[findall(isnan.(ocean_mask))].=0.0
#end

if plot_basin == true
    #ocean_mask[findall(ocean_mask.==0.0)].=NaN #replace()
    fig = Mkie.Figure(resolution = (900,600), backgroundcolor = :grey95)
    ax = Mkie.Axis(fig[1,1],xlabel="longitude",ylabel="latitude",title=basin_name* " (shown in red)")

    mx=maximum(basins)

    #used to plot basins
    for ff in 1:length(Γ.RAC)
        col=(basins[ff][:].==basinID)
        kk=findall(col.>0.0)
        !isempty(kk) ? Mkie.scatter!(ax,Γ.XC[ff][kk],Γ.YC[ff][kk],color=:red,markersize=2.0) : nothing
        #kk=findall((col.==0.0).*(!isnan).(ocean_mask[ff][:]))
        !isempty(kk) ? Mkie.scatter!(ax,Γ.XC[ff][kk],Γ.YC[ff][kk],color=basins[ff][kk],
            colorrange=(0.0,mx),markersize=2.0,colormap=:lisbon) : nothing
    end
    Mkie.Colorbar(fig[1,2], colormap=:lisbon, colorrange=(0.0, mx), height = Mkie.Relative(0.65))
    Mkie.save("testGRID_LLC90_plotted_"*basin_name*".png", fig)
    #ocean_mask[findall(isnan.(ocean_mask))].=0.0
end