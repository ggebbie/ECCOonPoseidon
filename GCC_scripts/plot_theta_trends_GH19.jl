using Revise, DrWatson, Statistics,
NCDatasets, Printf, 
DataFrames, LaTeXStrings,
Plots, Distances
include("GH19_helperfuncs.jl")

ds_EQ  = NCDataset("/home/ameza/GH19.jl/data/Theta_EQ-0015.nc")
ds_OPT = NCDataset("/home/ameza/GH19.jl/data/Theta_OPT-0015.nc")

year  = reverse(ds_OPT["year"][:]); nt = length(year)
lon = ds_OPT["longitude"][:]
lat = ds_OPT["latitude"][:]
depth = ds_OPT["depth"][:]; nz = length(depth)

theta_OPT = reverse(ds_OPT["theta"][:, :, :, :], dims =1) #reverse time for niceness
theta_EQ = reverse(ds_EQ["theta"][:, :, :, :], dims =1)

#check to see if I am plotting things correctly 
contourf(lon, lat, theta_OPT[1, 1, :, :])

wet_mask = (!isnan).(theta_OPT[1, :, :, :])

#coordinate meshgrid
LONS = lon' .* ones(length(lat))
LATS = lat .* ones(length(lon))'

#area of interest
PAC_msk = (-56 .<= LATS .<= 60) .&& (140 .<= LONS .<= 260)

#volumes for each cell. volume = 0 if not in area of interst
volumes = cell_volumes(depth, lon, lat)
mask_volume = similar(volumes)
[mask_volume[k, :, :] .= volumes[k, :, :] .* PAC_msk .* wet_mask[k, :, :] for k = 1:nz]

# contourf(lon, lat, mask_volume[-1, :, :])


#weight the data
ΔTs = []; ps = []
WOCE_times = findall(1872 .< year .< 1876)
Challenger_times = findall(1989 .< year .< 2001)

data_labels = ["EQ-0015", "OPT-0015"]
for (i, data) in enumerate([theta_EQ, theta_OPT])
    weighted_temp = zeros(nt, nz)

    #fill NaNs
    filled_data = copy(data)
    filled_data[isnan.(filled_data)] .= 0.0

    #volume weighted average 
    for tt in 1:nt, k in 1:nz
        weighted_temp[tt, k] =  sum(filled_data[tt, k, :, :] .* mask_volume[k, :, :]) / sum(mask_volume[k, :, :])
    end

    weighted_temp_anom =  weighted_temp' .- mean( weighted_temp', dims = 2)
    p = contourf(year, depth, weighted_temp_anom, yflip = true, 
            levels = -0.5:.1:0.5, c = :balance, title = data_labels[i], 
            xlabel = "years", ylabel = "depth [m]")
    #simple temperature difference
    ΔT = mean(weighted_temp_anom[:, Challenger_times], dims = 2) - weighted_temp_anom[:, WOCE_times]

    push!(ΔTs, ΔT)
    push!(ps, p)
end

# plot(ps...)

# plot(ΔTs[1], depth, yflip = true, label = data_labels[1], color = "blue", linestyle = :dash)
# plot!(ΔTs[2], depth, yflip = true, label = data_labels[2], 
#       color = "blue", linestyle = :solid, xlabel = "ΔT[K]", 
#       ylabel = "Depth [m]",  xticks = (-0.4:0.1:0.5), xlims = (-0.35, 0.5))
# vline!([0], c = "grey", label = nothing, 
# title = " Vertical profiles of temperature change \n between WOCE and Challenger",
# legend=:bottomright)

using JLD2
# close("all")

jldsave(datadir("OPT-0015_GH19.jld2"); ΔT_GH19 = ΔTs[2], depth_GH19 = depth)

