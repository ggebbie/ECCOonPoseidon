
include("../../src/intro.jl")
include("GH19_helperfuncs.jl")

using Revise, DrWatson, Statistics,
NCDatasets, Printf, 
DataFrames, LaTeXStrings, Distances, JLD2
import PyPlot as plt

ds_EQ  = NCDataset("/home/ameza/GH19.jl/data/Theta_EQ-0015.nc")
ds_OPT = NCDataset("/home/ameza/GH19.jl/data/Theta_OPT-0015.nc")

year  = reverse(ds_OPT["year"][:]); nt = length(year)
lon = ds_OPT["longitude"][:]
lat = ds_OPT["latitude"][:]
depth = ds_OPT["depth"][:]; nz = length(depth)

theta_OPT = reverse(ds_OPT["theta"][:, :, :, :], dims =1) #reverse time for niceness
theta_EQ = reverse(ds_EQ["theta"][:, :, :, :], dims =1)

#check to see if I am plotting things correctly 
import PyCall
fig, ax = plt.subplots()
ax.contourf(lon, lat, theta_OPT[1, 1, :, :])
fig
wet_mask = (!isnan).(theta_OPT[1, :, :, :])

#coordinate meshgrid
LONS = lon' .* ones(length(lat))
LATS = lat .* ones(length(lon))'

#area of interest NPAC 
PAC_msk = 0.0 .* LATS
# if region == "NPAC"
#     PAC_msk .= (23 .<= LATS .<= 56) .&& (120 .<= LONS .<= 290)
#     PAC_msk .= Float32.(PAC_msk); PAC_msk[iszero.(PAC_msk)] .= NaN
#     cutout = (15 .<= LATS .<= 56) .&& (260 .<= LONS .<= 290)
#     PAC_msk .= PAC_msk .- cutout; PAC_msk .= abs.(PAC_msk)
#     cutout = (-90 .<= LATS .<= -20) .&& (120 .<= LONS .<= 150)
#     PAC_msk .= PAC_msk .- cutout; PAC_msk .= abs.(PAC_msk)
# elseif region == "PAC"
#area of interest PAC 

# PAC_msk = (-56 .<= LATS .<= 60) .&& (140 .<= LONS .<= 290)
# PAC_msk = Float32.(PAC_msk); PAC_msk[iszero.(PAC_msk)] .= NaN
# cutout = (15 .<= LATS .<= 56) .&& (260 .<= LONS .<= 290)
# PAC_msk = PAC_msk .- cutout; PAC_msk = abs.(PAC_msk)
# cutout = (-90 .<= LATS .<= -20) .&& (120 .<= LONS .<= 150)
# PAC_msk = PAC_msk .- cutout; PAC_msk = abs.(PAC_msk)

# PAC_msk[isnan.(PAC_msk)] .= 0.0

PAC_msk = (-56 .<= LATS .<= 60) .&& (140 .<= LONS .<= 260)

fig, ax = plt.subplots()
ax.contourf(lon, lat, theta_OPT[1, 1, :, :] .* PAC_msk)
fig
#volumes for each cell. volume = 0 if not in area of interst
volumes = GH19_cell_volumes(depth, lon, lat)
mask_volume = similar(volumes)
[mask_volume[k, :, :] .= volumes[k, :, :] .* PAC_msk .* wet_mask[k, :, :] for k = 1:nz]

fig, ax = plt.subplots()
ax.contourf(lon, lat,  mask_volume[10, :, :])
fig


#weight the data
ΔTs = []; ps = []
WOCE_times = findall(1872 .< year .< 1876)[1]
Challenger_times = findall(1989 .< year .< 2001)[end]
println(year[WOCE_times] - year[Challenger_times] )
data_labels = ["EQ-0015", "OPT-0015"]
E,F = trend_matrices(year[WOCE_times:Challenger_times])

for (i, data) in enumerate([theta_EQ, theta_OPT])
    weighted_temp = zeros(nt, nz)
    β = zeros(nz)

    #fill NaNs
    filled_data = copy(data)
    filled_data[isnan.(filled_data)] .= 0.0

    #volume weighted average 
    for tt in 1:nt, k in 1:nz
        weighted_temp[tt, k] =  sum(filled_data[tt, k, :, :] .* mask_volume[k, :, :]) / sum(mask_volume[k, :, :])
    end
    for (i, tt) in enumerate(WOCE_times:Challenger_times), k in 1:nz
        β[k] += F[2,i] * weighted_temp[tt, k] 
    end
    # p = contourf(year, depth, weighted_temp_anom, yflip = true, 
    #         levels = -0.5:.1:0.5, c = :balance, title = data_labels[i], 
    #         xlabel = "years", ylabel = "depth [m]")
    #simple temperature difference

    push!(ΔTs, deepcopy(β))
    # push!(ps, p)
end

# plot(ps...)
# plot(ΔTs[1], depth, yflip = true, label = data_labels[1], color = "blue", linestyle = :dash)
# plot!(ΔTs[2], depth, yflip = true, label = data_labels[2], 
#       color = "blue", linestyle = :solid, xlabel = "ΔT[K]", 
#       ylabel = "Depth [m]",  xticks = (-0.4:0.1:0.5), xlims = (-0.35, 0.5))
# vline!([0], c = "grey", label = nothing, 
# title = " Vertical profiles of temperature change \n between WOCE and Challenger",
# legend=:bottomright)
# close("all")
fig, ax = plt.subplots()
ax.plot(100 * ΔTs[2], depth); 
ax.set_xlabel(L"^\circ" * "C per century")
ax.set_ylabel("Depth [km]")

# ax.tick_params(bottom=true, left=true)
# ax.set_xlim(-0.15, 0.15)
ax.invert_yaxis()
fig

jldsave(datadir("OPT-0015_GH19_PAC.jld2"); ΔT_GH19 = ΔTs[2], depth_GH19 = depth)

