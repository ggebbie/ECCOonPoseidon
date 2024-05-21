include("../../src/intro.jl")
include("GH19_helperfuncs.jl")

using Revise, DrWatson, Statistics,ECCOonPoseidon, 
NCDatasets, Printf, MeshArrays, MITgcmTools, 
DataFrames, LaTeXStrings, Distances, JLD2, PyCall, Interpolations
using TMI
import PyPlot as plt

using ECCOonPoseidon
ds_EQ  = NCDataset("/home/ameza/GH19.jl/data/Theta_EQ-0015.nc")
ds_OPT = NCDataset("/home/ameza/GH19.jl/data/Theta_OPT-0015.nc")

regions = NCDataset("/home/ameza/ECCOonPoseidon/data/regions_180x90.nc")


include(srcdir("config_exp.jl"))
using ECCOonPoseidon

(ϕ,λ) = ECCOonPoseidon.latlonC(γ)


region = "PAC"
PAC_msk = ECCOonPoseidon.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region)
PAC_msk_arr = convert2array(PAC_msk);
λ_arr = convert2array(λ);
λ_arr[λ_arr .< 0.0] .+= 360
ϕ_arr = convert2array(ϕ);

ϕ_arr = ϕ_arr[:, 1:270];  ϕ_arr = vcat(ϕ_arr[14:end, :], ϕ_arr[1:13, :])
λ_arr = λ_arr[:, 1:270];  λ_arr = vcat(λ_arr[14:end, :], λ_arr[1:13, :])
PAC_msk_arr = PAC_msk_arr[:, 1:270]
PAC_msk_arr = vcat(PAC_msk_arr[14:end, :], PAC_msk_arr[1:13, :])

year  = reverse(ds_OPT["year"][:]); nt = length(year)
lon = ds_OPT["longitude"][:]
lat = ds_OPT["latitude"][:]
depth = ds_OPT["depth"][:]; nz = length(depth)

theta_OPT = reverse(ds_OPT["theta"][:, :, :, :], dims =1) #reverse time for niceness
theta_EQ = reverse(ds_EQ["theta"][:, :, :, :], dims =1)

#coordinate meshgrid
LONS = lon' .* ones(length(lat))
LATS = lat .* ones(length(lon))'

#area of interest NPAC 


#check to see if I am plotting things correctly 


fig, ax = plt.subplots()
ax.contourf(lon, lat, theta_OPT[1, 1, :, :])
fig
wet_mask = (!isnan).(theta_OPT[1, :, :, :])

PAC_msk = regions["NPAC"][:, :]' .+  regions["TROPPAC"][:, :]'
PAC_msk[isnan.(PAC_msk)] .= 0

fig, ax = plt.subplots()
ax.contourf(lon, lat, theta_OPT[1, 1, :, :] .* PAC_msk)
fig

volumes = GH19_cell_volumes(depth, lon, lat)
mask_volume = similar(volumes)
[mask_volume[k, :, :] .= volumes[k, :, :] .* PAC_msk .* wet_mask[k, :, :] for k = 1:nz]

fig, ax = plt.subplots()
ax.contourf(lon, lat,  mask_volume[30, :, :])
fig


#weight the data
ΔTs = []; ps = []
WOCE_times = findall(1872 .< year .< 1876)[1]
# WOCE_times = findall(1989 .< year .< 2017)[1]

Challenger_times = findall(1989 .< year .< 1993)[end]
println(year[WOCE_times] - year[Challenger_times] )
data_labels = ["EQ-0015", "OPT-0015"]
E,F = ECCOonPoseidon.trend_matrices(year[WOCE_times:Challenger_times])

for (i, data) in enumerate([theta_EQ, theta_OPT])
    filled_data = copy(data)
    filled_data[isnan.(filled_data)] .= 0.0
    weighted_temp = zeros(2, nz)
    for (i, tt) in enumerate([WOCE_times, Challenger_times]), k in 1:nz
        weighted_temp[i, k] =  sum(filled_data[tt, k, :, :] .* mask_volume[k, :, :]) / sum(mask_volume[k, :, :])
    end

    push!(ΔTs, (weighted_temp[end, :] .- weighted_temp[1, :]) ./ (year[WOCE_times] - year[Challenger_times]))
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

jldsave(datadir("OPT-0015_GH19_PAC4.jld2"); ΔT_GH19 = ΔTs[2], depth_GH19 = depth)

