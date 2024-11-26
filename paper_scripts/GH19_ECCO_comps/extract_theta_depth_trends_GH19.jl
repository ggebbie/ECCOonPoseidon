include("../../src/intro.jl")
include("GH19_helperfuncs.jl")

using Revise, DrWatson, Statistics,ECCOonPoseidon, 
NCDatasets, Printf, MeshArrays, MITgcmTools, 
DataFrames, LaTeXStrings, Distances, JLD2, PyCall, Interpolations
import PyPlot as plt

using ECCOonPoseidon
ds_EQ  = NCDataset("/home/ameza/GH19.jl/data/Theta_EQ-0015.nc")
# ds_OPT = NCDataset("/home/ameza/GH19.jl/data/Theta_OPT-0015.nc")

ds_OPT_anom = NCDataset("/home/ameza/GH19.jl/data/Theta_anom_OPT-0015.nc")
ds_OPT  = NCDataset("/home/ameza/GH19.jl/data/Theta_OPT-0015.nc")

regions = NCDataset("/home/ameza/ECCOonPoseidon/data/regions_180x90.nc")


include(srcdir("config_exp.jl"))
using ECCOonPoseidon

(ϕ,λ) = ECCOonPoseidon.latlonC(γ)


region = "PAC"

year  = reverse(ds_OPT["year"][:]); nt = length(year)
lon = ds_OPT["longitude"][:]
lat = ds_OPT["latitude"][:]
depth = ds_OPT["depth"][:]; 
nz = length(depth)

theta_OPT = reverse(ds_OPT["theta"][:, :, :, :], dims =1) #reverse time for niceness
theta_EQ = reverse(ds_EQ["theta"][:, :, :, :], dims =1)

#coordinate meshgrid
LONS = lon' .* ones(length(lat))
LATS = lat .* ones(length(lon))'

fig, ax = plt.subplots()
ax.contourf(lon, lat, theta_EQ[1, 1, :, :])
fig

PAC_msk = regions["NPAC"][:, :]' .+  regions["TROPPAC"][:, :]' .+ regions["SUBANTPAC"][:, :]'
PAC_msk = PAC_msk .* (1 .* (LATS .> -45))
AUS = 1. .* (LATS .< -30) .* (LONS .< 145); AUS[AUS .== 1.0] .= 1 * NaN
CH =  1. .* (LATS .> 0) .* (LONS .> 105) .* (LONS .< 200); 
UP = 1. .* (LATS .< 65); UP[UP .== 0.0] .= 1 * NaN
PAC_msk  = (PAC_msk .+ AUS .+ CH) .* UP
PAC_msk[isnan.(PAC_msk)] .= 0
PAC_msk[PAC_msk  .>= 1.0] .= 1

fig, ax = plt.subplots()
ax.contourf(lon, lat, theta_OPT[1, 1, :, :] .* PAC_msk)
fig

volumes = GH19_cell_volumes(depth, lon, lat)
mask_volume = 0 .* similar(volumes)
wet_mask = (!isnan).(theta_OPT[1, :, :, :])
[mask_volume[k, :, :] .= (volumes[k, :, :] .* PAC_msk .* wet_mask[k, :, :]) for k = 1:nz]
mask_volume = Float32.(mask_volume)
fig, ax = plt.subplots()
ax.contourf(lon, lat,  volumes[30, :, :])
fig


#weight the data
ΔTs = []; ps = []
Challenger_times= findall(1872 .< year .< 1876)[end]
WOCE_times = findall(1989 .< year .< 2017)[1]
findall(1989 .< year .< 1993)

println(year[WOCE_times] - year[Challenger_times] )
data_labels = ["EQ-0015", "OPT-0015"]
# E,F = ECCOonPoseidon.trend_matrices(year[Challenger_times:WOCE_times])

for (i, data) in enumerate([theta_EQ, theta_OPT])
    filled_data =1 .* data
    filled_data[isnan.(filled_data)] .= 0.0
    print(minimum(filled_data))

    weighted_temp = zeros(2, nz)
    for (i, tt) in enumerate([WOCE_times, Challenger_times]), k in 1:nz
        weighted_temp[i, k] =  sum(filled_data[tt, k, :, :] .* mask_volume[k, :, :]) / sum(mask_volume[k, :, :])
    end

    push!(ΔTs, (weighted_temp[1, :] .- weighted_temp[end, :]) ./ (year[WOCE_times] - year[Challenger_times]))
end

fig, ax = plt.subplots()
ax.plot(ΔTs[1][20:end] .* 115, depth[20:end]); 
ax.plot(ΔTs[2][20:end].* 115, depth[20:end], label = "Anom"); 
ax.set_xlabel(L"^\circ" * "C per century")
ax.set_ylabel("Depth [km]")
ax.legend()
# ax.tick_params(bottom=true, left=true)
# ax.set_xlim(-0.15, 0.15)
ax.invert_yaxis()
fig

jldsave(datadir("OPT-0015_GH19_PAC_FINAL.jld2"); ΔT_GH19 = ΔTs[2], depth_GH19 = depth)

