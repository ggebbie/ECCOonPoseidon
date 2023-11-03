include("../../src/intro.jl")
include("GH19_helperfuncs.jl")

using Revise, DrWatson, Statistics,ECCOonPoseidon, 
NCDatasets, Printf, MeshArrays, MITgcmTools, 
DataFrames, LaTeXStrings, Distances, JLD2, PyCall, Interpolations
import PyPlot as plt

ds_EQ  = NCDataset("/home/ameza/GH19.jl/data/Theta_EQ-0015.nc")
ds_OPT = NCDataset("/home/ameza/GH19.jl/data/Theta_OPT-0015.nc")

include(srcdir("config_exp.jl"))
(ϕ,λ) = ECCOonPoseidon.latlonC(γ)


region = "PAC"
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; 
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

PAC_msk = (-38 .<= LATS .<= 64) .&& (115 .<= LONS .<= 300)
not_sel2 = Float32.((17 .<= LATS .<= 25) .&& (257 .<= LONS .<= 325)); not_sel2[isone.(not_sel2)] .= NaN
not_sel = Float32.((0 .<= LATS .<= 60) .&& (280 .<= LONS .<= 325)); not_sel[isone.(not_sel)] .= NaN
not_sel3 = Float32.((25 .<= LATS .<= 40) .&& (255 .<= LONS .<= 325)); not_sel3[isone.(not_sel3)] .= NaN
not_sel4 = Float32.((-60 .<= LATS .<= -30) .&& (0 .<= LONS .<= 140)); not_sel4[isone.(not_sel4)] .= NaN

PAC_msk = PAC_msk .+ not_sel .+ not_sel2 .+ not_sel3 .+ not_sel4
PAC_msk[isnan.(PAC_msk)] .= 0

fig, ax = plt.subplots()
ax.contourf(lon, lat, theta_OPT[1, 1, :, :] .* PAC_msk)
fig
#volumes for each cell. volume = 0 if not in area of interst

volumes = GH19_cell_volumes(depth, lon, lat)
mask_volume = similar(volumes)
[mask_volume[k, :, :] .= volumes[k, :, :] .* PAC_msk .* wet_mask[k, :, :] for k = 1:nz]

fig, ax = plt.subplots()
ax.contourf(lon, lat,  mask_volume[33, :, :])
fig


#weight the data
ΔTs = []; ps = []
WOCE_times = findall(1872 .< year .< 1876)[1]
Challenger_times = findall(1989 .< year .< 2017)[end]
println(year[WOCE_times] - year[Challenger_times] )
data_labels = ["EQ-0015", "OPT-0015"]
E,F = ECCOonPoseidon.trend_matrices(year[WOCE_times:Challenger_times])

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

